#!/usr/bin/env python3
"""Build a local reference library for one genus.

Turns the complete RefSeq genomes of a genus (downloaded with NCBI
`datasets`) into the versioned folder the app reads:

    <out>/
      manifest.json      version, thresholds, tool versions, counts, checksums
      library.sqlite     assemblies, replicons, genes, variant groups, names
      references/<acc>.fna.gz, .gff.gz
                         every assembly's genome and annotation as NCBI
                         serves them (references, pinned records, plasmids,
                         neighbourhoods are read from these)
      blast/genomes.*    nucleotide BLAST db of the kept replicons
      blast/groups_nt.*  nucleotide BLAST db of one gene per variant group
      blast/groups.*     protein BLAST db of one protein per variant group

Steps (see tools/library/README.md for the why):
  1. near-identical chromosomes are clustered with skani (greedy, best
     genome first) and one per cluster is kept;
  2. plasmids are clustered the same way, but must also align over most of
     their length, so a shared backbone with other cargo stays separate;
  3. the CDS and pseudogenes of every kept replicon go into a gene table;
  4. the proteins are clustered with MMseqs2 into variant groups (cadA1 and
     cadA2 apart, alleles of one together);
  5. the AMRFinderPlus and VFDB proteins are matched to the groups, so a
     group annotated "hypothetical protein" is still found by its name;
     groups with no name at all take the names of named groups they
     resemble (>= 50 % protein identity): the cadA of Tn5422 is annotated
     only "heavy metal translocating P-type ATPase";
  6. BLAST databases are made and every file is checksummed.

Only the Python standard library is used. External tools (on PATH):
skani, mmseqs, makeblastdb, blastp.

Usage:
  build_library.py --datasets DIR --catalogs DIR --genus Listeria \
      --out ~/straincompass/library/listeria/2026-10-02
"""

import argparse
import collections
import datetime
import gzip
import hashlib
import json
import os
import shutil
import sqlite3
import subprocess
import sys
import time
import urllib.request
from pathlib import Path

# Strains always kept as their own representative, whatever they cluster
# with: the references people name in papers and gene lists.
DEFAULT_KEEP = ["EGD-e", "10403S", "EGD", "F2365", "Clip11262", "Scott A", "LL195", "J1-220"]

SCHEMA = """
CREATE TABLE meta(key TEXT PRIMARY KEY, value TEXT);
CREATE TABLE assemblies(
  accession TEXT PRIMARY KEY, organism TEXT, strain TEXT,
  kept INTEGER, represented_by TEXT,
  assembly_name TEXT, paired_accession TEXT);
CREATE TABLE replicons(
  accession TEXT PRIMARY KEY, assembly TEXT, kind TEXT, name TEXT,
  length INTEGER, circular INTEGER, kept INTEGER,
  represented_by TEXT, represents INTEGER);
CREATE TABLE groups(
  id INTEGER PRIMARY KEY, label TEXT, product TEXT,
  rep_protein TEXT, rep_gene INTEGER, protein TEXT, nt TEXT,
  n_genes INTEGER, n_genomes INTEGER, n_plasmid INTEGER, n_chromosome INTEGER);
CREATE TABLE genes(
  id INTEGER PRIMARY KEY, replicon TEXT, start INTEGER, end INTEGER,
  strand TEXT, name TEXT, locus_tag TEXT, old_locus_tag TEXT,
  product TEXT, protein_id TEXT, pseudo INTEGER, group_id INTEGER,
  nt TEXT);
-- kind: gene (annotated name), locus (locus tag or old locus tag),
-- catalog (a curated entry matched by sequence), homolog (an unnamed group
-- resembling a named one; identity is the protein identity to it)
CREATE TABLE names(
  name TEXT COLLATE NOCASE, group_id INTEGER, kind TEXT, n INTEGER,
  identity REAL, via_group INTEGER);
CREATE TABLE catalog_hits(
  group_id INTEGER, source TEXT, symbol TEXT, product TEXT,
  identity REAL, coverage REAL);
"""

INDEXES = """
CREATE INDEX genes_group ON genes(group_id);
CREATE INDEX genes_place ON genes(replicon, start);
CREATE INDEX names_name ON names(name COLLATE NOCASE);
CREATE INDEX hits_group ON catalog_hits(group_id);
CREATE INDEX assemblies_paired ON assemblies(paired_accession);
"""


def log(msg):
    print(f"[{datetime.datetime.now():%H:%M:%S}] {msg}", file=sys.stderr, flush=True)


def run(cmd, **kw):
    log("$ " + " ".join(str(c) for c in cmd))
    subprocess.run([str(c) for c in cmd], check=True, **kw)


def tool_version(cmd):
    try:
        out = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
        text = (out.stdout or out.stderr).strip().splitlines()
        return text[0] if text else "?"
    except Exception:
        return "?"


def read_fasta(path):
    """{id: sequence} of a (possibly gzipped) FASTA file."""
    opener = gzip.open if str(path).endswith(".gz") else open
    seqs, cur, buf = {}, None, []
    with opener(path, "rt") as f:
        for line in f:
            if line.startswith(">"):
                if cur is not None:
                    seqs[cur] = "".join(buf)
                cur, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if cur is not None:
        seqs[cur] = "".join(buf)
    return seqs


def write_fasta(f, name, seq, desc=""):
    f.write(f">{name}{' ' + desc if desc else ''}\n")
    for i in range(0, len(seq), 80):
        f.write(seq[i : i + 80] + "\n")


COMPLEMENT = str.maketrans("ACGTRYKMBVDHNacgtrykmbvdhn", "TGCAYRMKVBHDNtgcayrmkvbhdn")


def revcomp(s):
    return s.translate(COMPLEMENT)[::-1]


def gff_attrs(field):
    out = {}
    for kv in field.split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            out[k] = (
                v.replace("%2C", ",").replace("%3B", ";").replace("%3D", "=").replace("%25", "%")
            )
    return out


# --- 1. assemblies -------------------------------------------------------


class Assembly:
    def __init__(self, rec, data_dir):
        self.accession = rec["accession"]
        org = rec.get("organism", {})
        self.organism = org.get("organismName", "")
        names = org.get("infraspecificNames", {})
        self.strain = names.get("strain") or names.get("isolate") or ""
        self.assembly_name = rec.get("assemblyInfo", {}).get("assemblyName", "")
        # the GenBank twin of a RefSeq assembly (GCA_ for GCF_)
        self.paired = rec.get("pairedAccession", "")
        info = rec.get("assemblyInfo", {})
        self.reference = info.get("refseqCategory") == "reference genome"
        cm = rec.get("checkmInfo") or {}
        self.quality = float(cm.get("completeness", 90)) - 5 * float(cm.get("contamination", 0))
        self.dir = data_dir / self.accession
        self.replicons = []  # filled by scan_replicons

    def files(self):
        fna = next(iter(sorted(self.dir.glob("*_genomic.fna"))), None)
        gff = self.dir / "genomic.gff"
        prot = self.dir / "protein.faa"
        return fna, gff, prot


def load_assemblies(data_dir):
    out = []
    with open(data_dir / "assembly_data_report.jsonl") as f:
        for line in f:
            a = Assembly(json.loads(line), data_dir)
            fna, gff, _ = a.files()
            if fna and gff.is_file():
                out.append(a)
            else:
                log(f"skipping {a.accession}: no genome or annotation file")
    return out


def priority(a, keep):
    named = any(k.lower() == a.strain.lower() for k in keep)
    return (not named, not a.reference, -a.quality, a.accession)


def scan_replicons(a):
    """The assembly's replicons from the GFF region lines."""
    _, gff, _ = a.files()
    reps = []
    with open(gff) as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9 or p[2] != "region" or p[3] != "1":
                continue
            at = gff_attrs(p[8])
            kind = "plasmid" if at.get("genome") == "plasmid" else "chromosome"
            name = at.get("plasmid-name", "") if kind == "plasmid" else ""
            reps.append(
                {
                    "accession": p[0],
                    "assembly": a.accession,
                    "kind": kind,
                    "name": name,
                    "length": int(p[4]),
                    "circular": at.get("Is_circular") == "true",
                }
            )
    a.replicons = reps


# --- 2. clustering with skani --------------------------------------------


def skani_edges(files, out_tsv, threads, small):
    """{(a, b): (ani, af_a, af_b)} over file stems, both directions."""
    lst = out_tsv.with_suffix(".list")
    lst.write_text("".join(f"{p}\n" for p in files))
    cmd = ["skani", "triangle", "-l", lst, "-t", threads, "--sparse", "-o", out_tsv]
    if small:
        # plasmids: --small-genomes' sensitivity, without its
        # --faster-small, which may skip comparisons
        cmd += ["-c", "30", "-m", "200"]
    run(cmd)
    edges = {}
    with open(out_tsv) as f:
        header = f.readline().rstrip("\n").split("\t")
        ix = {h: i for i, h in enumerate(header)}
        for line in f:
            p = line.rstrip("\n").split("\t")
            r = Path(p[ix["Ref_file"]]).stem
            q = Path(p[ix["Query_file"]]).stem
            ani = float(p[ix["ANI"]])
            af_r = float(p[ix["Align_fraction_ref"]])
            af_q = float(p[ix["Align_fraction_query"]])
            edges[(r, q)] = (ani, af_r, af_q)
            edges[(q, r)] = (ani, af_q, af_r)
    return edges


def greedy(order, edges, min_ani, min_af):
    """{member: representative}, walking `order` best first."""
    reps, rep_of = [], {}
    for x in order:
        hit = None
        for r in reps:
            e = edges.get((x, r))
            if e and e[0] >= min_ani and e[1] >= min_af and e[2] >= min_af:
                hit = r
                break
        if hit is None:
            reps.append(x)
            rep_of[x] = x
        else:
            rep_of[x] = hit
    return rep_of


# --- 3. genes ------------------------------------------------------------


def parse_genes(a, kept, genome):
    """CDS and pseudogenes on the assembly's kept replicons.

    Multi-part features (one ID on several lines) are joined; features
    across the origin of a circular replicon have an end beyond its length.
    """
    _, gff, _ = a.files()
    old_tag = {}
    feats = collections.OrderedDict()
    with open(gff) as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9 or p[0] not in kept:
                continue
            at = gff_attrs(p[8])
            if p[2] in ("gene", "pseudogene"):
                if "old_locus_tag" in at and "locus_tag" in at:
                    old_tag[at["locus_tag"]] = at["old_locus_tag"].split(",")[0]
                continue
            if p[2] != "CDS":
                continue
            fid = at.get("ID", f"{p[0]}:{p[3]}")
            part = (int(p[3]), int(p[4]))
            if fid in feats:
                feats[fid]["parts"].append(part)
            else:
                feats[fid] = {"seqid": p[0], "strand": p[6], "parts": [part], "at": at}
    genes = []
    for ft in feats.values():
        at, seqid = ft["at"], ft["seqid"]
        s = genome[seqid]
        parts = sorted(ft["parts"])
        chunks = []
        for lo, hi in parts:
            if hi > len(s):  # across the origin
                chunks.append(s[lo - 1 :] + s[: hi - len(s)])
            else:
                chunks.append(s[lo - 1 : hi])
        nt = "".join(chunks).upper()
        if ft["strand"] == "-":
            nt = revcomp(nt)
        pseudo = at.get("pseudo") == "true"
        prot = at.get("protein_id", "")
        inferred = ""
        if pseudo:
            # "similar to AA sequence:RefSeq:WP_010958663.1": the protein
            # this broken copy was once
            inf = at.get("inference", "")
            if "RefSeq:WP_" in inf:
                inferred = inf.split("RefSeq:")[-1].split(",")[0]
        tag = at.get("locus_tag", "")
        genes.append(
            {
                "replicon": seqid,
                "start": parts[0][0],
                "end": parts[-1][1],
                "strand": ft["strand"],
                "name": at.get("gene", ""),
                "locus_tag": tag,
                "old_locus_tag": old_tag.get(tag, ""),
                "product": at.get("product", ""),
                "protein_id": prot,
                "pseudo": pseudo,
                "inferred": inferred,
                "nt": nt,
            }
        )
    return genes


# --- 5. catalogs ---------------------------------------------------------


def catalog_records(cat_dir):
    """[(source, symbol, product, id, seq)] of the curated protein sets."""
    out = []
    amr = cat_dir / "AMRProt.fa"
    if amr.is_file():
        for h, seq in read_fasta_with_headers(amr):
            f = h.split("|")
            if len(f) >= 10:
                out.append(("AMRFinderPlus", f[3], f[-1].replace("_", " "), f[0], seq))
    for vf in sorted(cat_dir.glob("VFDB_setA_pro.fas*")):
        for h, seq in read_fasta_with_headers(vf):
            # VFG037176(gb|WP_001081735) (plc1) phospholipase C [...] [organism]
            vid = h.split()[0]
            rest = h[len(vid) :].strip()
            if not rest.startswith("("):
                continue
            gene, _, after = rest[1:].partition(")")
            product = after.split("[")[0].strip()
            out.append(("VFDB", gene.strip(), product, vid, seq))
        break
    return out


def read_fasta_with_headers(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    h, buf = None, []
    with opener(path, "rt", errors="replace") as f:
        for line in f:
            if line.startswith(">"):
                if h is not None:
                    yield h, "".join(buf)
                h, buf = line[1:].strip(), []
            else:
                buf.append(line.strip())
    if h is not None:
        yield h, "".join(buf)


def seed_records(seeds_tsv, cache_dir):
    """[(source, symbol, product, id, seq)] of the seed list: the named
    proteins of GenBank records, fetched once and kept in `cache_dir`."""
    if not seeds_tsv.is_file():
        return []
    wanted = collections.defaultdict(set)
    for line in seeds_tsv.read_text().splitlines():
        if line.strip() and not line.startswith("#"):
            gene, acc = line.split("\t")[:2]
            wanted[acc.strip()].add(gene.strip())
    out = []
    for acc, genes in sorted(wanted.items()):
        cached = cache_dir / f"{acc}.faa"
        if not cached.is_file():
            url = (
                "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
                f"?db=nuccore&id={acc}&rettype=fasta_cds_aa&retmode=text"
            )
            log(f"fetching seed record {acc}")
            with urllib.request.urlopen(url, timeout=120) as r:
                cached.write_bytes(r.read())
            time.sleep(0.4)
        for h, seq in read_fasta_with_headers(cached):
            # >lcl|L28104.1_prot_AAA25275.1_1 [gene=cadA] [protein=ATPase] ...
            tags = dict(t.split("=", 1) for t in h.split(" [")[1:] if "=" in t)
            tags = {k: v.rstrip("]") for k, v in tags.items()}
            gene = tags.get("gene", "")
            if gene in genes:
                out.append((f"GenBank {acc}", gene, tags.get("protein", ""), h.split()[0], seq))
        missing = genes - {o[1] for o in out if o[0] == f"GenBank {acc}"}
        if missing:
            log(f"seed record {acc} names no {', '.join(sorted(missing))}")
    return out


def catalog_names(symbol):
    """Names a catalog symbol answers to, as catalog.rs reads them:
    cadA_Lm -> cadA, SMR_efflux_bcrB -> bcrB, hbp1/svpA -> both."""
    names = set()
    for alias in symbol.split("/"):
        alias = alias.strip()
        if not alias:
            continue
        names.add(alias)
        head, _, tail = alias.partition("_")
        if tail and len(tail) <= 3:
            names.add(head)
        names.add(alias.rsplit("_", 1)[-1])
    return {n for n in names if len(n) >= 3}


# --- main ----------------------------------------------------------------


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--datasets", required=True, type=Path, help="unzipped `datasets` download (the folder holding ncbi_dataset/)")
    ap.add_argument("--catalogs", required=True, type=Path, help="folder with AMRProt.fa, version.txt and VFDB_setA_pro.fas.gz")
    ap.add_argument("--seeds", type=Path, default=Path(__file__).with_name("seeds.tsv"), help="GenBank records naming genes RefSeq leaves unnamed")
    ap.add_argument("--genus", required=True)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--work", type=Path, help="scratch folder (default: <out>.work)")
    ap.add_argument("--threads", type=int, default=max(1, (os.cpu_count() or 4) // 2), help="default: half the cores")
    ap.add_argument("--keep", default=",".join(DEFAULT_KEEP), help="strains always kept as their own representative")
    ap.add_argument("--chrom-ani", type=float, default=99.9)
    ap.add_argument("--chrom-af", type=float, default=90.0)
    ap.add_argument("--plasmid-ani", type=float, default=99.9)
    ap.add_argument("--plasmid-af", type=float, default=95.0)
    ap.add_argument("--group-id", type=float, default=0.9)
    ap.add_argument("--group-cov", type=float, default=0.8)
    ap.add_argument("--catalog-id", type=float, default=80.0)
    ap.add_argument("--catalog-cov", type=float, default=80.0)
    ap.add_argument("--homolog-id", type=float, default=50.0)
    ap.add_argument("--homolog-cov", type=float, default=80.0)
    args = ap.parse_args()
    # the lowest CPU priority, inherited by skani, mmseqs and BLAST: the
    # build takes minutes either way, and the machine stays usable
    os.nice(19)

    for t in ("skani", "mmseqs", "makeblastdb", "blastp"):
        if not shutil.which(t):
            sys.exit(f"{t} is not on PATH")
    out = args.out.expanduser()
    if out.exists():
        sys.exit(f"{out} exists; a library version is never overwritten")
    work = (args.work or out.with_name(out.name + ".work")).expanduser()
    work.mkdir(parents=True, exist_ok=True)
    build = out.with_name(out.name + ".partial")
    if build.exists():
        shutil.rmtree(build)
    (build / "blast").mkdir(parents=True)
    data_dir = args.datasets.expanduser() / "ncbi_dataset" / "data"
    keep = [k.strip() for k in args.keep.split(",") if k.strip()]

    # 1. assemblies and their replicons
    assemblies = load_assemblies(data_dir)
    assemblies.sort(key=lambda a: priority(a, keep))
    log(f"{len(assemblies)} assemblies")
    chrom_dir, plas_dir = work / "chromosomes", work / "plasmids"
    chrom_dir.mkdir(exist_ok=True)
    plas_dir.mkdir(exist_ok=True)
    replicons = {}
    for a in assemblies:
        scan_replicons(a)
        for r in a.replicons:
            replicons[r["accession"]] = r
        cf = chrom_dir / f"{a.accession}.fna"
        if cf.exists() and all((plas_dir / f"{r['accession']}.fna").exists() for r in a.replicons if r["kind"] == "plasmid"):
            continue
        fna, _, _ = a.files()
        genome = read_fasta(fna)
        with open(cf, "w") as f:
            for r in a.replicons:
                if r["kind"] == "chromosome" and r["accession"] in genome:
                    write_fasta(f, r["accession"], genome[r["accession"]])
        for r in a.replicons:
            if r["kind"] == "plasmid" and r["accession"] in genome:
                with open(plas_dir / f"{r['accession']}.fna", "w") as f:
                    write_fasta(f, r["accession"], genome[r["accession"]])
    n_plasmids = sum(1 for r in replicons.values() if r["kind"] == "plasmid")
    log(f"{len(replicons)} replicons, {n_plasmids} of them plasmids")

    # 2. near-duplicates
    order = [a.accession for a in assemblies]
    edges = skani_edges([chrom_dir / f"{x}.fna" for x in order], work / "chromosomes.tsv", args.threads, False)
    keep_set = {a.accession for a in assemblies if any(k.lower() == a.strain.lower() for k in keep)}
    rep_of = {}
    reps = []
    for x in order:
        hit = None
        if x not in keep_set:
            for r in reps:
                e = edges.get((x, r))
                if e and e[0] >= args.chrom_ani and e[1] >= args.chrom_af and e[2] >= args.chrom_af:
                    hit = r
                    break
        if hit is None:
            reps.append(x)
            rep_of[x] = x
        else:
            rep_of[x] = hit
    log(f"chromosomes: {len(reps)} representatives for {len(order)} assemblies")

    a_rank = {a.accession: i for i, a in enumerate(assemblies)}
    plasmids = sorted(
        (r for r in replicons.values() if r["kind"] == "plasmid"),
        key=lambda r: (a_rank[r["assembly"]], r["accession"]),
    )
    p_order = [r["accession"] for r in plasmids if (plas_dir / f"{r['accession']}.fna").exists()]
    p_edges = skani_edges([plas_dir / f"{x}.fna" for x in p_order], work / "plasmids.tsv", args.threads, True) if p_order else {}
    p_rep_of = greedy(p_order, p_edges, args.plasmid_ani, args.plasmid_af)
    log(f"plasmids: {len(set(p_rep_of.values()))} distinct of {len(p_order)}")

    for r in replicons.values():
        if r["kind"] == "chromosome":
            r["represented_by"] = r["accession"] if rep_of[r["assembly"]] == r["assembly"] else ""
            r["kept"] = rep_of[r["assembly"]] == r["assembly"]
        else:
            rep = p_rep_of.get(r["accession"], r["accession"])
            r["represented_by"] = rep
            r["kept"] = rep == r["accession"]
    # a dropped chromosome is represented by the kept one of its cluster
    first_chrom = {}
    for r in replicons.values():
        if r["kind"] == "chromosome" and r["kept"]:
            first_chrom.setdefault(r["assembly"], r["accession"])
    for r in replicons.values():
        if r["kind"] == "chromosome" and not r["kept"]:
            r["represented_by"] = first_chrom.get(rep_of[r["assembly"]], "")
    represents = collections.Counter(r["represented_by"] for r in replicons.values())
    for r in replicons.values():
        r["represents"] = represents[r["accession"]] if r["kept"] else 0
    kept = {acc for acc, r in replicons.items() if r["kept"]}

    # 3. genes of the kept replicons, and the replicons themselves
    db_path = build / "library.sqlite"
    db = sqlite3.connect(db_path)
    db.executescript(SCHEMA)
    proteins = {}  # protein id -> sequence
    first_gene = {}  # protein id -> gene id of its first (best genome) copy
    gene_rows = []
    plain = open(work / "genomes.fna", "w")
    refs = build / "references"
    refs.mkdir()
    for a in assemblies:
        # every assembly, kept or not: a reference asked for by accession
        # is the genome asked for, not its representative
        fna, gff, prot = a.files()
        for src, dst in ((fna, refs / f"{a.accession}.fna.gz"), (gff, refs / f"{a.accession}.gff.gz")):
            with open(src, "rb") as fi, gzip.open(dst, "wb", compresslevel=6) as fo:
                shutil.copyfileobj(fi, fo)
        mine = {r["accession"] for r in a.replicons if r["accession"] in kept}
        if not mine:
            continue
        genome = read_fasta(fna)
        for r in a.replicons:
            if r["accession"] not in mine:
                continue
            desc = f"{r['kind']}|{r['name']}|{a.organism}"
            write_fasta(plain, r["accession"], genome[r["accession"]], desc)
        prots = read_fasta(prot) if prot.is_file() else {}
        for g in parse_genes(a, mine, genome):
            gid = len(gene_rows) + 1
            g["id"] = gid
            pid = g["protein_id"]
            if pid and pid in prots:
                proteins.setdefault(pid, prots[pid])
                first_gene.setdefault(pid, gid)
            gene_rows.append(g)
    plain.close()
    log(f"{len(gene_rows)} genes, {len(proteins)} distinct proteins")

    # 4. variant groups
    prot_fa = work / "proteins.faa"
    with open(prot_fa, "w") as f:
        for pid, seq in proteins.items():
            write_fasta(f, pid, seq)
    clu = work / "clu"
    run(
        [
            "mmseqs", "easy-cluster", prot_fa, clu, work / "mmseqs_tmp",
            "--min-seq-id", args.group_id, "-c", args.group_cov, "--cov-mode", "0",
            "--threads", args.threads, "-v", "1",
        ],
        stdout=subprocess.DEVNULL,
    )
    group_of_prot, rep_prots = {}, []
    with open(f"{clu}_cluster.tsv") as f:
        for line in f:
            rep, mem = line.split()
            if rep not in group_of_prot:
                group_of_prot[rep] = len(rep_prots) + 1
                rep_prots.append(rep)
            group_of_prot[mem] = group_of_prot[rep]
    for g in gene_rows:
        pid = g["protein_id"] or g["inferred"]
        g["group_id"] = group_of_prot.get(pid)
    log(f"{len(rep_prots)} variant groups")

    by_group = collections.defaultdict(list)
    for g in gene_rows:
        if g["group_id"] and not g["pseudo"]:
            by_group[g["group_id"]].append(g)
    groups = []
    for gid, rep in enumerate(rep_prots, start=1):
        members = by_group.get(gid, [])
        names = collections.Counter(m["name"] for m in members if m["name"])
        products = collections.Counter(m["product"] for m in members if m["product"])
        rep_gene = gene_rows[first_gene[rep] - 1]
        weight = collections.Counter()
        for rc in {m["replicon"] for m in members}:
            weight[replicons[rc]["kind"]] += replicons[rc]["represents"]
        groups.append(
            {
                "id": gid,
                "label": names.most_common(1)[0][0] if names else "",
                "product": products.most_common(1)[0][0] if products else "",
                "rep_protein": rep,
                "rep_gene": rep_gene["id"],
                "protein": proteins[rep],
                "nt": rep_gene["nt"],
                "n_genes": len(members),
                "n_genomes": weight["plasmid"] + weight["chromosome"],
                "n_plasmid": weight["plasmid"],
                "n_chromosome": weight["chromosome"],
                "names": names,
                "tags": collections.Counter(
                    t for m in members for t in (m["locus_tag"], m["old_locus_tag"]) if t
                ),
            }
        )

    groups_faa = work / "groups.faa"
    groups_fna = work / "groups.fna"
    with open(groups_faa, "w") as fa, open(groups_fna, "w") as fn:
        for g in groups:
            desc = f"{g['label'] or '-'} {g['product']}"
            write_fasta(fa, f"g{g['id']}", g["protein"], desc)
            write_fasta(fn, f"g{g['id']}", g["nt"], desc)

    # 5. curated catalogs onto the groups
    cats = catalog_records(args.catalogs.expanduser())
    cats += seed_records(args.seeds.expanduser(), args.catalogs.expanduser())
    hits = []
    if cats:
        cat_fa = work / "catalogs.faa"
        meta = {}
        with open(cat_fa, "w") as f:
            for i, (src, sym, prod, cid, seq) in enumerate(cats):
                write_fasta(f, f"c{i}", seq)
                meta[f"c{i}"] = (src, sym, prod)
        run(["makeblastdb", "-in", groups_faa, "-dbtype", "prot", "-out", work / "groups_tmp"], stdout=subprocess.DEVNULL)
        res = work / "catalogs.tsv"
        run(
            [
                "blastp", "-query", cat_fa, "-db", work / "groups_tmp", "-evalue", "1e-10",
                "-max_target_seqs", "20", "-num_threads", args.threads,
                "-outfmt", "6 qseqid sseqid pident length qlen slen", "-out", res,
            ]
        )
        best = {}
        with open(res) as f:
            for line in f:
                q, s, pid, length, qlen, slen = line.split()
                pid, length, qlen, slen = float(pid), int(length), int(qlen), int(slen)
                cov = 100.0 * length / max(qlen, slen)
                if pid < args.catalog_id or cov < args.catalog_cov:
                    continue
                key = (int(s[1:]), meta[q][0], meta[q][1])
                if key not in best or best[key][0] < pid:
                    best[key] = (pid, cov, meta[q][2])
        for (gid, src, sym), (pid, cov, prod) in best.items():
            hits.append((gid, src, sym, prod, round(pid, 1), round(cov, 1)))
        log(f"{len(hits)} catalog matches on {len({h[0] for h in hits})} groups")

    # unnamed groups take the names of the named groups they resemble
    named = {g["id"]: set(g["names"]) for g in groups if g["names"]}
    for gid, _, sym, _, _, _ in hits:
        named.setdefault(gid, set()).update(catalog_names(sym))
    homologs = []
    res = work / "homologs.tsv"
    run(
        [
            "mmseqs", "easy-search", groups_faa, groups_faa, res, work / "mmseqs_tmp",
            "--min-seq-id", args.homolog_id / 100, "-c", args.homolog_cov / 100, "--cov-mode", "0",
            "-e", "1e-10", "--threads", args.threads, "-v", "1",
            "--format-output", "query,target,fident,qcov,tcov",
        ],
        stdout=subprocess.DEVNULL,
    )
    best_named = {}
    with open(res) as f:
        for line in f:
            q, t, fid, qcov, tcov = line.split()
            qi, ti = int(q[1:]), int(t[1:])
            if qi == ti or qi in named or ti not in named:
                continue
            ident = 100 * float(fid)
            if qi not in best_named or best_named[qi][0] < ident:
                best_named[qi] = (ident, ti)
    for qi, (ident, ti) in best_named.items():
        homologs += [(n, qi, "homolog", 1, round(ident, 1), ti) for n in sorted(named[ti])]
    log(f"{len(best_named)} unnamed groups named after a resembling group")

    # tables
    db.executemany(
        "INSERT INTO assemblies VALUES (?,?,?,?,?,?,?)",
        [
            (a.accession, a.organism, a.strain, rep_of[a.accession] == a.accession, rep_of[a.accession], a.assembly_name, a.paired)
            for a in assemblies
        ],
    )
    db.executemany(
        "INSERT INTO replicons VALUES (?,?,?,?,?,?,?,?,?)",
        [
            (r["accession"], r["assembly"], r["kind"], r["name"], r["length"], r["circular"], r["kept"], r["represented_by"], r["represents"])
            for r in replicons.values()
        ],
    )
    db.executemany(
        "INSERT INTO genes VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)",
        (
            (
                g["id"], g["replicon"], g["start"], g["end"], g["strand"], g["name"], g["locus_tag"],
                g["old_locus_tag"], g["product"], g["protein_id"] or g["inferred"], g["pseudo"],
                g["group_id"], g["nt"] if g["pseudo"] else None,
            )
            for g in gene_rows
        ),
    )
    db.executemany(
        "INSERT INTO groups VALUES (?,?,?,?,?,?,?,?,?,?,?)",
        [
            (g["id"], g["label"], g["product"], g["rep_protein"], g["rep_gene"], g["protein"], g["nt"],
             g["n_genes"], g["n_genomes"], g["n_plasmid"], g["n_chromosome"])
            for g in groups
        ],
    )
    name_rows = []
    for g in groups:
        name_rows += [(n, g["id"], "gene", c, 100.0, None) for n, c in g["names"].items()]
        name_rows += [(t, g["id"], "locus", c, 100.0, None) for t, c in g["tags"].items()]
    for gid, src, sym, _, ident, _ in hits:
        name_rows += [(n, gid, "catalog", 1, ident, None) for n in catalog_names(sym)]
    name_rows += homologs
    db.executemany("INSERT INTO names VALUES (?,?,?,?,?,?)", name_rows)
    db.executemany("INSERT INTO catalog_hits VALUES (?,?,?,?,?,?)", hits)
    db.executescript(INDEXES)

    version = out.name
    amr_version = (args.catalogs / "version.txt").read_text().strip() if (args.catalogs / "version.txt").is_file() else ""
    meta = {
        "genus": args.genus,
        "version": version,
        "built": datetime.date.today().isoformat(),
        "amrfinder_db": amr_version,
    }
    db.executemany("INSERT INTO meta VALUES (?,?)", list(meta.items()))
    db.commit()
    db.execute("VACUUM")
    db.close()

    # 6. BLAST databases
    title = f"{args.genus} reference library {version}"
    run(["makeblastdb", "-in", work / "genomes.fna", "-dbtype", "nucl", "-parse_seqids", "-out", build / "blast" / "genomes", "-title", title], stdout=subprocess.DEVNULL)
    run(["makeblastdb", "-in", groups_fna, "-dbtype", "nucl", "-parse_seqids", "-out", build / "blast" / "groups_nt", "-title", title], stdout=subprocess.DEVNULL)
    run(["makeblastdb", "-in", groups_faa, "-dbtype", "prot", "-parse_seqids", "-out", build / "blast" / "groups", "-title", title], stdout=subprocess.DEVNULL)

    files = {}
    for p in sorted(build.rglob("*")):
        if p.is_file():
            files[str(p.relative_to(build))] = {"sha256": sha256(p), "bytes": p.stat().st_size}
    n_kept_plasmids = sum(1 for r in replicons.values() if r["kept"] and r["kind"] == "plasmid")
    manifest = {
        "format": 1,
        **meta,
        "source": "NCBI RefSeq, complete genomes",
        "counts": {
            "assemblies": len(assemblies),
            "representative_assemblies": len(reps),
            "plasmids": n_plasmids,
            "distinct_plasmids": n_kept_plasmids,
            "replicons_kept": len(kept),
            "genes": len(gene_rows),
            "pseudogenes": sum(1 for g in gene_rows if g["pseudo"]),
            "groups": len(groups),
            "catalog_matches": len(hits),
            "groups_named_by_resemblance": len(best_named),
        },
        "thresholds": {
            "chromosome_ani": args.chrom_ani,
            "chromosome_af": args.chrom_af,
            "plasmid_ani": args.plasmid_ani,
            "plasmid_af": args.plasmid_af,
            "group_identity": args.group_id,
            "group_coverage": args.group_cov,
            "catalog_identity": args.catalog_id,
            "catalog_coverage": args.catalog_cov,
            "homolog_identity": args.homolog_id,
            "homolog_coverage": args.homolog_cov,
        },
        "always_kept": keep,
        "tools": {
            "skani": tool_version(["skani", "--version"]),
            "mmseqs": tool_version(["mmseqs", "version"]),
            "makeblastdb": tool_version(["makeblastdb", "-version"]),
        },
        "files": files,
    }
    (build / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    build.rename(out)
    log(f"done: {out}")
    print(json.dumps(manifest["counts"], indent=2))


if __name__ == "__main__":
    main()
