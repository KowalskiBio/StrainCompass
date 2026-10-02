# Reference libraries

A reference library is a local, versioned database of one genus' complete
genomes, with every gene sorted into variant groups. The app asks the
library first and NCBI only for what the library cannot answer:

| Feature | Library first | NCBI still used for |
|---|---|---|
| Reference by accession | any of the 835 genomes, by GCF_ or GCA_ | genomes outside the library |
| Gene panel by name | every variant of the name in the genus (Listeria has several cadA genes; the most common, Tn5422's, is not named cadA in RefSeq's annotation) | names neither the library nor the curated catalogs know |
| Gene pinned to a record ("qacH (HF565366.1)") | the library's copy of the record (GenBank accessions through their RefSeq twin) | records outside the library, or genes its annotation does not name |
| Genes beside a panel gene's source record | the record's stored annotation | records outside the library |
| Where a gene usually occurs | the genus search, in a second | "search all bacteria" |
| Plasmid for the whole-element comparison | the library's copy | plasmids outside the library |
| Naming novel genes in gained regions | the library's protein families, then SwissProt | the on-demand nr search per region |

Without a library everything works as before, through the curated catalogs
and NCBI.

## Layout

```
<library root>/<genus, lower case>/
  current -> 2026-10-02         the version the app reads (a symlink)
  2026-10-02/
    manifest.json               version, sources, thresholds, tool versions,
                                counts, size + sha256 of every file
    library.sqlite              assemblies, replicons, genes, groups, names
    references/<GCF>.fna.gz     every assembly's genome and annotation, as
    references/<GCF>.gff.gz     NCBI serves them (kept and merged alike)
    blast/genomes.*             nucleotide BLAST db of the kept replicons
    blast/groups_nt.*           one gene per variant group (DNA)
    blast/groups.*              one protein per variant group
```

The root is `STRAINCOMPASS_LIBRARY_DIR`, else `library/` beside the data
directory (on the VM: `~/straincompass/library`). Never put it inside the
data directory: every deploy backup copies all of it.

The app refuses a library whose files do not match the manifest sizes (a
half-copied folder) and falls back to NCBI; Settings shows why.

## How it is built

`build_library.py` (Python standard library only) runs these steps:

1. **Assemblies.** The complete RefSeq genomes of the genus, from NCBI
   `datasets`. Each sequence is a chromosome or a plasmid (from the GFF
   region lines), and genomes are ranked: named reference strains first
   (EGD-e, 10403S, ...), then NCBI's reference genome, then CheckM quality.
2. **Near-duplicate chromosomes** are found with skani (MinHash sketches,
   then ANI on the shared k-mers) and clustered greedily down the ranking: a
   genome joins the first kept one it matches at >= 99.9 % ANI over >= 90 %
   of both, else it is kept. Each kept replicon records how many it stands
   for, so counts stay counts of genomes.
3. **Plasmids** are clustered the same way, but must also align over >= 95 %
   of both: plasmids sharing a backbone with other cargo stay apart.
4. **Genes.** Every CDS and pseudogene of the kept replicons goes into the
   gene table, with its old locus tag (lmo0444). The locus tags of the
   merged genomes are indexed too (lm4b_02324 of CLIP 80459), each tied to
   its gene's group through the protein accession RefSeq gives identical
   copies. A pseudogene is tied to
   the protein it once was through its `inference` (groundwork for telling
   relics apart).
5. **Variant groups.** The proteins are clustered with MMseqs2 at 90 %
   identity over 80 % coverage: cadA of Tn5422 and of pLI100 (69 % alike)
   apart, alleles of one gene together.
6. **Names.** A group answers to its annotated names and locus tags; to the
   AMRFinderPlus and VFDB entries matching its protein (>= 80 % over 80 %);
   to the genes of the seed records (`seeds.tsv`, GenBank records naming
   genes RefSeq leaves unnamed); and, when it has no name at all, to the
   names of a named group it resembles (>= 50 % protein identity). The app
   offers a resembling group as a variant only from 70 %: the pLI100 cadA
   (76 % to Tn5422's) yes, internalins resembling inlA (50-68 %) no.
7. **References.** Every assembly's genome and annotation is stored gzipped,
   kept or merged: a genome asked for by accession is that genome, not its
   representative.
8. **Indexes.** `makeblastdb` packs the replicons and group genes into BLAST
   databases; every file is checksummed into the manifest.

## Building a new version

On a workstation (not the VM). The build runs at the lowest CPU priority
and by default on half the cores; Listeria (835 genomes) takes about
fifteen minutes and makes a ~1.4 GB folder (about 1 GB of it the stored
genomes).

```sh
conda env create -f tools/library/environment.yml     # once
conda activate sclib

mkdir -p ~/sc-library-build && cd ~/sc-library-build
datasets download genome taxon Listeria --assembly-source refseq \
  --assembly-level complete --include genome,gff3,cds,protein,seq-report \
  --dehydrated --filename listeria.zip
unzip -q listeria.zip -d listeria && datasets rehydrate --directory listeria

mkdir -p catalogs && cd catalogs
curl -O https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/latest/AMRProt.fa
curl -O https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/latest/version.txt
curl -O http://www.mgc.ac.cn/VFs/Down/VFDB_setA_pro.fas.gz
cd ..

python3 <repo>/tools/library/build_library.py --datasets listeria \
  --catalogs catalogs --genus Listeria --out out/listeria/$(date +%F)
```

A version folder is never overwritten; the build writes to `<out>.partial`
and renames it when complete.

## Installing on the server

Copy the version folder next to the existing ones, then switch the link:

```sh
# files unchanged since the installed version are hard-linked, not copied
rsync -a --link-dest=../2026-10-02 out/listeria/2026-10-03 \
  proxmox1:straincompass/library/listeria/
ssh proxmox1 'ln -sfn 2026-10-03 straincompass/library/listeria/current'
```

The stored genomes are written without a timestamp, so they are byte-for-
byte the same between builds of the same genomes and link instead of copy
(builds up to 2026-10-02.2 wrote timestamped gzips and share nothing).

The app reads the library per request, so no restart is needed. Keep the
previous version until the new one has been checked; switching back is the
same `ln -sfn`.

## One name, several genes

A gene name is not an identity. EGD-e calls both the virulence regulator
PrfA (lmo0200) and peptide chain release factor 1 (lmo2543) "prfA"; other
genomes call several shikimate dehydrogenases "aroE". The panel therefore:

- takes, among the reference's genes of a name, the one a curated entry of
  that name is (VFDB's prfA is lmo0200), else the first in the genome, and
  says which it took and how to ask for the other (its locus tag);
- adds as variants only library groups related to the gene taken (a blastx
  match, or the same product); unrelated groups sharing the name are listed,
  with a pinned entry to search them under a name of their own;
- leaves out groups the reference genome itself also carries: paralogs of
  the gene, not other versions of it.
