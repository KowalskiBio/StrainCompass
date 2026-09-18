export type Call = "PRESENT" | "PARTIAL" | "ABSENT";

export interface Project {
  id: number;
  name: string;
  organism: string;
  created_at: string;
  n_runs: number;
  n_queries: number;
  has_reference: boolean;
  has_panel: boolean;
  usage_bytes: number;
}

export interface ProjectFile {
  id: number;
  role: "reference_fasta" | "reference_gff" | "query" | "panel";
  display_name: string;
  size: number;
  created_at: string;
}

export interface RunQuery {
  file_id: number;
  name: string;
}

export interface Run {
  id: number;
  project_id: number;
  status: "queued" | "running" | "succeeded" | "failed";
  step: string | null;
  error: string | null;
  created_at: string;
  started_at: string | null;
  finished_at: string | null;
  queries: RunQuery[];
  has_panel: boolean;
  /** Runs computed before gained regions existed do not carry them. */
  has_gained: boolean;
}

export interface RunWithLogs {
  run: Run;
  logs: string[];
}

export interface GeneCoverageRow {
  locus_tag: string;
  symbol: string;
  protein_id: string;
  biotype: string;
  seqid: string;
  start: number;
  end: number;
  length: number;
  cov_bp: number;
  cov_pct: number;
  call: Call;
  best_identity: number;
  mismatches: number;
  indels: number;
}

export interface GapRow {
  seqid: string;
  start: number;
  end: number;
  length: number;
  n_genes: number;
  genes: string[];
}

/** Where a gained region sits relative to the reference. */
export type GainedAnchor = "between" | "flank" | "unanchored";

export interface GainedOrf {
  /** 1-based inclusive, in query contig coordinates. */
  start: number;
  end: number;
  strand: number;
  /** The ORF runs off an edge of the region, so it is probably truncated. */
  partial: boolean;
  confidence: number;
}

/**
 * One stretch of a query genome with no alignment to the reference.
 *
 * "No alignment" is weaker than "not in the reference" - see the engine's
 * gained module - so nothing here should be labelled as simply new.
 */
export interface GainedRow {
  qry_seqid: string;
  start: number;
  end: number;
  length: number;
  gc_pct: number;
  at_contig_end: boolean;
  anchor: GainedAnchor;
  anchor_seqid: string;
  anchor_start: number;
  anchor_end: number;
  flanks_disagree: boolean;
  left_gene: string;
  right_gene: string;
  /** null means gene prediction did not run, which is not a count of zero. */
  n_orfs: number | null;
  n_orfs_complete: number | null;
  orfs: GainedOrf[];
}

/** One search hit of a gained region against the reference: blastn for
 * the nucleotide tiers, tblastx for the translated tier. */
export interface GainedBlastHit {
  ref_seqid: string;
  /** 1-based inclusive, ordered low..high regardless of the hit's strand. */
  ref_start: number;
  ref_end: number;
  /** Nucleotide identity for blastn, amino-acid identity for tblastx. */
  identity: number;
  /** Bases for blastn, aligned residues for tblastx. */
  length: number;
  /** Interval on the region itself, not the whole query contig. */
  qry_start: number;
  qry_end: number;
  evalue: number;
  bitscore: number;
  /** Reference genes the hit overlaps, display-formatted
   * ("LM4B_RS13070 (rlmN)"). A hit means little until it says what it
   * hits. */
  genes: string[];
}

/**
 * The reference back-check of one gained region: the region's sequence
 * searched against the reference genome three ways - a nucleotide
 * search split into a strong and a weak tier, a translated search
 * (tblastx) that runs only when the strong tier is empty, and the
 * longest exact match, cutoff-free. Empty everything is the closest
 * available evidence of absence.
 */
export interface GainedVerify {
  qry_seqid: string;
  start: number;
  end: number;
  /** Nucleotide hits with E <= 1e-5, best first. */
  hits: GainedBlastHit[];
  /** Nucleotide hits with 1e-5 < E <= 10: below the verdicts'
   * threshold, shown rather than silently dropped. Usually noise. */
  weak_hits: GainedBlastHit[];
  /** Translated (tblastx) hits. null means the search did not run,
   * which is not the same as having looked and found none. */
  tx_hits: GainedBlastHit[] | null;
  /** Plain-language reason the translated search did not run. */
  tx_note: string | null;
  /** Longest run of bases the region shares, exactly, with any position
   * of the reference. Chance alone gives ~log4(region * reference). */
  longest_exact_bp: number;
  /** Where that match sits; empty seqid when nothing is shared at all. */
  longest_exact_seqid: string;
  longest_exact_start: number;
  longest_exact_end: number;
  longest_exact_qry_start: number;
  longest_exact_qry_end: number;
}

/** The best reference-protein match of one predicted gene inside a
 * gained region: the only name an unannotated query genome can be
 * given locally. */
export interface OrfMatch {
  locus_tag: string;
  /** RefSeq protein accession, for the NCBI protein link. */
  protein_id: string;
  /** Gene symbol or product, whichever the annotation carries. */
  label: string;
  /** Amino-acid identity over the aligned part, percent. */
  identity: number;
  /** Share of the predicted gene covered, percent. */
  coverage: number;
  evalue: number;
}

/** One predicted gene of a gained region with its nucleotide sequence,
 * so novel ones can be linked out to NCBI BLAST. */
export interface IdentifiedOrf {
  start: number;
  end: number;
  strand: number;
  /** null means no similar protein in the reference - expected for
   * the true gains. */
  match: OrfMatch | null;
  seq: string;
}

/** The gene-level answer to "what is gained here": each predicted gene
 * of one region, named by its best match among the reference's own
 * proteins (blastx). */
export interface GainedIdentify {
  qry_seqid: string;
  start: number;
  end: number;
  orfs: IdentifiedOrf[];
  region_seq: string;
}

export interface PanelRow {
  gene_id: string;
  qlen: number;
  cov_pct: number;
  identity: number;
  best_evalue: string;
  call: Call;
}

export interface MatrixRow {
  locus_tag: string;
  symbol: string;
  biotype: string;
  seqid: string;
  start: number;
  end: number;
  calls: Call[];
  cov_pcts: number[];
}

export interface Page<T> {
  rows: T[];
  total: number;
}

export interface ParamSpec {
  name: string;
  label: string;
  help: string;
  layer: "align" | "postprocess";
  kind:
    | { kind: "int"; default: number; min: number; max: number }
    | { kind: "float"; default: number; min: number; max: number }
    | { kind: "optional_int"; default: number | null; min: number; max: number }
    | { kind: "bool"; default: boolean };
}

export interface RunParams {
  min_gap: number;
  present_cov: number;
  partial_cov: number;
  blast_cov: number;
  blast_pid: number;
  blast_evalue: number;
  nucmer_minmatch: number | null;
  nucmer_breaklen: number | null;
  dnadiff: boolean;
}

export interface WgaGene {
  locus_tag: string;
  symbol: string;
  biotype: string;
  seqid: string;
  start: number;
  end: number;
  strand: number;
  product: string;
}

export interface WgaBlock {
  ref_seqid: string;
  ref_start: number;
  ref_end: number;
  qry_seqid: string;
  qry_start: number;
  qry_end: number;
  qry_rev: boolean;
  identity: number;
}

export interface WgaQuery {
  query_id: number;
  query_name: string;
  blocks: WgaBlock[];
  /** Presence call per gene, aligned with WgaData.genes. */
  calls?: Call[];
  /** Coverage percent per gene, aligned with WgaData.genes. */
  cov_pcts?: number[];
  /** Best block identity per gene, aligned with WgaData.genes. */
  identities?: number[];
}

export interface WgaData {
  reference: [string, number][];
  genes: WgaGene[];
  queries: WgaQuery[];
}

/** SNP at a reference position: bases are ASCII char codes. */
export interface SnpEvent {
  pos: number;
  r: number;
  q: number;
}

/** Reference bases missing from the query (deletion in the query). */
export interface DelEvent {
  pos: number;
  len: number;
}

/** Query bases inserted after the reference position `pos`. */
export interface InsEvent {
  pos: number;
  seq: string;
}

export interface AlignmentEvents {
  snps: SnpEvent[];
  dels: DelEvent[];
  ins: InsEvent[];
}

export interface AlignmentQuery {
  query_id: number;
  query_name: string;
  blocks: WgaBlock[];
  events: Record<string, AlignmentEvents>;
}

export interface AlignmentData {
  reference: [string, number][];
  queries: AlignmentQuery[];
}

/** Wire form of AlignmentEvents: parallel arrays instead of one
 * object per event (a divergent query carries >100k SNPs). api.ts
 * hydrates these into the interfaces above once, off the draw path. */
export interface AlignmentEventsWire {
  snp_pos: number[];
  snp_ref: number[];
  snp_qry: number[];
  del_pos: number[];
  del_len: number[];
  ins_pos: number[];
  ins_seq: string[];
}

export interface AlignmentQueryWire {
  query_id: number;
  query_name: string;
  blocks: WgaBlock[];
  events: Record<string, AlignmentEventsWire>;
}

export interface AlignmentDataWire {
  reference: [string, number][];
  queries: AlignmentQueryWire[];
}

/** Reference bases of a window, for the alignment viewer letters mode. */
export interface RefseqWindow {
  seqid: string;
  start: number;
  end: number;
  seq: string;
}

export interface GeneBlock {
  ref_start: number;
  ref_end: number;
  qry_start: number;
  qry_end: number;
  qry_rev: boolean;
  identity: number;
  ref_seq: string;
  qry_seq: string;
}

export interface GeneQueryAlignment {
  query_id: number;
  query_name: string;
  call: Call;
  cov_pct: number;
  best_identity: number;
  mismatches: number;
  indels: number;
  blocks: GeneBlock[];
  unaligned: [number, number][];
  premature_stops: { codon_index: number; aa_position: number }[];
}

export interface GeneDetail {
  locus_tag: string;
  symbol: string;
  protein_id: string;
  biotype: string;
  seqid: string;
  start: number;
  end: number;
  strand: number;
  length: number;
  reference_seq: string;
  queries: GeneQueryAlignment[];
}

export interface RunFile {
  name: string;
  friendly: string;
  size: number;
  human_size: string;
}

export interface TableQuery {
  query_id?: number;
  page?: number;
  page_size?: number;
  sort_by?: string;
  sort_dir?: string;
  search?: string;
  call?: string;
  cols?: string;
}

/** NCBI's nucleotide record for a sequence, scrolled/zoomed to one range. */
export function nuccoreRangeUrl(seqid: string, start: number, end: number): string {
  return `https://www.ncbi.nlm.nih.gov/nuccore/${encodeURIComponent(seqid)}?report=graph&from=${start}&to=${end}`;
}

export function formatSize(bytes: number): string {
  if (bytes < 1024) return `${bytes} B`;
  if (bytes < 1024 * 1024) return `${(bytes / 1024).toFixed(1)} KB`;
  if (bytes < 1024 * 1024 * 1024) return `${(bytes / 1024 / 1024).toFixed(1)} MB`;
  return `${(bytes / 1024 / 1024 / 1024).toFixed(1)} GB`;
}

export function formatDuration(fromIso: string | null, toIso: string | null): string {
  if (!fromIso) return "-";
  const from = new Date(fromIso).getTime();
  const to = toIso ? new Date(toIso).getTime() : Date.now();
  const s = Math.max(0, Math.round((to - from) / 1000));
  if (s < 60) return `${s}s`;
  const m = Math.floor(s / 60);
  return `${m}m ${s % 60}s`;
}

export function formatDate(iso: string): string {
  const d = new Date(iso.endsWith("Z") ? iso : `${iso}Z`);
  if (isNaN(d.getTime())) return iso;
  return d.toLocaleString(undefined, {
    year: "numeric",
    month: "short",
    day: "numeric",
    hour: "2-digit",
    minute: "2-digit",
  });
}
