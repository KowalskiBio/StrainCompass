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
}

export interface WgaData {
  reference: [string, number][];
  genes: WgaGene[];
  queries: WgaQuery[];
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
