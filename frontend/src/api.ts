import type {
  AlignmentData,
  AlignmentDataWire,
  GeneDetail,
  GainedRow,
  GapRow,
  GainedIdentify,
  GainedVerify,
  GeneCoverageRow,
  MatrixRow,
  Page,
  PanelRow,
  Project,
  ProjectFile,
  RefseqWindow,
  Run,
  RunFile,
  RunParams,
  RunQuery,
  RunWithLogs,
  TableQuery,
  WgaData,
} from "./types";

const BASE = "/api";

async function request<T>(path: string, init?: RequestInit): Promise<T> {
  const resp = await fetch(`${BASE}${path}`, init);
  if (!resp.ok) {
    let msg = `The server answered with an error (${resp.status}).`;
    try {
      const body = await resp.json();
      if (typeof body.error === "string") {
        msg = body.error;
        try {
          const parsed = JSON.parse(body.error);
          if (parsed && typeof parsed === "object") {
            msg = Object.entries(parsed)
              .map(([field, m]) => `${field}: ${m}`)
              .join(" ");
          }
        } catch {
          /* plain message */
        }
      }
    } catch {
      /* keep default */
    }
    throw new Error(msg);
  }
  // Not every endpoint answers with JSON: the delete routes used to reply with
  // a bare "deleted", and resp.json() then threw AFTER the server had already
  // done the work. The caller saw "Unexpected token 'd'" and assumed the
  // delete had failed, when it had succeeded. Only parse when the server says
  // it sent JSON.
  const contentType = resp.headers.get("content-type") ?? "";
  if (resp.status === 204 || !contentType.includes("application/json")) {
    return undefined as T;
  }
  return resp.json() as Promise<T>;
}

/**
 * Alignment payloads are big (per-base variant events of every query)
 * and immutable for a finished run, so one fetch serves every consumer
 * on the page: the strain map's deep-zoom variant layer and the
 * alignment viewer share it, and toggling between them never refetches.
 * The in-flight promise is cached too, so simultaneous callers coalesce.
 * Failures are evicted, letting a later call retry.
 */
const alignmentCache = new Map<number, Promise<AlignmentData>>();

/**
 * Every gained region of one query, as one list rather than pages. The
 * gained table and the strain map's gained layer want exactly the same
 * rows, so one fetch serves both and switching between them is free.
 * Small - a few hundred rows at worst - so no eviction is needed.
 */
const gainedCache = new Map<string, Promise<GainedRow[]>>();

function gainedAll(runId: number, queryId: number): Promise<GainedRow[]> {
  const key = `${runId}:${queryId}`;
  const hit = gainedCache.get(key);
  if (hit) return hit;
  const p = (async () => {
    const rows: GainedRow[] = [];
    for (let page = 0; ; page++) {
      const d = await request<Page<GainedRow>>(
        `/runs/${runId}/gained${qs({ query_id: queryId, page, page_size: 1000 })}`,
      );
      rows.push(...d.rows);
      if (d.rows.length === 0 || rows.length >= d.total) break;
    }
    return rows;
  })();
  // Let a later call retry rather than caching the failure forever.
  p.catch(() => gainedCache.delete(key));
  gainedCache.set(key, p);
  return p;
}
/** Payloads are tens of MB each for big runs: keep the cache shallow. */
const ALIGNMENT_CACHE_MAX = 4;

/** Turn the columnar wire form into the per-event objects the
 * consumers already index. Runs once per fetch, inside the cached
 * promise, so toggling views never repeats it. */
function hydrateAlignment(w: AlignmentDataWire): AlignmentData {
  return {
    reference: w.reference,
    queries: w.queries.map((q) => ({
      query_id: q.query_id,
      query_name: q.query_name,
      blocks: q.blocks,
      events: Object.fromEntries(
        Object.entries(q.events).map(([seqid, ev]) => [
          seqid,
          {
            snps: ev.snp_pos.map((pos, i) => ({
              pos,
              r: ev.snp_ref[i],
              q: ev.snp_qry[i],
            })),
            dels: ev.del_pos.map((pos, i) => ({ pos, len: ev.del_len[i] })),
            ins: ev.ins_pos.map((pos, i) => ({ pos, seq: ev.ins_seq[i] })),
          },
        ]),
      ),
    })),
  };
}

function alignmentCached(runId: number): Promise<AlignmentData> {
  let p = alignmentCache.get(runId);
  if (p) {
    // refresh recency
    alignmentCache.delete(runId);
    alignmentCache.set(runId, p);
    return p;
  }
  p = request<AlignmentDataWire>(`/runs/${runId}/alignment`)
    .then(hydrateAlignment)
    .catch((e) => {
      alignmentCache.delete(runId);
      throw e;
    });
  alignmentCache.set(runId, p);
  while (alignmentCache.size > ALIGNMENT_CACHE_MAX) {
    const oldest = alignmentCache.keys().next().value!;
    alignmentCache.delete(oldest);
  }
  return p;
}

/** Gene detail payloads are immutable for a finished run, and the strain
 * map's click-to-align panel, the results-table preview and the gene
 * dialog all ask for the same gene; coalesce them into one fetch. */
const geneDetailCache = new Map<string, Promise<GeneDetail>>();

function geneDetailCached(runId: number, locus: string): Promise<GeneDetail> {
  const key = `${runId}:${locus}`;
  let p = geneDetailCache.get(key);
  if (!p) {
    p = request<GeneDetail>(`/runs/${runId}/gene/${encodeURIComponent(locus)}`);
    p.catch(() => geneDetailCache.delete(key));
    geneDetailCache.set(key, p);
  }
  return p;
}

/** Reference back-checks are immutable for a finished run and cost a
 * blast search on the server, so each region is fetched at most once
 * per session and repeated opens of its verdict panel coalesce. */
const gainedVerifyCache = new Map<string, Promise<GainedVerify>>();
const gainedIdentifyCache = new Map<string, Promise<GainedIdentify>>();
const gainedSeqCache = new Map<string, Promise<Map<string, string>>>();

/** Gene identification is the same trade as the back-check - immutable
 * for a finished run, a blast search per call - so it coalesces the
 * same way. */
export function gainedIdentifyCached(
  runId: number,
  queryId: number,
  seqid: string,
  start: number,
  end: number,
): Promise<GainedIdentify> {
  const key = `${runId}:${queryId}:${seqid}:${start}`;
  let p = gainedIdentifyCache.get(key);
  if (!p) {
    p = request<GainedIdentify>(
      `/runs/${runId}/gained/identify?query_id=${queryId}` +
        `&seqid=${encodeURIComponent(seqid)}&start=${start}&end=${end}`,
    );
    p.catch(() => gainedIdentifyCache.delete(key));
    gainedIdentifyCache.set(key, p);
  }
  return p;
}

function gainedVerifyCached(
  runId: number,
  queryId: number,
  seqid: string,
  start: number,
  end: number,
): Promise<GainedVerify> {
  const key = `${runId}:${queryId}:${seqid}:${start}`;
  let p = gainedVerifyCache.get(key);
  if (!p) {
    p = request<GainedVerify>(
      `/runs/${runId}/gained/verify?query_id=${queryId}` +
        `&seqid=${encodeURIComponent(seqid)}&start=${start}&end=${end}`,
    );
    p.catch(() => gainedVerifyCache.delete(key));
    gainedVerifyCache.set(key, p);
  }
  return p;
}

export const api = {
  listProjects: () => request<Project[]>("/projects"),
  createProject: (name: string, organism?: string) =>
    request<Project>("/projects", {
      method: "POST",
      headers: { "content-type": "application/json" },
      body: JSON.stringify(organism ? { name, organism } : { name }),
    }),
  getProject: (id: number) => request<Project>(`/projects/${id}`),
  renameProject: (id: number, name: string) =>
    request<Project>(`/projects/${id}`, {
      method: "PATCH",
      headers: { "content-type": "application/json" },
      body: JSON.stringify({ name }),
    }),
  deleteProject: (id: number) =>
    request<void>(`/projects/${id}`, { method: "DELETE" }),

  listFiles: (projectId: number) =>
    request<ProjectFile[]>(`/projects/${projectId}/files`),
  uploadReference: (projectId: number, fasta: File, gff: File) => {
    const form = new FormData();
    form.append("fasta", fasta);
    form.append("gff", gff);
    return request<ProjectFile[]>(`/projects/${projectId}/reference`, {
      method: "POST",
      body: form,
    });
  },
  fetchReferenceFromNcbi: (projectId: number, accession: string) =>
    request<unknown>(`/projects/${projectId}/reference/ncbi`, {
      method: "POST",
      headers: { "content-type": "application/json" },
      body: JSON.stringify({ accession }),
    }),
  uploadQueries: (projectId: number, files: File[]) => {
    const form = new FormData();
    for (const f of files) form.append("files", f);
    return request<ProjectFile[]>(`/projects/${projectId}/queries`, {
      method: "POST",
      body: form,
    });
  },
  uploadPanel: (projectId: number, file: File) => {
    const form = new FormData();
    form.append("file", file);
    return request<ProjectFile>(`/projects/${projectId}/panel`, {
      method: "POST",
      body: form,
    });
  },
  uploadPanelIds: (projectId: number, file: File) => {
    const form = new FormData();
    form.append("file", file);
    return request<{ file: ProjectFile; found: string[]; from_ncbi: string[]; missing: string[] }>(
      `/projects/${projectId}/panel/from_ids`,
      { method: "POST", body: form },
    );
  },
  buildPanelFromText: (projectId: number, text: string) =>
    request<{ file: ProjectFile; found: string[]; from_ncbi: string[]; missing: string[] }>(
      `/projects/${projectId}/panel/from_text`,
      {
        method: "POST",
        headers: { "content-type": "application/json" },
        body: JSON.stringify({ text }),
      },
    ),
  deleteFile: (projectId: number, fileId: number) =>
    request<void>(`/projects/${projectId}/files/${fileId}`, {
      method: "DELETE",
    }),

  startRun: (projectId: number, queryIds: number[], params: RunParams | null) =>
    request<Run>(`/projects/${projectId}/runs`, {
      method: "POST",
      headers: { "content-type": "application/json" },
      body: JSON.stringify({ query_ids: queryIds, params }),
    }),
  listRuns: (projectId: number) =>
    request<Run[]>(`/projects/${projectId}/runs`),
  getRun: (runId: number) => request<RunWithLogs>(`/runs/${runId}`),
  renameRun: (runId: number, name: string) =>
    request<Run>(`/runs/${runId}/name`, {
      method: "PUT",
      headers: { "content-type": "application/json" },
      body: JSON.stringify({ name }),
    }),
  deleteRun: (runId: number) =>
    request<void>(`/runs/${runId}`, { method: "DELETE" }),
  getRunParams: (runId: number) =>
    request<{ params: RunParams; defaults: RunParams; schema: ParamSpecLike[] }>(
      `/runs/${runId}/params`,
    ),
  getPresets: () =>
    request<{ schema: ParamSpecLike[]; presets: string[] }>("/presets"),

  genesCoverage: (runId: number, q: TableQuery) =>
    request<Page<GapRow | GeneCoverageRow>>(
      `/runs/${runId}/genes_coverage${qs(q)}`,
    ) as unknown as Promise<Page<GeneCoverageRow>>,
  unalignedGaps: (runId: number, q: TableQuery) =>
    request<Page<GapRow>>(`/runs/${runId}/unaligned_gaps${qs(q)}`),
  gained: (runId: number, q: TableQuery) =>
    request<Page<GainedRow>>(`/runs/${runId}/gained${qs(q)}`),
  gainedAll,
  panelRecheck: (runId: number, q: TableQuery) =>
    request<Page<PanelRow>>(`/runs/${runId}/panel_recheck${qs(q)}`),
  matrix: (runId: number, q: TableQuery) =>
    request<Page<MatrixRow>>(`/runs/${runId}/matrix${qs(q)}`),
  wga: (runId: number) => request<WgaData>(`/runs/${runId}/wga`),
  alignment: alignmentCached,
  refseq: (runId: number, seqid: string, start: number, end: number) =>
    request<RefseqWindow>(
      `/runs/${runId}/refseq?seqid=${encodeURIComponent(seqid)}&start=${start}&end=${end}`,
    ),
  geneDetail: geneDetailCached,
  gainedVerify: gainedVerifyCached,
  gainedIdentify: gainedIdentifyCached,
  /** The sequence of every gained region of one query, keyed by
   * "seqid:start-end": one immutable request per table view, backing
   * the copy-to-clipboard column. */
  gainedSequences: (runId: number, queryId: number) => {
    const key = `${runId}:${queryId}`;
    let p = gainedSeqCache.get(key);
    if (!p) {
      p = request<{ seqid: string; start: number; end: number; seq: string }[]>(
        `/runs/${runId}/gained/sequences?query_id=${queryId}`,
      ).then((rows) => {
        const m = new Map<string, string>();
        for (const r of rows) m.set(`${r.seqid}:${r.start}-${r.end}`, r.seq);
        return m;
      });
      p.catch(() => gainedSeqCache.delete(key));
      gainedSeqCache.set(key, p);
    }
    return p;
  },
  listRunFiles: (runId: number) =>
    request<{ run_id: number; status: string; files: RunFile[] }>(
      `/runs/${runId}/files`,
    ),

  exportUrl: (runId: number, table: string, q: TableQuery, format: string) =>
    `${BASE}/runs/${runId}/export/${table}${qs({ ...q, page: undefined, page_size: undefined })}&format=${format}`,
  geneExportUrl: (runId: number, locus: string, format: string) =>
    `${BASE}/runs/${runId}/gene/${encodeURIComponent(locus)}/export?format=${format}`,
  runFileUrl: (runId: number, name: string, download: boolean) =>
    `${BASE}/runs/${runId}/files/${encodeURIComponent(name)}${download ? "?download=1" : ""}`,

  getSettings: () =>
    request<{ has_ncbi_api_key: boolean; ncbi_api_key: string | null }>(
      "/settings",
    ),
  putNcbiKey: (key: string) =>
    request<{ status: string; masked: string }>("/settings/ncbi_api_key", {
      method: "PUT",
      headers: { "content-type": "application/json" },
      body: JSON.stringify({ key }),
    }),
  deleteNcbiKey: () =>
    request<{ status: string }>("/settings/ncbi_api_key", {
      method: "DELETE",
    }),
  usage: (projectId: number) =>
    request<{ bytes: number; human: string }>(`/projects/${projectId}/usage`),
};

export interface ParamSpecLike {
  name: string;
  label: string;
  help: string;
  layer: "align" | "postprocess";
  kind: Record<string, unknown> & { kind: string };
}

function qs(q: TableQuery): string {
  const parts: string[] = [];
  if (q.query_id !== undefined) parts.push(`query_id=${q.query_id}`);
  if (q.page !== undefined) parts.push(`page=${q.page}`);
  if (q.page_size !== undefined) parts.push(`page_size=${q.page_size}`);
  if (q.sort_by) parts.push(`sort_by=${encodeURIComponent(q.sort_by)}`);
  if (q.sort_dir) parts.push(`sort_dir=${encodeURIComponent(q.sort_dir)}`);
  if (q.search) parts.push(`search=${encodeURIComponent(q.search)}`);
  if (q.call) parts.push(`call=${encodeURIComponent(q.call)}`);
  if (q.cols) parts.push(`cols=${encodeURIComponent(q.cols)}`);
  return parts.length ? `?${parts.join("&")}` : "";
}

export type { RunQuery };
