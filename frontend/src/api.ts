import type {
  GeneDetail,
  GapRow,
  GeneCoverageRow,
  MatrixRow,
  Page,
  PanelRow,
  Project,
  ProjectFile,
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
  return resp.json() as Promise<T>;
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
    request<string>(`/projects/${id}`, { method: "DELETE" }),

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
    return request<{ file: ProjectFile; found: string[]; missing: string[] }>(
      `/projects/${projectId}/panel/from_ids`,
      { method: "POST", body: form },
    );
  },
  deleteFile: (projectId: number, fileId: number) =>
    request<string>(`/projects/${projectId}/files/${fileId}`, {
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
  deleteRun: (runId: number) =>
    request<string>(`/runs/${runId}`, { method: "DELETE" }),
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
  panelRecheck: (runId: number, q: TableQuery) =>
    request<Page<PanelRow>>(`/runs/${runId}/panel_recheck${qs(q)}`),
  matrix: (runId: number, q: TableQuery) =>
    request<Page<MatrixRow>>(`/runs/${runId}/matrix${qs(q)}`),
  wga: (runId: number) => request<WgaData>(`/runs/${runId}/wga`),
  geneDetail: (runId: number, locus: string) =>
    request<GeneDetail>(
      `/runs/${runId}/gene/${encodeURIComponent(locus)}`,
    ),
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
