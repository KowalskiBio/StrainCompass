import { useCallback, useEffect, useMemo, useState } from "react";
import { useNavigate, useParams, useSearchParams } from "react-router-dom";
import { api } from "../api";
import type { Project, ProjectFile, Run, RunParams } from "../types";
import { formatSize, formatDuration, formatDate } from "../types";
import {
  Button,
  ErrorBox,
  Spinner,
  Tabs,
} from "../components/ui";
import { FilesPanel, RunDrawer } from "../components/RunDrawer";
import { InputWizard } from "../components/InputWizard";
import { ResultsTables, type TableKind } from "../components/ResultsTables";
import { GenomeView } from "../components/GenomeView";
import { GeneMsaDialog } from "../components/GeneMsaDialog";

type Tab = "inputs" | "runs" | "table" | "genome";

export default function ProjectPage() {
  const { id } = useParams();
  const projectId = Number(id);
  const navigate = useNavigate();
  const [params, setParams] = useSearchParams();

  const [project, setProject] = useState<Project | null>(null);
  const [files, setFiles] = useState<ProjectFile[]>([]);
  const [runs, setRuns] = useState<Run[]>([]);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);

  const [wizardOpen, setWizardOpen] = useState(false);
  const [renameOpen, setRenameOpen] = useState(false);
  const [deleteOpen, setDeleteOpen] = useState(false);
  const [drawerRunId, setDrawerRunId] = useState<number | null>(null);

  const tab = (params.get("tab") as Tab) ?? "inputs";
  const runId = params.get("run") ? Number(params.get("run")) : null;
  const gene = params.get("gene");

  // Applies every key in `updates` to the URL in one history entry. Calling
  // setParams multiple times in a row (once per key) is unsafe: each call
  // starts from the searchParams snapshot of when it was made, so a later
  // call in the same tick can clobber an earlier one's change instead of
  // building on it - this is how switching a table's query used to silently
  // revert an unrelated "which sub-table is active" param.
  const setParams2 = useCallback(
    (updates: Record<string, string | null>) => {
      setParams(
        (prev) => {
          const next = new URLSearchParams(prev);
          for (const [key, value] of Object.entries(updates)) {
            if (value === null || value === "") next.delete(key);
            else next.set(key, value);
          }
          return next;
        },
        { replace: true },
      );
    },
    [setParams],
  );

  const setParam = useCallback(
    (key: string, value: string | null) => setParams2({ [key]: value }),
    [setParams2],
  );

  const reload = useCallback(() => {
    if (!Number.isFinite(projectId)) return;
    api
      .getProject(projectId)
      .then((p) => {
        setProject(p);
        setError(null);
      })
      .catch((e) => setError((e as Error).message));
    api
      .listFiles(projectId)
      .then((f) => setFiles(f))
      .catch(() => {});
    api
      .listRuns(projectId)
      .then((r) => {
        setRuns(r);
        // pick the newest succeeded run if none is selected
        if (!params.get("run") && r.length > 0) {
          const best =
            r.find((x) => x.status === "succeeded") ?? r[0];
          setParam("run", String(best.id));
        }
      })
      .catch(() => {})
      .finally(() => setLoading(false));
  }, [projectId, params, setParam]);

  useEffect(() => {
    reload();
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [projectId]);

  // open the wizard automatically for a fresh project without a reference
  useEffect(() => {
    if (!loading && project && !project.has_reference && files.length === 0) {
      setWizardOpen(true);
    }
  }, [loading, project, files.length]);

  // the genome-view gene deep link is consumed on mount, then dropped so
  // re-visiting the tab does not jump back to that gene forever
  useEffect(() => {
    if (tab === "genome" && params.get("ggene")) {
      const t = setTimeout(() => setParam("ggene", null), 100);
      return () => clearTimeout(t);
    }
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [tab]);

  const selectedRun = useMemo(
    () => runs.find((r) => r.id === runId) ?? null,
    [runs, runId],
  );

  // poll while any run is queued or running
  useEffect(() => {
    const active = runs.some(
      (r) => r.status === "queued" || r.status === "running",
    );
    if (!active) return;
    const t = setInterval(() => {
      api
        .listRuns(projectId)
        .then(setRuns)
        .catch(() => {});
    }, 2500);
    return () => clearInterval(t);
  }, [runs, projectId]);

  if (loading)
    return (
      <div className="flex items-center gap-3 text-zinc-500 py-24 justify-center dark:text-zinc-400">
        <Spinner /> Loading the project...
      </div>
    );
  if (error)
    return (
      <div className="max-w-3xl mx-auto px-4 py-8">
        <ErrorBox message={error} />
      </div>
    );
  if (!project)
    return (
      <div className="max-w-3xl mx-auto px-4 py-8">
        <ErrorBox message="This project does not exist (anymore)." />
      </div>
    );

  const refFasta = files.find((f) => f.role === "reference_fasta");
  const queries = files.filter((f) => f.role === "query");
  const canRun = Boolean(refFasta) && queries.length > 0;

  function openTab(t: string) {
    setParam("tab", t === "inputs" ? null : (t as Tab));
  }

  function startRun(queryIds: number[], runParams: RunParams) {
    api
      .startRun(projectId, queryIds, runParams)
      .then((run) => {
        reload();
        setDrawerRunId(run.id);
        setParam("run", String(run.id));
      })
      .catch((e) => setError((e as Error).message));
  }

  return (
    <div className="max-w-[1600px] mx-auto px-4 py-6">
      {/* header */}
      <div className="flex flex-wrap items-start justify-between gap-4 mb-4">
        <div className="min-w-0">
          <div className="flex items-center gap-3">
            <h1 className="text-2xl font-semibold tracking-tight truncate">
              {project.name}
            </h1>
            <span className="text-xs text-zinc-400 font-mono dark:text-zinc-500">#{project.id}</span>
          </div>
          <p className="text-zinc-500 mt-1 text-[15px] dark:text-zinc-400">
            {project.organism && project.organism !== "bacteria"
              ? `${project.organism} - `
              : ""}
            {project.has_reference ? "reference ready" : "no reference yet"} -{" "}
            {project.n_queries} quer{project.n_queries === 1 ? "y" : "ies"} -{" "}
            {project.n_runs} run{project.n_runs === 1 ? "" : "s"} -{" "}
            {formatSize(project.usage_bytes)} on the server
          </p>
        </div>
        <div className="flex items-center gap-2">
          <Button variant="secondary" onClick={() => setRenameOpen(true)}>
            Rename
          </Button>
          <Button variant="danger" onClick={() => setDeleteOpen(true)}>
            Delete
          </Button>
          {canRun && (
            <Button
              onClick={() => {
                setWizardOpen(true);
              }}
            >
              New comparison
            </Button>
          )}
        </div>
      </div>

      {error && (
        <div className="mb-4">
          <ErrorBox message={error} />
        </div>
      )}

      <Tabs
        tabs={[
          { key: "inputs", label: "Inputs" },
          { key: "runs", label: `Runs (${runs.length})` },
          {
            key: "table",
            label: "Results table",
            disabled: !selectedRun || selectedRun.status !== "succeeded",
          },
          {
            key: "genome",
            label: "Genome view",
            disabled: !selectedRun || selectedRun.status !== "succeeded",
          },
        ]}
        active={tab}
        onChange={openTab}
      />

      <div className="py-6">
        {tab === "inputs" && (
          <InputsTab
            projectId={projectId}
            files={files}
            onChanged={reload}
            onOpenWizard={() => setWizardOpen(true)}
            canRun={canRun}
          />
        )}
        {tab === "runs" && (
          <RunsTab
            runs={runs}
            selectedRun={selectedRun}
            onSelect={(r) => {
              setParam("run", String(r.id));
              reload();
            }}
            onDelete={(r) => {
              api.deleteRun(r.id).then(reload).catch((e) => setError((e as Error).message));
            }}
          />
        )}
        {tab === "table" && selectedRun && (
          <ResultsTables
            run={selectedRun}
            initialTable={
              (params.get("table") as TableKind) ?? "genes_coverage"
            }
            initialQueryId={
              params.get("q") ? Number(params.get("q")) : undefined
            }
            onStateChange={(table, queryId) => {
              setParams2({
                table: table === "genes_coverage" ? null : table,
                q: queryId ? String(queryId) : null,
              });
            }}
            onOpenGene={(locus) => setParam("gene", locus)}
          />
        )}
        {tab === "genome" && selectedRun && (
          <GenomeView
            run={selectedRun}
            initialGene={params.get("ggene")}
            initialRange={
              params.get("gstart") && params.get("gend")
                ? {
                    seqid: params.get("gseq") ?? "",
                    start: Number(params.get("gstart")),
                    end: Number(params.get("gend")),
                  }
                : null
            }
            onOpenGene={(locus) => setParam("gene", locus)}
            onRangeChange={(r) => {
              setParams2({
                gseq: r.seqid,
                gstart: String(Math.round(r.start)),
                gend: String(Math.round(r.end)),
              });
            }}
          />
        )}
      </div>

      {wizardOpen && (
        <InputWizard
          open
          onClose={() => setWizardOpen(false)}
          projectId={projectId}
          files={files}
          onFilesChanged={reload}
          onCompare={startRun}
        />
      )}

      <RunDrawer
        runId={drawerRunId}
        onClose={() => setDrawerRunId(null)}
        onFinished={() => reload()}
      />

      <GeneMsaDialog
        runId={selectedRun?.id ?? 0}
        locus={gene}
        onClose={() => setParam("gene", null)}
        onShowInGenome={(locus) => {
          setParam("gene", null);
          setParam("tab", "genome");
          setParam("ggene", locus);
        }}
      />

      {renameOpen && (
        <RenameDialog
          project={project}
          onClose={() => setRenameOpen(false)}
          onDone={() => {
            setRenameOpen(false);
            reload();
          }}
        />
      )}
      {deleteOpen && (
        <DeleteProjectDialog
          project={project}
          onClose={() => setDeleteOpen(false)}
          onDone={() => {
            setDeleteOpen(false);
            navigate("/");
          }}
        />
      )}
    </div>
  );
}

function InputsTab({
  projectId,
  files,
  onChanged,
  onOpenWizard,
  canRun,
}: {
  projectId: number;
  files: ProjectFile[];
  onChanged: () => void;
  onOpenWizard: () => void;
  canRun: boolean;
}) {
  const refFasta = files.find((f) => f.role === "reference_fasta");
  const refGff = files.find((f) => f.role === "reference_gff");
  const queries = files.filter((f) => f.role === "query");
  const panel = files.find((f) => f.role === "panel");

  return (
    <div className="space-y-6 max-w-4xl">
      <section>
        <h2 className="text-lg font-medium mb-2">Reference genome</h2>
        {refFasta && refGff ? (
          <ul className="divide-y divide-zinc-100 border border-zinc-200 rounded-xl bg-white dark:divide-zinc-800 dark:border-zinc-800 dark:bg-zinc-900">
            {[refFasta, refGff].map((f) => (
              <FileRow key={f.id} file={f} projectId={projectId} onChanged={onChanged} />
            ))}
          </ul>
        ) : (
          <EmptyInline
            text="No reference genome yet."
            action="Add reference"
            onClick={onOpenWizard}
          />
        )}
      </section>

      <section>
        <h2 className="text-lg font-medium mb-2">
          Query genomes ({queries.length})
        </h2>
        {queries.length > 0 ? (
          <ul className="divide-y divide-zinc-100 border border-zinc-200 rounded-xl bg-white dark:divide-zinc-800 dark:border-zinc-800 dark:bg-zinc-900">
            {queries.map((f) => (
              <FileRow key={f.id} file={f} projectId={projectId} onChanged={onChanged} />
            ))}
          </ul>
        ) : (
          <EmptyInline
            text="No query genomes yet."
            action="Add queries"
            onClick={onOpenWizard}
          />
        )}
      </section>

      <section>
        <h2 className="text-lg font-medium mb-2">
          Gene panel (optional, for the strict recheck)
        </h2>
        {panel ? (
          <ul className="divide-y divide-zinc-100 border border-zinc-200 rounded-xl bg-white dark:divide-zinc-800 dark:border-zinc-800 dark:bg-zinc-900">
            <FileRow file={panel} projectId={projectId} onChanged={onChanged} />
          </ul>
        ) : (
          <p className="text-zinc-500 text-[15px] dark:text-zinc-400">
            Not used. Add a panel from the setup dialog if you want certain
            genes double-checked with a precise search.
          </p>
        )}
      </section>

      <div className="pt-2">
        <Button size="lg" variant={canRun ? "primary" : "secondary"} onClick={onOpenWizard}>
          {canRun ? "Run a new comparison" : "Set up inputs"}
        </Button>
      </div>
    </div>
  );
}

function FileRow({
  file,
  projectId,
  onChanged,
}: {
  file: ProjectFile;
  projectId: number;
  onChanged: () => void;
}) {
  const [busy, setBusy] = useState(false);
  return (
    <li className="flex items-center justify-between gap-4 px-4 py-3">
      <div className="min-w-0">
        <p className="text-[15px] font-medium truncate">{file.display_name}</p>
        <p className="text-xs text-zinc-400 mt-0.5 dark:text-zinc-500">
          {formatSize(file.size)} - added {formatDate(file.created_at)}
        </p>
      </div>
      <button
        className="text-sm text-zinc-400 hover:text-red-600 h-9 px-2 rounded-md hover:bg-red-50 shrink-0 dark:text-zinc-500 dark:hover:text-red-400 dark:hover:bg-red-950/40"
        disabled={busy}
        onClick={async () => {
          if (
            !window.confirm(
              `Remove ${file.display_name} from this project?`,
            )
          )
            return;
          setBusy(true);
          try {
            await api.deleteFile(projectId, file.id);
            onChanged();
          } finally {
            setBusy(false);
          }
        }}
      >
        Remove
      </button>
    </li>
  );
}

function EmptyInline({
  text,
  action,
  onClick,
}: {
  text: string;
  action: string;
  onClick: () => void;
}) {
  return (
    <div className="flex items-center justify-between gap-4 border border-dashed border-zinc-300 rounded-xl px-4 py-3 dark:border-zinc-700">
      <p className="text-zinc-500 text-[15px] dark:text-zinc-400">{text}</p>
      <Button variant="secondary" onClick={onClick}>
        {action}
      </Button>
    </div>
  );
}

function RunsTab({
  runs,
  selectedRun,
  onSelect,
  onDelete,
}: {
  runs: Run[];
  selectedRun: Run | null;
  onSelect: (r: Run) => void;
  onDelete: (r: Run) => void;
}) {
  const [paramsOf, setParamsOf] = useState<Record<number, RunParams>>({});
  const [showFiles, setShowFiles] = useState<number | null>(null);

  useEffect(() => {
    runs
      .filter((r) => r.status === "succeeded" || r.status === "failed")
      .forEach((r) => {
        if (paramsOf[r.id]) return;
        api
          .getRunParams(r.id)
          .then((p) => setParamsOf((prev) => ({ ...prev, [r.id]: p.params })))
          .catch(() => {});
      });
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [runs]);

  if (runs.length === 0)
    return (
      <p className="text-zinc-500 text-[15px] dark:text-zinc-400">
        No comparisons have been run yet. Set up the inputs and run the first
        one.
      </p>
    );

  return (
    <div className="space-y-4">
      <ul className="divide-y divide-zinc-100 border border-zinc-200 rounded-xl bg-white dark:divide-zinc-800 dark:border-zinc-800 dark:bg-zinc-900">
        {runs.map((r) => (
          <li key={r.id} className="px-4 py-3">
            <div className="flex flex-wrap items-center justify-between gap-3">
              <div className="flex items-center gap-3 min-w-0">
                <StatusDot status={r.status} />
                <div>
                  <p className="text-[15px] font-medium">
                    Run #{r.id}
                    {selectedRun?.id === r.id && (
                      <span className="text-xs text-zinc-400 font-normal dark:text-zinc-500">
                        {" "}
                        (currently shown)
                      </span>
                    )}
                  </p>
                  <p className="text-xs text-zinc-400 mt-0.5 dark:text-zinc-500">
                    {formatDate(r.created_at)}
                    {r.finished_at
                      ? ` - took ${formatDuration(r.started_at, r.finished_at)}`
                      : r.status === "running"
                        ? " - running"
                        : ""}
                    {" - "}
                    {r.queries.length} quer{r.queries.length === 1 ? "y" : "ies"}
                  </p>
                </div>
              </div>
              <div className="flex items-center gap-2 shrink-0">
                {r.status === "succeeded" && (
                  <Button
                    variant={
                      selectedRun?.id === r.id ? "secondary" : "primary"
                    }
                    onClick={() => onSelect(r)}
                  >
                    {selectedRun?.id === r.id ? "Selected" : "Show results"}
                  </Button>
                )}
                {r.status === "succeeded" && (
                  <Button
                    variant="secondary"
                    onClick={() => setShowFiles(showFiles === r.id ? null : r.id)}
                  >
                    Files
                  </Button>
                )}
                <Button
                  variant="ghost"
                  onClick={() => {
                    if (window.confirm(`Delete run #${r.id} and its results?`))
                      onDelete(r);
                  }}
                >
                  Delete
                </Button>
              </div>
            </div>
            {r.status === "failed" && r.error && (
              <p className="text-sm text-red-700 mt-2 dark:text-red-400">{r.error}</p>
            )}
            {paramsOf[r.id] && (
              <p className="text-xs text-zinc-400 mt-2 font-mono dark:text-zinc-500">
                present &gt;= {paramsOf[r.id].present_cov}% coverage, partial
                &gt; {paramsOf[r.id].partial_cov}%, gaps &gt;=
                {paramsOf[r.id].min_gap} bp
                {paramsOf[r.id].dnadiff ? ", dnadiff on" : ""}
              </p>
            )}
          </li>
        ))}
      </ul>
      {showFiles !== null && (
        <div className="border border-zinc-200 rounded-xl bg-white p-4 dark:border-zinc-800 dark:bg-zinc-900">
          <h3 className="font-medium mb-3">Files of run #{showFiles}</h3>
          <FilesPanel runId={showFiles} />
        </div>
      )}
    </div>
  );
}

function StatusDot({ status }: { status: string }) {
  const color =
    status === "succeeded"
      ? "bg-emerald-500"
      : status === "failed"
        ? "bg-red-500"
        : status === "running"
          ? "bg-blue-500 animate-pulse"
          : "bg-zinc-400";
  return <span className={`w-2.5 h-2.5 rounded-full shrink-0 ${color}`} />;
}

function RenameDialog({
  project,
  onClose,
  onDone,
}: {
  project: Project;
  onClose: () => void;
  onDone: () => void;
}) {
  const [name, setName] = useState(project.name);
  const [error, setError] = useState<string | null>(null);
  return (
    <div className="fixed inset-0 z-50 flex items-start justify-center bg-zinc-900/40 p-4">
      <div className="bg-white rounded-xl shadow-xl border border-zinc-200 w-full max-w-md mt-24 dark:bg-zinc-900 dark:border-zinc-800">
        <div className="px-6 py-4 border-b border-zinc-200 dark:border-zinc-800">
          <h2 className="text-lg font-semibold">Rename project</h2>
        </div>
        <div className="p-6 space-y-4">
          {error && <ErrorBox message={error} />}
          <input
            autoFocus
            value={name}
            onChange={(e) => setName(e.target.value)}
            className="w-full h-11 px-3 rounded-lg border border-zinc-300 text-[15px] focus:border-zinc-500 outline-none dark:border-zinc-700 dark:bg-zinc-900"
            onKeyDown={(e) => {
              if (e.key === "Enter" && name.trim()) {
                api
                  .renameProject(project.id, name.trim())
                  .then(onDone)
                  .catch((e2) => setError((e2 as Error).message));
              }
            }}
          />
          <div className="flex justify-end gap-3">
            <Button variant="ghost" onClick={onClose}>
              Cancel
            </Button>
            <Button
              disabled={!name.trim()}
              onClick={() => {
                api
                  .renameProject(project.id, name.trim())
                  .then(onDone)
                  .catch((e2) => setError((e2 as Error).message));
              }}
            >
              Save
            </Button>
          </div>
        </div>
      </div>
    </div>
  );
}

function DeleteProjectDialog({
  project,
  onClose,
  onDone,
}: {
  project: Project;
  onClose: () => void;
  onDone: () => void;
}) {
  const [confirmName, setConfirmName] = useState("");
  const [error, setError] = useState<string | null>(null);
  const [busy, setBusy] = useState(false);
  return (
    <div className="fixed inset-0 z-50 flex items-start justify-center bg-zinc-900/40 p-4">
      <div className="bg-white rounded-xl shadow-xl border border-zinc-200 w-full max-w-md mt-24 dark:bg-zinc-900 dark:border-zinc-800">
        <div className="px-6 py-4 border-b border-zinc-200 dark:border-zinc-800">
          <h2 className="text-lg font-semibold">Delete project</h2>
        </div>
        <div className="p-6 space-y-4">
          {error && <ErrorBox message={error} />}
          <p className="text-[15px] text-zinc-600 dark:text-zinc-400">
            This deletes the project <strong>{project.name}</strong>, its{" "}
            {project.n_runs} run{project.n_runs === 1 ? "" : "s"} and all
            uploaded files ({formatSize(project.usage_bytes)}). This cannot be
            undone.
          </p>
          <p className="text-sm text-zinc-500 dark:text-zinc-400">
            Type the project name to confirm:
          </p>
          <input
            autoFocus
            value={confirmName}
            onChange={(e) => setConfirmName(e.target.value)}
            className="w-full h-11 px-3 rounded-lg border border-zinc-300 text-[15px] focus:border-zinc-500 outline-none dark:border-zinc-700 dark:bg-zinc-900"
          />
          <div className="flex justify-end gap-3">
            <Button variant="ghost" onClick={onClose}>
              Cancel
            </Button>
            <Button
              variant="danger"
              disabled={confirmName !== project.name || busy}
              onClick={() => {
                setBusy(true);
                api
                  .deleteProject(project.id)
                  .then(onDone)
                  .catch((e) => {
                    setError((e as Error).message);
                    setBusy(false);
                  });
              }}
            >
              {busy ? "Deleting..." : "Delete permanently"}
            </Button>
          </div>
        </div>
      </div>
    </div>
  );
}
