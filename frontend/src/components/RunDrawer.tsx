import { useEffect, useRef, useState } from "react";
import { api } from "../api";
import type { Run, RunFile } from "../types";
import { formatDuration } from "../types";
import { Spinner } from "./ui";

/**
 * The run drawer: live progress, elapsed time and tool logs for a running
 * (or finished) run. Polls while the run is queued or running.
 */
export function RunDrawer({
  runId,
  onClose,
  onFinished,
}: {
  runId: number | null;
  onClose: () => void;
  onFinished?: (run: Run) => void;
}) {
  const [run, setRun] = useState<Run | null>(null);
  const [logs, setLogs] = useState<string[]>([]);
  const logRef = useRef<HTMLDivElement>(null);
  const finishedRef = useRef(false);

  useEffect(() => {
    if (runId === null) return;
    setRun(null);
    setLogs([]);
    finishedRef.current = false;
    let stop = false;
    async function poll() {
      while (!stop) {
        try {
          const data = await api.getRun(runId!);
          setRun(data.run);
          setLogs(data.logs);
          if (
            (data.run.status === "succeeded" || data.run.status === "failed") &&
            !finishedRef.current
          ) {
            finishedRef.current = true;
            onFinished?.(data.run);
          }
          if (data.run.status === "succeeded" || data.run.status === "failed") {
            return;
          }
        } catch {
          /* keep polling */
        }
        await new Promise((r) => setTimeout(r, 1500));
      }
    }
    poll();
    return () => {
      stop = true;
    };
  }, [runId, onFinished]);

  useEffect(() => {
    logRef.current?.scrollTo({ top: logRef.current.scrollHeight });
  }, [logs]);

  if (runId === null) return null;
  const running = run?.status === "queued" || run?.status === "running";

  return (
    <div className="fixed bottom-0 right-0 left-0 sm:left-auto sm:right-6 sm:bottom-6 z-40 w-full sm:w-[420px] bg-white border border-zinc-200 rounded-t-xl sm:rounded-xl shadow-2xl overflow-hidden">
      <div className="flex items-center justify-between px-5 h-14 border-b border-zinc-200">
        <div className="flex items-center gap-3 min-w-0">
          {running ? (
            <Spinner className="text-zinc-500" />
          ) : run?.status === "succeeded" ? (
            <span className="text-emerald-600">{"\u2713"}</span>
          ) : run?.status === "failed" ? (
            <span className="text-red-600">{"\u2715"}</span>
          ) : null}
          <span className="font-medium truncate">
            {run
              ? run.status === "succeeded"
                ? "Comparison finished"
                : run.status === "failed"
                  ? "Comparison failed"
                  : "Comparing genomes"
              : "Starting..."}
          </span>
          {run && running && (
            <span className="text-sm text-zinc-400 whitespace-nowrap">
              {formatDuration(run.started_at, null)} elapsed
            </span>
          )}
        </div>
        <button
          onClick={onClose}
          className="w-9 h-9 grid place-items-center rounded-md text-zinc-400 hover:text-zinc-900 hover:bg-zinc-100"
          aria-label="Close"
        >
          <svg width="14" height="14" viewBox="0 0 16 16" fill="none">
            <path d="M3 3l10 10M13 3L3 13" stroke="currentColor" strokeWidth="1.8" strokeLinecap="round" />
          </svg>
        </button>
      </div>
      <div className="px-5 py-3 max-h-56 overflow-y-auto thin-scroll" ref={logRef}>
        {run?.step && running && (
          <p className="text-[15px] text-zinc-800 mb-2 flex items-center gap-2">
            <span className="w-1.5 h-1.5 rounded-full bg-zinc-900 animate-pulse" />
            {run.step}
          </p>
        )}
        {run?.status === "failed" && run.error && (
          <p className="text-[15px] text-red-700 mb-2">{run.error}</p>
        )}
        {logs.length === 0 && (
          <p className="text-sm text-zinc-400">Waiting for the first log line...</p>
        )}
        {logs.map((l, i) => (
          <p key={i} className="text-xs text-zinc-500 font-mono leading-relaxed">
            {l}
          </p>
        ))}
      </div>
    </div>
  );
}

/** The Files panel: every artifact of a run, viewable and downloadable. */
export function FilesPanel({ runId }: { runId: number }) {
  const [files, setFiles] = useState<RunFile[]>([]);
  const [error, setError] = useState<string | null>(null);

  useEffect(() => {
    api
      .listRunFiles(runId)
      .then((f) => setFiles(f.files))
      .catch((e) => setError((e as Error).message));
  }, [runId]);

  if (error) return <p className="text-red-700 text-[15px]">{error}</p>;
  if (files.length === 0)
    return <p className="text-zinc-500 text-[15px]">No files yet.</p>;

  return (
    <ul className="divide-y divide-zinc-100 border border-zinc-200 rounded-lg">
      {files.map((f) => (
        <li key={f.name} className="flex items-center justify-between px-4 py-3">
          <div className="min-w-0 mr-4">
            <p className="font-medium text-[15px]">{f.friendly}</p>
            <p className="text-xs text-zinc-400 mt-0.5">
              {f.name} - {f.human_size}
            </p>
          </div>
          <div className="flex gap-2 shrink-0">
            <a
              href={api.runFileUrl(runId, f.name, false)}
              target="_blank"
              rel="noreferrer"
              className="h-10 px-3 inline-flex items-center rounded-lg border border-zinc-300 text-[15px] hover:bg-zinc-100"
            >
              View
            </a>
            <a
              href={api.runFileUrl(runId, f.name, true)}
              className="h-10 px-3 inline-flex items-center rounded-lg border border-zinc-300 text-[15px] hover:bg-zinc-100"
            >
              Download
            </a>
          </div>
        </li>
      ))}
    </ul>
  );
}
