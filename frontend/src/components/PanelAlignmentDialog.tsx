import { useEffect, useState } from "react";
import { api } from "../api";
import type { Call, PanelAlignmentView, Run } from "../types";
import { GeneAlignmentTiles } from "./GeneAlignmentPanel";
import { CallBadge, ErrorBox, Modal, Spinner } from "./ui";

/**
 * Why does a panel gene differ between strains? The gene aligned base by
 * base to its hit in one strain (left), next to another strain of the
 * user's choice (right, a Present one by default), both against the same
 * panel sequence, so a Partial call can be read against a full one.
 */
export function PanelAlignmentDialog({
  run,
  geneId,
  queryId,
  onClose,
}: {
  run: Run;
  geneId: string;
  /** The strain the dialog was opened from; null when closed. */
  queryId: number | null;
  onClose: () => void;
}) {
  const [left, setLeft] = useState<PanelAlignmentView | null>(null);
  const [leftError, setLeftError] = useState<string | null>(null);
  const [calls, setCalls] = useState<Map<number, Call> | null>(null);
  const [rightId, setRightId] = useState<number | undefined>(undefined);

  useEffect(() => {
    if (queryId === null) return;
    let cancelled = false;
    setLeft(null);
    setLeftError(null);
    api
      .panelAlignment(run.id, queryId, geneId)
      .then((v) => !cancelled && setLeft(v))
      .catch((e) => !cancelled && setLeftError((e as Error).message));
    return () => {
      cancelled = true;
    };
  }, [run.id, queryId, geneId]);

  // every strain's call, to offer a Present one on the right first
  useEffect(() => {
    if (queryId === null) return;
    let cancelled = false;
    setCalls(null);
    api
      .panelMatrix(run.id, { search: geneId, page_size: 1000 })
      .then((p) => {
        if (cancelled) return;
        const row = p.rows.find((r) => r.gene_id === geneId);
        const m = new Map<number, Call>();
        run.queries.forEach((q, i) => row && m.set(q.file_id, row.calls[i]));
        setCalls(m);
        const present = run.queries.find(
          (q) => q.file_id !== queryId && m.get(q.file_id) === "PRESENT",
        );
        const other = run.queries.find((q) => q.file_id !== queryId);
        setRightId((present ?? other)?.file_id);
      })
      .catch(() => {
        if (cancelled) return;
        setRightId(run.queries.find((q) => q.file_id !== queryId)?.file_id);
      });
    return () => {
      cancelled = true;
    };
  }, [run.id, run.queries, queryId, geneId]);

  return (
    <Modal
      open={queryId !== null}
      onClose={onClose}
      extraWide
      title={
        <span className="flex items-baseline gap-3 flex-wrap">
          Investigate <span className="font-mono">{geneId}</span>
          {left && (
            <span className="text-sm font-normal text-zinc-500 dark:text-zinc-400">
              aligned to {left.variant_source || left.variant} ({left.panel_len.toLocaleString("en-US")}{" "}
              bp)
            </span>
          )}
        </span>
      }
    >
      <div className="space-y-4">
        <p className="text-sm text-zinc-600 dark:text-zinc-400">
          The panel{"’"}s {geneId} sequence (top row) aligned base by base to where each strain
          matches it (bottom row). Highlighted letters are mismatches, dashes are insertions or
          deletions, and dotted rows are parts of the gene with no match in that strain.
        </p>
        <div className="grid grid-cols-1 xl:grid-cols-2 gap-4">
          <Side
            title={run.queries.find((q) => q.file_id === queryId)?.name ?? ""}
            view={left}
            error={leftError}
            loading={queryId !== null && !left && !leftError}
          />
          <RightSide
            run={run}
            geneId={geneId}
            leftId={queryId}
            rightId={rightId}
            onRight={setRightId}
            calls={calls}
            variant={left?.variant}
          />
        </div>
      </div>
    </Modal>
  );
}

function RightSide({
  run,
  geneId,
  leftId,
  rightId,
  onRight,
  calls,
  variant,
}: {
  run: Run;
  geneId: string;
  leftId: number | null;
  rightId: number | undefined;
  onRight: (id: number) => void;
  calls: Map<number, Call> | null;
  /** Align to the same panel sequence as the left side. */
  variant: string | undefined;
}) {
  const [view, setView] = useState<PanelAlignmentView | null>(null);
  const [error, setError] = useState<string | null>(null);

  useEffect(() => {
    if (rightId === undefined || variant === undefined) return;
    let cancelled = false;
    setView(null);
    setError(null);
    api
      .panelAlignment(run.id, rightId, geneId, variant)
      .then((v) => !cancelled && setView(v))
      .catch((e) => !cancelled && setError((e as Error).message));
    return () => {
      cancelled = true;
    };
  }, [run.id, rightId, geneId, variant]);

  const label = (c: Call | undefined) =>
    c === "PRESENT" ? "Present" : c === "PARTIAL" ? "Partial" : c === "ABSENT" ? "Absent" : "";

  return (
    <Side
      title={
        <label className="flex items-center gap-2">
          <span className="sr-only">Compare with</span>
          <select
            value={rightId ?? ""}
            onChange={(e) => onRight(Number(e.target.value))}
            className="h-9 px-2 rounded-lg border border-zinc-300 bg-white text-sm font-medium dark:border-zinc-700 dark:bg-zinc-900"
          >
            {run.queries
              .filter((q) => q.file_id !== leftId)
              .map((q) => (
                <option key={q.file_id} value={q.file_id}>
                  {q.name}
                  {calls?.get(q.file_id) ? ` (${label(calls.get(q.file_id))})` : ""}
                </option>
              ))}
          </select>
        </label>
      }
      view={view}
      error={error}
      loading={rightId !== undefined && !view && !error}
    />
  );
}

function Side({
  title,
  view,
  error,
  loading,
}: {
  title: React.ReactNode;
  view: PanelAlignmentView | null;
  error: string | null;
  loading: boolean;
}) {
  return (
    <div className="border border-zinc-200 rounded-lg overflow-hidden min-w-0 dark:border-zinc-800">
      <div className="flex flex-wrap items-center gap-x-3 gap-y-2 px-4 py-3 bg-zinc-50 border-b border-zinc-200 dark:bg-zinc-800/60 dark:border-zinc-800">
        <span className="font-medium truncate max-w-72">{title}</span>
        {view && <CallBadge call={view.call} />}
      </div>
      <div className="p-4 space-y-3">
        {error && <ErrorBox message={error} />}
        {loading && (
          <div className="flex items-center gap-3 text-zinc-500 py-6 justify-center dark:text-zinc-400">
            <Spinner /> Aligning...
          </div>
        )}
        {view && (
          <>
            <dl className="grid grid-cols-2 sm:grid-cols-4 gap-x-4 gap-y-2 text-sm">
              <Stat label="Gene covered">{view.panel_coverage.toFixed(0)} %</Stat>
              <Stat label="Identity">{view.identity.toFixed(1)} %</Stat>
              <Stat label="Mismatches">{view.mismatches.toLocaleString("en-US")}</Stat>
              <Stat label="Gap bases">{view.gap_bases.toLocaleString("en-US")}</Stat>
            </dl>
            <p className="text-xs text-zinc-500 font-mono dark:text-zinc-400">
              {view.contig}:{view.contig_start.toLocaleString("en-US")}-
              {view.contig_end.toLocaleString("en-US")} ({view.strand < 0 ? "-" : "+"} strand)
            </p>
            {view.match_note && (
              <p className="rounded-md border border-amber-200 bg-amber-50 px-3 py-2 text-sm text-amber-900 dark:border-amber-900 dark:bg-amber-950/40 dark:text-amber-200">
                {view.match_note}
              </p>
            )}
            <div className="space-y-4 thin-scroll max-h-[60vh] overflow-y-auto">
              <GeneAlignmentTiles
                refName="the panel gene"
                gene={{ reference_seq: view.panel_seq, start: 1, end: view.panel_len, strand: 1 }}
                q={{
                  query_id: view.query_id,
                  query_name: view.query_name,
                  call: view.call,
                  cov_pct: view.panel_coverage,
                  best_identity: view.identity,
                  mismatches: view.mismatches,
                  indels: view.gap_bases,
                  premature_stops: [],
                  blocks: [
                    {
                      ref_start: view.panel_start,
                      ref_end: view.panel_end,
                      qry_start: view.contig_start,
                      qry_end: view.contig_end,
                      qry_rev: view.strand < 0,
                      identity: view.identity,
                      ref_seq: view.panel_row,
                      qry_seq: view.strain_row,
                    },
                  ],
                  unaligned: [
                    ...(view.panel_start > 1 ? [[1, view.panel_start - 1] as [number, number]] : []),
                    ...(view.panel_end < view.panel_len
                      ? [[view.panel_end + 1, view.panel_len] as [number, number]]
                      : []),
                  ],
                }}
              />
            </div>
          </>
        )}
      </div>
    </div>
  );
}

function Stat({ label, children }: { label: string; children: React.ReactNode }) {
  return (
    <div>
      <dt className="text-xs text-zinc-500 dark:text-zinc-400">{label}</dt>
      <dd className="tabular-nums">{children}</dd>
    </div>
  );
}
