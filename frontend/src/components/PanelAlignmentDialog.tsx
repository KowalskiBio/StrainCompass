import { useEffect, useState } from "react";
import { api } from "../api";
import type { Call, PanelAlignmentView, PanelContext, Run } from "../types";
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
  // the run's presence thresholds, for "Present means ..." on each side
  const [limits, setLimits] = useState<{ pid: number; cov: number }>({ pid: 90, cov: 90 });
  useEffect(() => {
    api
      .getRunParams(run.id)
      .then((p) => setLimits({ pid: p.params.blast_pid, cov: p.params.blast_cov }))
      .catch(() => {});
  }, [run.id]);

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
        {left && (
          <SequenceActions
            label={`${geneId} (panel sequence, ${left.panel_len.toLocaleString("en-US")} bp)`}
            fastaHeader={`${geneId} ${left.variant_source || left.variant}`}
            seq={left.panel_seq}
            blastTitle={`Open NCBI BLAST with the panel's ${geneId} sequence`}
          />
        )}
        <p className="text-sm text-zinc-600 dark:text-zinc-400">
          The panel{"’"}s {geneId} sequence (top row) aligned base by base to where each strain
          matches it (bottom row). Highlighted letters are mismatches, dashes are insertions or
          deletions, and dotted rows are parts of the gene with no match in that strain.
        </p>
        <div className="grid grid-cols-1 xl:grid-cols-2 gap-4">
          <Side
            run={run}
            geneId={geneId}
            limits={limits}
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
            limits={limits}
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
  limits,
}: {
  run: Run;
  geneId: string;
  leftId: number | null;
  rightId: number | undefined;
  onRight: (id: number) => void;
  calls: Map<number, Call> | null;
  /** Align to the same panel sequence as the left side. */
  variant: string | undefined;
  limits: { pid: number; cov: number };
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
      run={run}
      geneId={geneId}
      limits={limits}
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
  run,
  geneId,
  limits,
  title,
  view,
  error,
  loading,
}: {
  run: Run;
  geneId: string;
  limits: { pid: number; cov: number };
  title: React.ReactNode;
  view: PanelAlignmentView | null;
  error: string | null;
  loading: boolean;
}) {
  // the strain's gene context names the related gene of a missing gene
  const [ctx, setCtx] = useState<PanelContext | null>(null);
  useEffect(() => {
    if (!view) return;
    let cancelled = false;
    setCtx(null);
    api
      .panelContext(run.id, view.query_id, geneId)
      .then((c) => !cancelled && setCtx(c))
      .catch(() => {});
    return () => {
      cancelled = true;
    };
  }, [run.id, view, geneId]);
  const strain = view?.query_name.replace(/\.(fasta|fa|fna)$/i, "") ?? "";
  const related = view && view.call !== "PRESENT";
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
            <Verdict view={view} geneId={geneId} ctx={ctx} />
            {related && (
              <p className="text-xs font-medium text-zinc-500 dark:text-zinc-400">
                How the related gene compares with {geneId}:
              </p>
            )}
            <dl className="grid grid-cols-2 sm:grid-cols-4 gap-x-4 gap-y-2 text-sm">
              <Stat label={`Share of ${geneId} matched`}>{view.panel_coverage.toFixed(0)} %</Stat>
              <Stat label="DNA identity">{view.identity.toFixed(1)} %</Stat>
              <Stat label="Mismatches">{view.mismatches.toLocaleString("en-US")}</Stat>
              <Stat label="Gap bases">{view.gap_bases.toLocaleString("en-US")}</Stat>
            </dl>
            <p className="text-xs text-zinc-500 dark:text-zinc-400">
              Present means at least {limits.pid} % identity over at least {limits.cov} % of the
              gene.
            </p>
            <p className="text-xs text-zinc-500 font-mono dark:text-zinc-400">
              {view.contig}:{view.contig_start.toLocaleString("en-US")}-
              {view.contig_end.toLocaleString("en-US")} ({view.strand < 0 ? "-" : "+"} strand)
            </p>
            <SequenceActions
              label="Matched stretch of this strain"
              fastaHeader={`${view.query_name.replace(/\.(fasta|fa|fna)$/i, "")} ${view.contig}:${view.contig_start}-${view.contig_end}(${view.strand < 0 ? "-" : "+"}) best match to ${view.gene_id}`}
              seq={view.strain_row.replace(/-/g, "")}
              blastTitle="Open NCBI BLAST with this stretch of your strain (it is sent to NCBI)"
              compact
            />
            {view.match_note && (
              <p className="rounded-md border border-amber-200 bg-amber-50 px-3 py-2 text-sm text-amber-900 dark:border-amber-900 dark:bg-amber-950/40 dark:text-amber-200">
                {view.match_note}
              </p>
            )}
            <p className="text-xs text-zinc-600 dark:text-zinc-400">
              Top row: <b>{geneId}</b> ({view.variant_source || "panel sequence"}). Bottom row:{" "}
              <b>{strain}</b>
              {related ? `, the related gene${ctx?.related_gene ? ` ${ctx.related_gene}` : ""}` : ""}.
            </p>
            <div className="space-y-4 thin-scroll max-h-[60vh] overflow-y-auto">
              <GeneAlignmentTiles
                refName={geneId}
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

/** One plain sentence: is the gene of interest in this strain? */
function Verdict({
  view,
  geneId,
  ctx,
}: {
  view: PanelAlignmentView;
  geneId: string;
  ctx: PanelContext | null;
}) {
  const cov = view.panel_coverage.toFixed(0);
  const id = view.identity.toFixed(0);
  if (view.call === "PRESENT")
    return (
      <p className="rounded-md border border-emerald-200 bg-emerald-50 px-3 py-2 text-sm text-emerald-900 dark:border-emerald-900 dark:bg-emerald-950/40 dark:text-emerald-200">
        <b>{geneId} is present</b>: {cov} % of the gene, {id} % identical.
      </p>
    );
  if (view.call === "PARTIAL")
    return (
      <p className="rounded-md border border-amber-200 bg-amber-50 px-3 py-2 text-sm text-amber-900 dark:border-amber-900 dark:bg-amber-950/40 dark:text-amber-200">
        <b>{geneId} is partly present</b>: {cov} % of the gene matches, {id} % identical.
      </p>
    );
  const name = ctx?.related_gene;
  const prot = ctx?.related_identity;
  return (
    <p className="rounded-md border border-zinc-200 bg-zinc-50 px-3 py-2 text-sm text-zinc-800 dark:border-zinc-700 dark:bg-zinc-800/60 dark:text-zinc-200">
      <b>{geneId} is absent.</b> Shown below is the closest related gene
      {name ? <>, <i>{name}</i></> : ""}
      {prot != null ? `, ${prot.toFixed(0)} % identical to ${geneId} as protein` : ""}. It is a
      different gene, not a copy of {geneId}.
    </p>
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

/** NCBI BLAST's web form, prefilled with a nucleotide query. */
function blastUrl(seq: string): string {
  return `https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome&QUERY=${encodeURIComponent(seq)}`;
}

async function copyText(text: string): Promise<boolean> {
  try {
    await navigator.clipboard.writeText(text);
    return true;
  } catch {
    // no clipboard API (plain http): the old selection route
    const ta = document.createElement("textarea");
    ta.value = text;
    ta.style.position = "fixed";
    ta.style.opacity = "0";
    document.body.appendChild(ta);
    ta.select();
    const ok = document.execCommand("copy");
    ta.remove();
    return ok;
  }
}

/** Copy a sequence (plain or FASTA) or send it to NCBI BLAST. */
function SequenceActions({
  label,
  fastaHeader,
  seq,
  blastTitle,
  compact,
}: {
  label: string;
  fastaHeader: string;
  seq: string;
  blastTitle: string;
  compact?: boolean;
}) {
  const [copied, setCopied] = useState<"seq" | "fasta" | null>(null);
  const fasta = `>${fastaHeader}\n${seq.match(/.{1,60}/g)?.join("\n") ?? ""}\n`;
  const copy = (what: "seq" | "fasta") =>
    copyText(what === "seq" ? seq : fasta).then((ok) => {
      if (!ok) return;
      setCopied(what);
      setTimeout(() => setCopied(null), 1500);
    });
  const btn = `${compact ? "h-8 px-2.5 text-xs" : "h-9 px-3 text-sm"} inline-flex items-center rounded-lg border border-zinc-300 bg-white font-medium hover:bg-zinc-100 dark:border-zinc-700 dark:bg-zinc-900 dark:hover:bg-zinc-800`;
  return (
    <div className="flex flex-wrap items-center gap-2">
      <span className={`${compact ? "text-xs" : "text-sm"} text-zinc-600 dark:text-zinc-400 mr-1`}>
        {label}
      </span>
      <button className={btn} onClick={() => copy("seq")} title="Copy the bases only">
        {copied === "seq" ? "Copied" : "Copy sequence"}
      </button>
      <button className={btn} onClick={() => copy("fasta")} title="Copy as FASTA, with a header line">
        {copied === "fasta" ? "Copied" : "Copy FASTA"}
      </button>
      <a className={btn} href={blastUrl(seq)} target="_blank" rel="noreferrer" title={blastTitle}>
        BLAST at NCBI
      </a>
    </div>
  );
}
