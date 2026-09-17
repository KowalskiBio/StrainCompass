import { useEffect, useRef, useState } from "react";
import { api } from "../api";
import type { GeneDetail, GeneQueryAlignment, WgaGene } from "../types";
import { nuccoreRangeUrl } from "../types";
import { CallBadge, Spinner } from "./ui";

/**
 * The alignment of one clicked map gene to the query the strain map is
 * colored by, shown beneath the map: reference row on top, query below,
 * mismatches and indels highlighted. The gene dialog (every query at
 * once) stays one click away.
 */
export function GeneAlignmentPanel({
  runId,
  gene,
  queryId,
  onClose,
  onOpenAll,
}: {
  runId: number;
  /** The clicked map gene: supplies the function annotation until the
   * detail payload arrives (it carries sequences, not the product). */
  gene: WgaGene;
  /** The query whose alignment is shown; the one coloring the map. */
  queryId: number | undefined;
  onClose: () => void;
  onOpenAll: (locus: string) => void;
}) {
  const locus = gene.locus_tag;
  const [detail, setDetail] = useState<GeneDetail | null>(null);
  const [error, setError] = useState<string | null>(null);
  const ref = useRef<HTMLDivElement>(null);

  useEffect(() => {
    let cancelled = false;
    setDetail(null);
    setError(null);
    api
      .geneDetail(runId, locus)
      .then((d) => {
        if (!cancelled) setDetail(d);
      })
      .catch((e) => {
        if (!cancelled) setError((e as Error).message);
      });
    return () => {
      cancelled = true;
    };
  }, [runId, locus]);

  // The map fills the window, so the panel usually opens below the fold.
  useEffect(() => {
    ref.current?.scrollIntoView({ behavior: "smooth", block: "nearest" });
  }, [locus]);

  const q =
    detail?.queries.find((x) => x.query_id === queryId) ??
    (detail && detail.queries.length > 0 ? detail.queries[0] : null);

  return (
    <div
      ref={ref}
      className="border border-zinc-200 rounded-xl bg-white dark:border-zinc-800 dark:bg-zinc-900"
    >
      {/* header: gene on the left, actions on the right */}
      <div className="flex flex-wrap items-start justify-between gap-3 px-4 py-3 border-b border-zinc-200 dark:border-zinc-800">
        {detail ? (
          <div className="min-w-0">
            <p className="font-semibold text-[15px]">
              {detail.locus_tag}
              {detail.symbol && (
                <span className="font-normal text-zinc-500 dark:text-zinc-400">
                  {" "}
                  ({detail.symbol})
                </span>
              )}
            </p>
            <p className="text-xs text-zinc-500 mt-0.5 dark:text-zinc-400">
              {detail.biotype} -{" "}
              <a
                href={nuccoreRangeUrl(detail.seqid, detail.start, detail.end)}
                target="_blank"
                rel="noreferrer"
                className="underline hover:text-zinc-900 dark:hover:text-zinc-100"
              >
                {detail.seqid}:{detail.start.toLocaleString("en-US")}-
                {detail.end.toLocaleString("en-US")}
              </a>{" "}
              ({detail.strand > 0 ? "+" : "-"} strand)
            </p>
            {gene.product && (
              <p className="text-xs text-zinc-600 mt-1 leading-snug dark:text-zinc-300">
                {gene.product}
              </p>
            )}          </div>
        ) : (
          <p className="font-semibold text-[15px]">{locus}</p>
        )}
        <div className="flex items-center gap-2 shrink-0">
          {detail && (
            <>
              <a
                href={api.geneExportUrl(runId, detail.locus_tag, "fasta")}
                className="h-9 px-3 inline-flex items-center rounded-lg border border-zinc-300 hover:bg-zinc-100 text-sm dark:border-zinc-700 dark:hover:bg-zinc-800"
              >
                FASTA
              </a>
              <a
                href={api.geneExportUrl(runId, detail.locus_tag, "clustal")}
                className="h-9 px-3 inline-flex items-center rounded-lg border border-zinc-300 hover:bg-zinc-100 text-sm dark:border-zinc-700 dark:hover:bg-zinc-800"
              >
                Clustal
              </a>
              <button
                onClick={() => onOpenAll(detail.locus_tag)}
                className="h-9 px-3 rounded-lg border border-zinc-300 hover:bg-zinc-100 text-sm dark:border-zinc-700 dark:hover:bg-zinc-800"
              >
                All queries
              </button>
            </>
          )}
          <button
            onClick={onClose}
            aria-label="Close"
            className="w-9 h-9 grid place-items-center rounded-lg text-zinc-400 hover:text-zinc-900 hover:bg-zinc-100 dark:text-zinc-500 dark:hover:text-zinc-100 dark:hover:bg-zinc-800"
          >
            <svg width="12" height="12" viewBox="0 0 16 16" fill="none">
              <path
                d="M3 3l10 10M13 3L3 13"
                stroke="currentColor"
                strokeWidth="1.8"
                strokeLinecap="round"
              />
            </svg>
          </button>
        </div>
      </div>

      {/* body: the query's alignment */}
      <div className="p-4">
        {error && <p className="text-red-700 dark:text-red-400">{error}</p>}
        {!detail && !error && (
          <div className="flex items-center gap-3 text-zinc-500 py-6 justify-center dark:text-zinc-400">
            <Spinner /> Preparing the alignment...
          </div>
        )}
        {detail && !q && (
          <p className="text-sm text-zinc-500 dark:text-zinc-400">
            This gene has no alignment data.
          </p>
        )}
        {detail && q && <QueryAlignment q={q} gene={detail} />}
        {detail && (
          <p className="text-xs text-zinc-400 mt-3 dark:text-zinc-500">
            Reference row on top, query below. Highlighted letters are
            mismatches; dashes mark insertions or deletions; a dotted query
            row marks reference stretches with no alignment to this query.
          </p>
        )}
      </div>
    </div>
  );
}

/** The clicked gene against one query: stats, then the alignment blocks. */
function QueryAlignment({ q, gene }: { q: GeneQueryAlignment; gene: TileGene }) {
  const alignedLen = q.blocks.reduce((a, b) => a + b.qry_seq.replace(/-/g, "").length, 0);
  return (
    <div className="space-y-4">
      <div className="flex flex-wrap items-center gap-3">
        <span className="font-medium truncate max-w-64" title={q.query_name}>
          {q.query_name}
        </span>
        <CallBadge call={q.call} />
        {q.premature_stops.length > 0 && (
          <span className="text-xs text-red-700 bg-red-50 border border-red-200 rounded-full px-2 py-0.5 font-medium dark:text-red-400 dark:bg-red-950/40 dark:border-red-900">
            {q.premature_stops.length} premature stop
            {q.premature_stops.length > 1 ? "s" : ""}
          </span>
        )}
        <span className="text-xs text-zinc-500 font-mono dark:text-zinc-400">
          {q.cov_pct.toFixed(1)}% covered
          {q.best_identity > 0
            ? `, ${q.best_identity.toFixed(1)}% identity`
            : ", no alignment"}
          , {q.mismatches} mismatches, {q.indels} indel bases
        </span>
      </div>
      {q.premature_stops.length > 0 && (
        <p className="text-xs text-red-700 dark:text-red-400">
          Premature stop codons at amino acid position
          {q.premature_stops.length > 1 ? "s" : ""}{" "}
          {q.premature_stops.map((s) => s.aa_position).join(", ")} (of{" "}
          {Math.floor(alignedLen / 3)}).
        </p>
      )}
      {q.blocks.length === 0 && q.unaligned.length === 0 && (
        <p className="text-sm text-zinc-400 dark:text-zinc-500">
          No part of this gene is aligned to this query.
        </p>
      )}
      <div className="space-y-4 max-h-96 overflow-y-auto thin-scroll">
        <GeneAlignmentTiles q={q} gene={gene} />
      </div>
    </div>
  );
}

const COLS = 60;

/** The gene facts the tiling needs, from GeneDetail. */
export type TileGene = Pick<GeneDetail, "reference_seq" | "start" | "end" | "strand">;

/** One piece of the gene walk: an aligned block, or a stretch the query
 * has no alignment to. */
type Tile =
  | { kind: "block"; block: GeneQueryAlignment["blocks"][number] }
  | { kind: "unaligned"; start: number; end: number; seq: string };

/**
 * Interleave the aligned blocks and the unaligned stretches of one query
 * into a single walk of the gene, so a partial gene reads as one picture:
 * what aligned, then what did not, in order.
 *
 * The gene's own reference sequence supplies the unaligned bases (the
 * blocks only carry what aligned), indexed through the gene-oriented
 * `reference_seq`: a minus-strand gene's sequence is the reverse complement,
 * so its tiles also run in descending reference order to read 5' -> 3'.
 */
export function GeneAlignmentTiles({ q, gene }: { q: GeneQueryAlignment; gene: TileGene }) {
  const tiles: Tile[] = q.blocks.map((b) => ({ kind: "block" as const, block: b }));
  for (const [s, e] of q.unaligned) {
    const from = gene.strand < 0 ? gene.end - e : s - gene.start;
    tiles.push({
      kind: "unaligned",
      start: s,
      end: e,
      seq: gene.reference_seq.slice(from, from + (e - s + 1)),
    });
  }
  tiles.sort((a, b) =>
    (a.kind === "block" ? a.block.ref_start : a.start) -
    (b.kind === "block" ? b.block.ref_start : b.start),
  );
  const ordered = gene.strand < 0 ? tiles.reverse() : tiles;
  return (
    <>
      {ordered.map((t, i) =>
        t.kind === "block" ? (
          <AlignmentBlock key={i} block={t.block} />
        ) : (
          <UnalignedBlock key={i} start={t.start} end={t.end} seq={t.seq} />
        ),
      )}
    </>
  );
}

/** A stretch of the reference with no alignment to the query: the
 * reference bases shown, the query row dotted out. */
function UnalignedBlock({
  start,
  end,
  seq,
}: {
  start: number;
  end: number;
  seq: string;
}) {
  return (
    <div>
      <p className="text-xs text-zinc-400 mb-1 font-mono dark:text-zinc-500">
        unaligned reference {start.toLocaleString("en-US")} -{" "}
        {end.toLocaleString("en-US")} ({(end - start + 1).toLocaleString("en-US")} bp),
        no alignment to this query
      </p>
      <div className="font-mono text-xs leading-5 overflow-x-auto thin-scroll">
        {chunk(seq.length, COLS).map(([, colStart]) => (
          <div key={colStart} className="whitespace-pre">
            <span className="text-zinc-300 select-none inline-block w-16 text-right pr-2 dark:text-zinc-700">
              {colStart + 1}
            </span>
            <span className="text-zinc-500 dark:text-zinc-500">
              {seq.slice(colStart, colStart + COLS)}
            </span>
            {"\n"}
            <span className="text-zinc-300 select-none inline-block w-16 text-right pr-2 dark:text-zinc-700">
              {" "}
            </span>
            <span className="text-zinc-300 select-none dark:text-zinc-700">
              {"\u00b7".repeat(Math.min(COLS, seq.length - colStart))}
            </span>
          </div>
        ))}
      </div>
    </div>
  );
}

/** One aligned block: the shared rendering the gene dialog also uses. */
export function AlignmentBlock({
  block,
}: {
  block: GeneQueryAlignment["blocks"][number];
}) {
  const ref = block.ref_seq;
  const qry = block.qry_seq;
  const nCols = Math.min(ref.length, qry.length);
  return (
    <div>
      <p className="text-xs text-zinc-500 mb-1 font-mono dark:text-zinc-400">
        block {block.ref_start.toLocaleString("en-US")} -{" "}
        {block.ref_end.toLocaleString("en-US")} in the reference
        {block.qry_rev ? ", query aligned on the reverse strand" : ""}, identity{" "}
        {block.identity.toFixed(1)}%
      </p>
      <div className="font-mono text-xs leading-5 overflow-x-auto thin-scroll">
        {chunk(nCols, COLS).map(([, colStart]) => (
          <div key={colStart} className="whitespace-pre">
            <span className="text-zinc-300 select-none inline-block w-16 text-right pr-2 dark:text-zinc-700">
              {colStart + 1}
            </span>
            <Row seq={ref.slice(colStart, colStart + COLS)} kind="ref" other={qry.slice(colStart, colStart + COLS)} />
            {"\n"}
            <span className="text-zinc-300 select-none inline-block w-16 text-right pr-2 dark:text-zinc-700">
              {" "}
            </span>
            <Row seq={qry.slice(colStart, colStart + COLS)} kind="qry" other={ref.slice(colStart, colStart + COLS)} />
          </div>
        ))}
      </div>
      <p className="text-xs text-zinc-400 mt-1 dark:text-zinc-500">{nCols} alignment columns</p>
    </div>
  );
}

/**
 * One monospace row. Mismatches are tinted red; gaps (dashes) amber.
 * The reference row is plain so the eye is drawn to query differences.
 */
function Row({ seq, kind, other }: { seq: string; kind: "ref" | "qry"; other: string }) {
  const out: React.ReactNode[] = [];
  for (let i = 0; i < seq.length; i++) {
    const c = seq[i];
    const o = other[i];
    let cls = "";
    if (c === "-" || o === "-") {
      cls = kind === "qry" ? "bg-amber-100 text-amber-900 dark:bg-amber-950/50 dark:text-amber-300" : "";
    } else if (kind === "qry" && c !== o) {
      cls = "bg-red-100 text-red-800 dark:bg-red-950/50 dark:text-red-300";
    }
    out.push(
      cls ? (
        <span key={i} className={cls}>
          {c}
        </span>
      ) : (
        <span key={i}>{c}</span>
      ),
    );
  }
  return <span className={kind === "qry" ? "text-zinc-800 dark:text-zinc-300" : "text-zinc-500 dark:text-zinc-500"}>{out}</span>;
}

function chunk(total: number, n: number): [number, number][] {
  const out: [number, number][] = [];
  for (let i = 0; i < total; i += n) out.push([i, i]);
  return out;
}
