import { useEffect, useState } from "react";
import { api } from "../api";
import type { GeneDetail } from "../types";
import { nuccoreRangeUrl } from "../types";
import { CallBadge, Modal, Spinner } from "./ui";
import { AlignmentBlock } from "./GeneAlignmentPanel";

/**
 * The gene alignment viewer: pairwise alignment of the reference gene
 * against every query, with mismatch and indel highlighting, a coordinate
 * ruler and premature stop flags. Also exports FASTA / Clustal.
 */
export function GeneMsaDialog({
  runId,
  locus,
  onClose,
  onShowInGenome,
}: {
  runId: number;
  locus: string | null;
  onClose: () => void;
  onShowInGenome: (locus: string) => void;
}) {
  const [detail, setDetail] = useState<GeneDetail | null>(null);
  const [error, setError] = useState<string | null>(null);

  useEffect(() => {
    if (!locus) {
      setDetail(null);
      setError(null);
      return;
    }
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

  return (
    <Modal
      open={Boolean(locus)}
      onClose={onClose}
      wide
      title={
        <span className="flex items-baseline gap-3 flex-wrap">
          Gene alignment
          {detail && (
            <span className="text-sm font-normal text-zinc-500 dark:text-zinc-400">
              {detail.locus_tag}
              {detail.symbol ? ` (${detail.symbol})` : ""} - {detail.biotype} -{" "}
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
              {detail.protein_id && (
                <>
                  {" - "}
                  <a
                    href={`https://www.ncbi.nlm.nih.gov/protein/${detail.protein_id}`}
                    target="_blank"
                    rel="noreferrer"
                    className="underline hover:text-zinc-900 dark:hover:text-zinc-100"
                  >
                    {detail.protein_id}
                  </a>
                </>
              )}
            </span>
          )}
        </span>
      }
    >
      {error && <p className="text-red-700 dark:text-red-400">{error}</p>}
      {!detail && !error && (
        <div className="flex items-center gap-3 text-zinc-500 py-8 justify-center dark:text-zinc-400">
          <Spinner /> Preparing the alignment...
        </div>
      )}
      {detail && (
        <div className="space-y-6">
          <div className="flex flex-wrap gap-2">
            <a
              href={api.geneExportUrl(runId, detail.locus_tag, "fasta")}
              className="h-11 px-4 inline-flex items-center rounded-lg border border-zinc-300 hover:bg-zinc-100 text-[15px] dark:border-zinc-700 dark:hover:bg-zinc-800"
            >
              Export FASTA
            </a>
            <a
              href={api.geneExportUrl(runId, detail.locus_tag, "clustal")}
              className="h-11 px-4 inline-flex items-center rounded-lg border border-zinc-300 hover:bg-zinc-100 text-[15px] dark:border-zinc-700 dark:hover:bg-zinc-800"
            >
              Export Clustal
            </a>
            <button
              onClick={() => onShowInGenome(detail.locus_tag)}
              className="h-11 px-4 inline-flex items-center rounded-lg border border-zinc-300 hover:bg-zinc-100 text-[15px] dark:border-zinc-700 dark:hover:bg-zinc-800"
            >
              Show in genome view
            </button>
          </div>

          {detail.queries.length === 0 && (
            <p className="text-zinc-500 dark:text-zinc-400">
              This gene has no alignment data (the run may not have finished).
            </p>
          )}

          {detail.queries.map((q) => (
            <QueryAlignment key={q.query_id} q={q} />
          ))}

          <p className="text-xs text-zinc-400 dark:text-zinc-500">
            Reference row on top, query below. Highlighted letters are
            mismatches; dashes mark insertions or deletions. Blocks appear in
            reference order; unaligned stretches between blocks are listed
            below the blocks.
          </p>
        </div>
      )}
    </Modal>
  );
}

function QueryAlignment({ q }: { q: GeneDetail["queries"][number] }) {
  const [expanded, setExpanded] = useState(true);
  const alignedLen = q.blocks.reduce((a, b) => a + b.qry_seq.replace(/-/g, "").length, 0);
  const stats = `${q.cov_pct.toFixed(1)}% covered${
    q.best_identity > 0 ? `, ${q.best_identity.toFixed(1)}% identity` : ", no alignment"
  }, ${q.mismatches} mismatches, ${q.indels} indel bases`;
  return (
    <div className="border border-zinc-200 rounded-lg overflow-hidden dark:border-zinc-800">
      <div className="flex flex-wrap items-center justify-between gap-2 px-4 py-3 bg-zinc-50 border-b border-zinc-200 dark:bg-zinc-800/60 dark:border-zinc-800">
        <div className="flex items-center gap-3 min-w-0">
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
        </div>
        <div className="flex items-center gap-3">
          <span className="text-xs text-zinc-500 font-mono dark:text-zinc-400">{stats}</span>
          <button
            className="h-9 px-2 rounded-md text-sm text-zinc-500 hover:text-zinc-900 hover:bg-zinc-100 dark:text-zinc-400 dark:hover:text-zinc-100 dark:hover:bg-zinc-800"
            onClick={() => setExpanded(!expanded)}
          >
            {expanded ? "Hide" : "Show"}
          </button>
        </div>
      </div>
      {expanded && (
        <div className="p-4 space-y-4 thin-scroll max-h-96 overflow-y-auto">
          {q.premature_stops.length > 0 && (
            <p className="text-xs text-red-700 dark:text-red-400">
              Premature stop codons at amino acid position
              {q.premature_stops.length > 1 ? "s" : ""}{" "}
              {q.premature_stops.map((s) => s.aa_position).join(", ")} (of{" "}
              {Math.floor(alignedLen / 3)}).
            </p>
          )}
          {q.blocks.length === 0 && (
            <p className="text-sm text-zinc-400 dark:text-zinc-500">
              No part of this gene is aligned to this query.
            </p>
          )}
          {q.blocks.map((b, i) => (
            <AlignmentBlock key={i} block={b} />
          ))}
          {q.unaligned.length > 0 && (
            <div className="text-xs text-zinc-500 dark:text-zinc-400">
              {q.unaligned.map(([s, e], i) => (
                <p
                  key={i}
                  className="font-mono bg-zinc-50 border border-zinc-200 rounded px-2 py-1 my-1 inline-block mr-2 dark:bg-zinc-800/60 dark:border-zinc-800"
                >
                  unaligned reference bases {s.toLocaleString("en-US")} -{" "}
                  {e.toLocaleString("en-US")} ({(e - s + 1).toLocaleString("en-US")} bp)
                </p>
              ))}
            </div>
          )}
        </div>
      )}
    </div>
  );
}
