import { useEffect, useMemo, useState } from "react";
import { api } from "../api";
import type { ElementReport, PanelContext, Run } from "../types";
import { CallBadge, ErrorBox, Modal, Spinner } from "./ui";

/**
 * A panel gene beyond its call: where it sits in one strain (plasmid,
 * chromosomal insertion, mobile-element neighbours), and whether the
 * whole element carrying it is missing from the strains that lack it.
 */
export function PanelGeneDialog({
  run,
  geneId,
  initialQueryId,
  onClose,
}: {
  run: Run;
  geneId: string | null;
  initialQueryId: number | undefined;
  onClose: () => void;
}) {
  const [queryId, setQueryId] = useState<number | undefined>(initialQueryId);
  useEffect(() => setQueryId(initialQueryId), [geneId, initialQueryId]);

  return (
    <Modal
      open={Boolean(geneId)}
      onClose={onClose}
      wide
      title={
        <span className="flex items-baseline gap-3 flex-wrap">
          Panel gene <span className="font-mono">{geneId}</span>
        </span>
      }
    >
      {geneId && (
        <div className="space-y-8">
          <WhereItSits
            run={run}
            geneId={geneId}
            queryId={queryId ?? run.queries[0]?.file_id}
            onQuery={setQueryId}
          />
          {run.queries.length > 1 && (
            <AcrossStrains run={run} geneId={geneId} preferredSource={queryId} />
          )}
        </div>
      )}
    </Modal>
  );
}

function fmtBp(n: number): string {
  if (n >= 1_000_000) return `${(n / 1e6).toFixed(2)} Mb`;
  if (n >= 10_000) return `${(n / 1e3).toFixed(0)} kb`;
  if (n >= 1_000) return `${(n / 1e3).toFixed(1)} kb`;
  return `${n} bp`;
}

const VERDICT_STYLE: Record<PanelContext["verdict"], string> = {
  plasmid:
    "border-violet-200 bg-violet-50 text-violet-900 dark:border-violet-900 dark:bg-violet-950/40 dark:text-violet-200",
  chromosome_insertion:
    "border-teal-200 bg-teal-50 text-teal-900 dark:border-teal-900 dark:bg-teal-950/40 dark:text-teal-200",
  chromosome_shared:
    "border-zinc-200 bg-zinc-50 text-zinc-800 dark:border-zinc-800 dark:bg-zinc-800/60 dark:text-zinc-200",
  unplaced:
    "border-amber-200 bg-amber-50 text-amber-900 dark:border-amber-900 dark:bg-amber-950/40 dark:text-amber-200",
  not_found:
    "border-zinc-200 bg-zinc-50 text-zinc-600 dark:border-zinc-800 dark:bg-zinc-800/60 dark:text-zinc-400",
};

const VERDICT_LABEL: Record<PanelContext["verdict"], string> = {
  plasmid: "Likely plasmid",
  chromosome_insertion: "Chromosomal insertion",
  chromosome_shared: "Chromosome",
  unplaced: "Unplaced contig",
  not_found: "Not found",
};

function SectionTitle({ children }: { children: React.ReactNode }) {
  return (
    <h3 className="text-base font-semibold text-zinc-900 dark:text-zinc-100">{children}</h3>
  );
}

function QuerySelect({
  run,
  value,
  onChange,
  only,
}: {
  run: Run;
  value: number | undefined;
  onChange: (id: number) => void;
  only?: Set<number>;
}) {
  return (
    <select
      value={value ?? ""}
      onChange={(e) => onChange(Number(e.target.value))}
      className="h-9 px-2 rounded-lg border border-zinc-300 bg-white text-sm dark:border-zinc-700 dark:bg-zinc-900"
    >
      {run.queries
        .filter((q) => !only || only.has(q.file_id))
        .map((q) => (
          <option key={q.file_id} value={q.file_id}>
            {q.name}
          </option>
        ))}
    </select>
  );
}

function WhereItSits({
  run,
  geneId,
  queryId,
  onQuery,
}: {
  run: Run;
  geneId: string;
  queryId: number | undefined;
  onQuery: (id: number) => void;
}) {
  const [ctx, setCtx] = useState<PanelContext | null>(null);
  const [error, setError] = useState<string | null>(null);

  useEffect(() => {
    if (queryId === undefined) return;
    let cancelled = false;
    setCtx(null);
    setError(null);
    api
      .panelContext(run.id, queryId, geneId)
      .then((c) => !cancelled && setCtx(c))
      .catch((e) => !cancelled && setError((e as Error).message));
    return () => {
      cancelled = true;
    };
  }, [run.id, queryId, geneId]);

  const alignedPct =
    ctx?.contig && ctx.contig.length > 0
      ? (100 * ctx.contig.aligned_bp) / ctx.contig.length
      : 0;

  return (
    <section className="space-y-3">
      <div className="flex items-center justify-between gap-3 flex-wrap">
        <SectionTitle>Where does it sit?</SectionTitle>
        {run.queries.length > 1 && (
          <label className="flex items-center gap-2 text-sm text-zinc-600 dark:text-zinc-400">
            Strain
            <QuerySelect run={run} value={queryId} onChange={onQuery} />
          </label>
        )}
      </div>
      {error && <ErrorBox message={error} />}
      {!ctx && !error && <Spinner className="text-zinc-400" />}
      {ctx && (
        <>
          <div className={`rounded-lg border px-4 py-3 text-[15px] ${VERDICT_STYLE[ctx.verdict]}`}>
            <p className="font-semibold">{VERDICT_LABEL[ctx.verdict]}</p>
            <p className="mt-1">{ctx.verdict_text}</p>
          </div>

          {ctx.contig && (
            <dl className="grid grid-cols-2 sm:grid-cols-4 gap-x-6 gap-y-3 text-sm">
              <Fact label="Call">
                <CallBadge call={ctx.call} />
                <span className="ml-2 text-zinc-500 dark:text-zinc-400">
                  {ctx.cov_pct.toFixed(0)} % cov., {ctx.identity.toFixed(1)} % id.
                </span>
              </Fact>
              <Fact label="Contig">
                <span className="font-mono">{ctx.contig.seqid}</span>
                <span className="ml-1 text-zinc-500 dark:text-zinc-400">({fmtBp(ctx.contig.length)})</span>
              </Fact>
              <Fact label="Hit">
                <span className="font-mono tabular-nums">
                  {ctx.hit_start.toLocaleString("en-US")}-{ctx.hit_end.toLocaleString("en-US")} (
                  {ctx.hit_strand < 0 ? "-" : "+"})
                </span>
              </Fact>
              <Fact label="Contig aligned to reference">
                <AlignedBar pct={alignedPct} />
              </Fact>
              {ctx.region && (
                <Fact label="Gained region (no reference counterpart)">
                  <span className="font-mono tabular-nums">
                    {ctx.region.start.toLocaleString("en-US")}-{ctx.region.end.toLocaleString("en-US")}
                  </span>
                  <span className="ml-1 text-zinc-500 dark:text-zinc-400">({fmtBp(ctx.region.length)})</span>
                </Fact>
              )}
              {ctx.region && ctx.region.anchor !== "unanchored" && (
                <Fact label="Flanking reference genes">
                  {ctx.region.left_gene || "-"} / {ctx.region.right_gene || "-"}
                </Fact>
              )}
              <Fact label="Genome">
                {fmtBp(ctx.genome_bp)}, {ctx.n_contigs} contigs
              </Fact>
            </dl>
          )}

          {ctx.contig && ctx.region && (
            <div>
              <p className="text-sm font-medium text-zinc-700 dark:text-zinc-300 mb-2">
                Genes within {fmtBp(ctx.window)} of the hit
              </p>
              {ctx.genes.length > 0 ? (
                <ContextGenes ctx={ctx} />
              ) : (
                <p className="text-sm text-zinc-500 dark:text-zinc-400">
                  {ctx.genes_note ||
                    "No genes were predicted around the hit. Genes are predicted and named only inside gained regions."}
                </p>
              )}
            </div>
          )}
        </>
      )}
    </section>
  );
}

function Fact({ label, children }: { label: string; children: React.ReactNode }) {
  return (
    <div>
      <dt className="text-xs text-zinc-500 dark:text-zinc-400">{label}</dt>
      <dd className="mt-0.5 flex items-center flex-wrap text-zinc-900 dark:text-zinc-100">{children}</dd>
    </div>
  );
}

function AlignedBar({ pct }: { pct: number }) {
  return (
    <span className="inline-flex items-center gap-2 whitespace-nowrap">
      <span className="w-20 h-2 rounded-full bg-zinc-200 overflow-hidden dark:bg-zinc-700">
        <span
          className="block h-full bg-zinc-500 dark:bg-zinc-400"
          style={{ width: `${Math.min(100, pct)}%` }}
        />
      </span>
      <span className="tabular-nums">{pct.toFixed(0)} %</span>
    </span>
  );
}

function ContextGenes({ ctx }: { ctx: PanelContext }) {
  return (
    <div className="border border-zinc-200 rounded-lg overflow-hidden dark:border-zinc-800">
      <table className="w-full text-sm">
        <thead className="bg-zinc-50 text-zinc-500 text-xs dark:bg-zinc-800/60 dark:text-zinc-400">
          <tr>
            <th className="text-left font-medium px-3 py-2">Position</th>
            <th className="text-left font-medium px-3 py-2">Name</th>
            <th className="text-right font-medium px-3 py-2">Distance</th>
          </tr>
        </thead>
        <tbody>
          {ctx.genes.map((g) => (
            <tr
              key={`${g.start}-${g.end}`}
              className={`border-t border-zinc-100 dark:border-zinc-800 ${
                g.is_hit ? "bg-emerald-50/60 dark:bg-emerald-950/20" : ""
              }`}
            >
              <td className="px-3 py-1.5 font-mono tabular-nums whitespace-nowrap">
                {g.start.toLocaleString("en-US")}-{g.end.toLocaleString("en-US")} ({g.strand < 0 ? "-" : "+"})
              </td>
              <td className="px-3 py-1.5">
                <span className="inline-flex items-center gap-2 flex-wrap">
                  {g.label || <span className="text-zinc-400 dark:text-zinc-500">unnamed</span>}
                  {g.is_hit && (
                    <span className="px-1.5 py-0.5 rounded text-xs bg-emerald-100 text-emerald-800 dark:bg-emerald-900/40 dark:text-emerald-300">
                      {ctx.gene_id}
                    </span>
                  )}
                  {g.mobile && (
                    <span className="px-1.5 py-0.5 rounded text-xs bg-violet-100 text-violet-800 dark:bg-violet-900/40 dark:text-violet-300">
                      mobile element
                    </span>
                  )}
                  {g.partial && (
                    <span className="text-xs text-zinc-400 dark:text-zinc-500">partial</span>
                  )}
                </span>
              </td>
              <td className="px-3 py-1.5 text-right tabular-nums text-zinc-500 dark:text-zinc-400">
                {g.is_hit ? "" : fmtBp(g.distance)}
              </td>
            </tr>
          ))}
        </tbody>
      </table>
    </div>
  );
}

function AcrossStrains({
  run,
  geneId,
  preferredSource,
}: {
  run: Run;
  geneId: string;
  preferredSource: number | undefined;
}) {
  // the strains carrying the gene, from the panel matrix row
  const [positives, setPositives] = useState<Set<number> | null>(null);
  const [source, setSource] = useState<number | undefined>(undefined);
  const [useAccession, setUseAccession] = useState(false);
  const [accession, setAccession] = useState("");
  const [report, setReport] = useState<ElementReport | null>(null);
  const [busy, setBusy] = useState(false);
  const [error, setError] = useState<string | null>(null);

  useEffect(() => {
    let cancelled = false;
    setPositives(null);
    setReport(null);
    setError(null);
    api
      .panelMatrix(run.id, { search: geneId, page: 0, page_size: 1000 })
      .then((p) => {
        if (cancelled) return;
        const row = p.rows.find((r) => r.gene_id === geneId);
        const pos = new Set<number>();
        row?.calls.forEach((c, i) => {
          if (c !== "ABSENT") pos.add(run.queries[i].file_id);
        });
        setPositives(pos);
      })
      .catch((e) => !cancelled && setError((e as Error).message));
    return () => {
      cancelled = true;
    };
  }, [run, geneId]);

  useEffect(() => {
    if (!positives) return;
    setSource(
      preferredSource !== undefined && positives.has(preferredSource)
        ? preferredSource
        : [...positives][0],
    );
  }, [positives, preferredSource]);

  async function compare() {
    if (source === undefined) return;
    setBusy(true);
    setError(null);
    setReport(null);
    try {
      setReport(
        await api.panelElement(run.id, geneId, source, useAccession ? accession.trim() : undefined),
      );
    } catch (e) {
      setError((e as Error).message);
    } finally {
      setBusy(false);
    }
  }

  const sourceSize = useMemo(
    () => report?.hits.find((h) => h.query_id === report.source_query_id)?.genome_bp,
    [report],
  );
  const noPositives = positives !== null && positives.size === 0;

  return (
    <section className="space-y-3">
      <SectionTitle>Is the whole element missing in the other strains?</SectionTitle>
      <p className="text-sm text-zinc-600 dark:text-zinc-400">
        Takes the DNA carrying {geneId} in one strain and measures how much of it every strain of
        this run contains. A strain lacking the whole element (e.g. a plasmid) differs from one
        lacking only the gene.
      </p>
      {noPositives ? (
        <p className="text-sm text-zinc-500 dark:text-zinc-400">
          {geneId} was not found in any strain of this run, so there is no element to compare. You
          can still compare a complete record from NCBI:
        </p>
      ) : null}
      <div className="flex flex-wrap items-end gap-4 text-sm">
        {!noPositives && (
          <fieldset className="space-y-2">
            <label className="flex items-center gap-2">
              <input
                type="radio"
                checked={!useAccession}
                onChange={() => setUseAccession(false)}
              />
              Contig carrying the gene in
              {positives ? (
                <QuerySelect run={run} value={source} onChange={setSource} only={positives} />
              ) : (
                <Spinner className="text-zinc-400" />
              )}
            </label>
            <label className="flex items-center gap-2">
              <input type="radio" checked={useAccession} onChange={() => setUseAccession(true)} />
              Complete record from NCBI (e.g. a plasmid)
            </label>
          </fieldset>
        )}
        {(useAccession || noPositives) && (
          <input
            value={accession}
            onChange={(e) => {
              setAccession(e.target.value);
              setUseAccession(true);
            }}
            placeholder="NZ_CP168866.1"
            className="h-9 px-3 rounded-lg border border-zinc-300 font-mono text-sm w-52 dark:border-zinc-700 dark:bg-zinc-900"
          />
        )}
        <button
          onClick={compare}
          disabled={
            busy ||
            (useAccession || noPositives
              ? !accession.trim()
              : source === undefined)
          }
          className="h-9 px-4 rounded-lg bg-zinc-900 text-white font-medium hover:bg-zinc-700 disabled:opacity-40 dark:bg-zinc-100 dark:text-zinc-900 dark:hover:bg-zinc-300"
        >
          {busy ? "Comparing..." : "Compare"}
        </button>
        {busy && (
          <span className="text-zinc-500 dark:text-zinc-400">
            Searching every strain; this takes a few seconds per strain.
          </span>
        )}
      </div>
      {error && <ErrorBox message={error} />}
      {report && (
        <>
          <p className="text-sm text-zinc-700 dark:text-zinc-300">
            Element:{" "}
            <span className="font-mono">{report.element_name}</span>
            {report.element_title && <> - {report.element_title}</>} ({fmtBp(report.element_len)})
          </p>
          <div className="border border-zinc-200 rounded-lg overflow-x-auto dark:border-zinc-800">
            <table className="w-full text-sm">
              <thead className="bg-zinc-50 text-zinc-500 text-xs dark:bg-zinc-800/60 dark:text-zinc-400">
                <tr>
                  <th className="text-left font-medium px-3 py-2">Strain</th>
                  <th className="text-left font-medium px-3 py-2">{geneId}</th>
                  <th className="text-left font-medium px-3 py-2">Element present</th>
                  <th className="text-left font-medium px-3 py-2">Reading</th>
                  <th className="text-right font-medium px-3 py-2">Pieces</th>
                  <th className="text-right font-medium px-3 py-2 whitespace-nowrap">Largest piece</th>
                  <th className="text-right font-medium px-3 py-2">Identity</th>
                  <th className="text-right font-medium px-3 py-2">Genome</th>
                </tr>
              </thead>
              <tbody>
                {report.hits.map((h) => {
                  const delta = sourceSize !== undefined ? h.genome_bp - sourceSize : 0;
                  const reading =
                    h.covered_pct >= 90
                      ? "present"
                      : h.covered_pct < 5
                        ? "absent entirely"
                        : `partly present (${fmtBp(h.covered_bp)})`;
                  return (
                    <tr key={h.query_id} className="border-t border-zinc-100 dark:border-zinc-800">
                      <td className="px-3 py-1.5 whitespace-nowrap">
                        {h.query_name}
                        {report.element_kind === "contig" &&
                          h.query_id === report.source_query_id && (
                            <span className="ml-1.5 text-xs text-zinc-400 dark:text-zinc-500">(source)</span>
                          )}
                      </td>
                      <td className="px-3 py-1.5">{h.call ? <CallBadge call={h.call} /> : "-"}</td>
                      <td className="px-3 py-1.5">
                        <AlignedBar pct={h.covered_pct} />
                      </td>
                      <td className="px-3 py-1.5 whitespace-nowrap">{reading}</td>
                      <td className="px-3 py-1.5 text-right tabular-nums">{h.pieces}</td>
                      <td className="px-3 py-1.5 text-right tabular-nums">
                        {h.largest_piece ? fmtBp(h.largest_piece) : "-"}
                      </td>
                      <td className="px-3 py-1.5 text-right tabular-nums whitespace-nowrap">
                        {h.covered_bp ? `${h.identity.toFixed(1)} %` : "-"}
                      </td>
                      <td
                        className="px-3 py-1.5 text-right tabular-nums whitespace-nowrap"
                        title={`${h.genome_bp.toLocaleString("en-US")} bp`}
                      >
                        {fmtBp(h.genome_bp)}
                        {sourceSize !== undefined && delta !== 0 && (
                          <span className="ml-1 text-xs text-zinc-500 dark:text-zinc-400">
                            ({delta > 0 ? "+" : "-"}
                            {fmtBp(Math.abs(delta))})
                          </span>
                        )}
                      </td>
                    </tr>
                  );
                })}
              </tbody>
            </table>
          </div>
          <p className="text-sm text-zinc-500 dark:text-zinc-400">
            A plasmid is often split over several contigs; the contig option compares only the piece
            that carries the gene. For the whole plasmid, compare a complete record from NCBI.
            Missing from an assembly is not proof of missing from the bacterium: a small plasmid can
            be lost in culture, DNA extraction or assembly. Confirm by PCR or by mapping the raw reads.
          </p>
        </>
      )}
    </section>
  );
}
