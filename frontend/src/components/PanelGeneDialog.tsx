import { Fragment, useEffect, useMemo, useState } from "react";
import { api } from "../api";
import type {
  ContextGene,
  ElementGene,
  ElementHit,
  ElementReport,
  GeneOrigin,
  OriginRecord,
  PanelContext,
  Run,
} from "../types";
import { CallBadge, ErrorBox, Modal, RelatedChip, Spinner } from "./ui";
import { PanelAlignmentDialog } from "./PanelAlignmentDialog";

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
  // a plasmid picked in the NCBI section, handed to the comparison
  const [preset, setPreset] = useState<{ acc: string; n: number } | null>(null);
  useEffect(() => setPreset(null), [geneId]);
  // the strain whose alignment the Investigate dialog shows
  const [investigate, setInvestigate] = useState<number | null>(null);
  useEffect(() => setInvestigate(null), [geneId]);

  return (
    <Modal
      open={Boolean(geneId)}
      onClose={onClose}
      extraWide
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
            onInvestigate={setInvestigate}
          />
          <GeneOriginSection
            run={run}
            geneId={geneId}
            onCompare={(acc) => setPreset({ acc, n: Date.now() })}
          />
          <AcrossStrains
            run={run}
            geneId={geneId}
            preferredSource={queryId}
            preset={preset}
          />
        </div>
      )}
      {geneId && (
        <PanelAlignmentDialog
          run={run}
          geneId={geneId}
          queryId={investigate}
          onClose={() => setInvestigate(null)}
        />
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
  onInvestigate,
}: {
  run: Run;
  geneId: string;
  queryId: number | undefined;
  onQuery: (id: number) => void;
  onInvestigate: (queryId: number) => void;
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

          {ctx.match_note && (
            <div className="rounded-lg border border-amber-200 bg-amber-50 px-4 py-3 text-[15px] text-amber-900 dark:border-amber-900 dark:bg-amber-950/40 dark:text-amber-200">
              <p className="font-semibold">
                {ctx.call === "PARTIAL" ? "Partial match" : "Closest match"}
              </p>
              <p className="mt-1">{ctx.match_note}</p>
            </div>
          )}

          {ctx.contig && queryId !== undefined && (
            <button
              onClick={() => onInvestigate(queryId)}
              className="h-10 px-4 rounded-lg border border-zinc-300 bg-white text-sm font-medium hover:bg-zinc-100 dark:border-zinc-700 dark:bg-zinc-900 dark:hover:bg-zinc-800"
              title="Align the panel gene base by base to this strain, next to another strain"
            >
              Investigate alignment
            </button>
          )}

          {ctx.contig && (
            <dl className="grid grid-cols-2 sm:grid-cols-4 gap-x-6 gap-y-3 text-sm">
              <Fact label="Call">
                <CallBadge call={ctx.call} />
                {ctx.call !== "ABSENT" && (
                  <span className="ml-2 text-zinc-500 dark:text-zinc-400">
                    {ctx.cov_pct.toFixed(0)} % cov., {ctx.identity.toFixed(1)} % id.
                  </span>
                )}
              </Fact>
              {ctx.related_identity != null && (
                <Fact label="Related gene">
                  <RelatedChip identity={ctx.related_identity} name={ctx.related_gene} />
                  {ctx.related_gene && (
                    <span className="ml-2 text-zinc-700 dark:text-zinc-300">{ctx.related_gene}</span>
                  )}
                </Fact>
              )}
              <Fact label="Contig">
                <span className="font-mono">{ctx.contig.seqid}</span>
                <span className="ml-1 text-zinc-500 dark:text-zinc-400">({fmtBp(ctx.contig.length)})</span>
              </Fact>
              <Fact label={ctx.call === "PRESENT" ? "Hit" : "Closest match"}>
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

          {ctx.contig && (
            <div>
              <p className="text-sm font-medium text-zinc-700 dark:text-zinc-300">
                Genes within {fmtBp(ctx.window)} of the {ctx.call === "PRESENT" ? "hit" : "closest match"}
              </p>
              <p className="text-xs text-zinc-500 dark:text-zinc-400 mb-2">
                Where this strain matches the reference, the reference{"\u2019"}s own genes are
                shown at their aligned place. In stretches the reference lacks, genes are predicted
                and named by similarity.
              </p>
              {ctx.genes.length > 0 ? (
                <ContextGenes ctx={ctx} />
              ) : (
                <p className="text-sm text-zinc-500 dark:text-zinc-400">
                  No genes were found around the hit.
                </p>
              )}
              {ctx.genes_note && (
                <p className="mt-2 text-xs text-zinc-500 dark:text-zinc-400">{ctx.genes_note}</p>
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

const SOURCE_LABEL: Record<string, string> = {
  annotation: "reference annotation",
  panel: "panel gene",
  reference: "predicted, like a reference gene",
  ncbi: "predicted, named by NCBI",
  "": "predicted, unnamed",
};

/** How closely a predicted gene matches the protein it is named after. */
function MatchStats({ g }: { g: ContextGene }) {
  return (
    <span
      className="text-xs text-zinc-500 tabular-nums dark:text-zinc-400"
      title="Amino-acid identity over the aligned part, and the share of this gene the alignment covers."
    >
      {g.match_identity!.toFixed(0)} % identity
      {g.match_coverage != null && <> over {g.match_coverage.toFixed(0)} % of the gene</>}
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
            <th className="text-left font-medium px-3 py-2">Source</th>
            <th className="text-right font-medium px-3 py-2">Distance</th>
          </tr>
        </thead>
        <tbody>
          {ctx.genes.map((g) => (
            <tr
              key={`${g.source}:${g.locus_tag}:${g.start}-${g.end}`}
              className={`border-t border-zinc-100 dark:border-zinc-800 ${
                g.is_hit ? "bg-emerald-50/60 dark:bg-emerald-950/20" : ""
              }`}
            >
              <td className="px-3 py-1.5 font-mono tabular-nums whitespace-nowrap">
                {g.start.toLocaleString("en-US")}-{g.end.toLocaleString("en-US")} ({g.strand < 0 ? "-" : "+"})
              </td>
              <td className="px-3 py-1.5 leading-6">
                {g.label ||
                  (!g.locus_tag && <span className="text-zinc-400 dark:text-zinc-500">unnamed</span>)}
                {g.locus_tag && g.locus_tag !== g.label && (
                  <> <span className="font-mono text-xs text-zinc-500 dark:text-zinc-400">{g.locus_tag}</span></>
                )}
                {g.is_hit && (
                  <> <span className="inline-block px-1.5 rounded text-xs bg-emerald-100 text-emerald-800 dark:bg-emerald-900/40 dark:text-emerald-300">
                      {g.source === "panel"
                        ? "hit"
                        : ctx.call === "PRESENT"
                          ? ctx.gene_id
                          : `${ctx.call === "PARTIAL" ? "partial match" : "closest match"} to ${ctx.gene_id}`}
                    </span></>
                )}
                {g.source !== "panel" && g.label && g.match_identity != null && (
                  <> <span className="text-xs text-zinc-500 dark:text-zinc-400">(<MatchStats g={g} />)</span></>
                )}
                {g.mobile && (
                  <> <span className="inline-block px-1.5 rounded text-xs bg-violet-100 text-violet-800 dark:bg-violet-900/40 dark:text-violet-300">
                      mobile element
                    </span></>
                )}
                {g.partial && (
                  <> <span
                      className="text-xs text-zinc-400 dark:text-zinc-500"
                      title="Only part of this gene is here: it runs past the edge of the aligned stretch or of the region the gene was predicted in."
                    >
                      (partial)
                    </span></>
                )}
                {g.source === "panel" && g.match_label && g.match_identity != null && (
                  <div className="text-xs text-zinc-500 dark:text-zinc-400">
                    most similar named protein: {g.match_label}, <MatchStats g={g} />
                  </div>
                )}
              </td>
              <td
                className="px-3 py-1.5 text-xs text-zinc-500 whitespace-nowrap dark:text-zinc-400"
                title={
                  g.source === "panel"
                    ? "Identified by matching the panel's curated sequence; the gene's start and end come from gene prediction (Prodigal)."
                    : undefined
                }
              >
                {g.source === "panel" && ctx.panel_source
                  ? ctx.panel_source
                  : (SOURCE_LABEL[g.source] ?? g.source)}
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
  preset,
}: {
  run: Run;
  geneId: string;
  preferredSource: number | undefined;
  preset: { acc: string; n: number } | null;
}) {
  // the strains carrying the gene, from the panel matrix row
  const [positives, setPositives] = useState<Set<number> | null>(null);
  const [source, setSource] = useState<number | undefined>(undefined);
  const [useAccession, setUseAccession] = useState(false);
  const [accession, setAccession] = useState("");
  const [report, setReport] = useState<ElementReport | null>(null);
  // the strain whose share of the element is shown in detail
  const [openRow, setOpenRow] = useState<number | null>(null);
  useEffect(() => setOpenRow(null), [report]);
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
        // only strains carrying the gene in full have its element; a
        // Partial match sits in a related gene, in other DNA
        row?.calls.forEach((c, i) => {
          if (c === "PRESENT") pos.add(run.queries[i].file_id);
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

  async function compare(presetAcc?: string) {
    const acc = presetAcc ?? (useAccession ? accession.trim() : undefined);
    // a record needs no strain carrying the gene; the strain then only
    // anchors the genome-size differences
    const src = source ?? (acc ? run.queries[0]?.file_id : undefined);
    if (src === undefined) return;
    setBusy(true);
    setError(null);
    setReport(null);
    try {
      setReport(await api.panelElement(run.id, geneId, src, acc));
    } catch (e) {
      setError((e as Error).message);
    } finally {
      setBusy(false);
    }
  }

  useEffect(() => {
    if (!preset) return;
    setUseAccession(true);
    setAccession(preset.acc);
    compare(preset.acc);
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [preset]);

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
          {geneId} is not present in full in any strain of this run, so there is no element to
          compare. You can still compare a complete record from NCBI:
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
              DNA carrying the gene in
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
          onClick={() => compare()}
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
            {report.element_kind === "region" && <>the insertion carrying {geneId}, </>}
            <span className="font-mono">{report.element_name}</span>
            {report.element_title && <> - {report.element_title}</>} ({fmtBp(report.element_len)})
          </p>
          {report.element_kind === "region" && (
            <p className="text-xs text-zinc-500 dark:text-zinc-400">
              {geneId} sits in a stretch the reference lacks, inside a contig that is otherwise
              shared, so that stretch is compared rather than the whole contig.
            </p>
          )}
          <p className="text-xs text-zinc-500 dark:text-zinc-400">
            Click a strain to see which parts of the element it holds
            {report.genes && report.genes.length > 0 ? ", gene by gene" : ""}.
          </p>
          <div className="border border-zinc-200 rounded-lg overflow-hidden dark:border-zinc-800">
            <table className="w-full text-sm">
              <thead className="bg-zinc-50 text-zinc-500 text-xs dark:bg-zinc-800/60 dark:text-zinc-400">
                <tr>
                  <th className="text-left font-medium px-3 py-2">Strain</th>
                  <th className="text-left font-medium px-3 py-2">{geneId}</th>
                  <th className="text-left font-medium px-3 py-2">Element present</th>
                  <th className="text-right font-medium px-3 py-2">Pieces</th>
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
                  const isOpen = openRow === h.query_id;
                  return (
                    <Fragment key={h.query_id}>
                      <tr
                        onClick={() => setOpenRow(isOpen ? null : h.query_id)}
                        aria-expanded={isOpen}
                        className={`border-t border-zinc-100 cursor-pointer hover:bg-zinc-50 dark:border-zinc-800 dark:hover:bg-zinc-800/50 ${isOpen ? "bg-zinc-50 dark:bg-zinc-800/50" : ""}`}
                      >
                        <td className="px-3 py-1.5 break-all">
                          <span className="inline-flex items-start gap-1.5">
                            <span
                              aria-hidden
                              className={`mt-0.5 text-[10px] text-zinc-400 transition-transform ${isOpen ? "rotate-90" : ""}`}
                            >
                              ▶
                            </span>
                            <span>
                              {h.query_name}
                              {report.element_kind !== "accession" &&
                                h.query_id === report.source_query_id && (
                                  <span className="ml-1.5 text-xs text-zinc-400 dark:text-zinc-500">
                                    (source)
                                  </span>
                                )}
                            </span>
                          </span>
                        </td>
                        <td className="px-3 py-1.5">
                          <span className="inline-flex flex-wrap items-center gap-1">
                            {h.call ? <CallBadge call={h.call} /> : "-"}
                            {h.related_identity != null && (
                              <RelatedChip identity={h.related_identity} />
                            )}
                          </span>
                        </td>
                        <td className="px-3 py-1.5">
                          <AlignedBar pct={h.covered_pct} />
                          <span className="block text-xs text-zinc-500 dark:text-zinc-400">
                            {reading}
                          </span>
                        </td>
                        <td className="px-3 py-1.5 text-right tabular-nums">
                          {h.pieces}
                          {h.largest_piece > 0 && (
                            <span className="block text-xs text-zinc-500 whitespace-nowrap dark:text-zinc-400">
                              largest {fmtBp(h.largest_piece)}
                            </span>
                          )}
                        </td>
                        <td className="px-3 py-1.5 text-right tabular-nums">
                          {h.covered_bp ? `${h.identity.toFixed(1)} %` : "-"}
                        </td>
                        <td
                          className="px-3 py-1.5 text-right tabular-nums"
                          title={`${h.genome_bp.toLocaleString("en-US")} bp`}
                        >
                          {fmtBp(h.genome_bp)}
                          {sourceSize !== undefined && delta !== 0 && (
                            <span className="block text-xs text-zinc-500 dark:text-zinc-400">
                              {delta > 0 ? "+" : "-"}
                              {fmtBp(Math.abs(delta))}
                            </span>
                          )}
                        </td>
                      </tr>
                      {isOpen && (
                        <tr className="bg-zinc-50/60 dark:bg-zinc-800/30">
                          <td colSpan={6} className="px-3 pt-2 pb-4">
                            <ElementDetail report={report} hit={h} />
                          </td>
                        </tr>
                      )}
                    </Fragment>
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

type GeneShare = "present" | "partly" | "missing";

/** How much of one element gene a strain's stretches cover, and how alike. */
function geneCover(g: ElementGene, hit: ElementHit): { pct: number; identity: number | null } {
  const len = g.end - g.start + 1;
  let covered = 0;
  let weighted = 0;
  for (const p of hit.covered ?? []) {
    const o = Math.min(g.end, p.end) - Math.max(g.start, p.start) + 1;
    if (o > 0) {
      covered += o;
      weighted += o * p.identity;
    }
  }
  return { pct: (100 * covered) / len, identity: covered > 0 ? weighted / covered : null };
}

function shareOf(pct: number): GeneShare {
  return pct >= 90 ? "present" : pct >= 10 ? "partly" : "missing";
}

const SHARE_FILL: Record<GeneShare, string> = {
  present: "fill-emerald-500 dark:fill-emerald-400",
  partly: "fill-amber-400 dark:fill-amber-500",
  missing: "fill-zinc-300 dark:fill-zinc-600",
};

const SHARE_TEXT: Record<GeneShare, string> = {
  present: "text-emerald-700 dark:text-emerald-400",
  partly: "text-amber-700 dark:text-amber-400",
  missing: "text-zinc-500 dark:text-zinc-400",
};

/** One strain's share of the element: where its stretches lie and, for an
 * annotated record, which of the record's genes they hold. */
function ElementDetail({ report, hit }: { report: ElementReport; hit: ElementHit }) {
  const [filter, setFilter] = useState<"all" | GeneShare>("all");
  const genes = report.genes ?? [];
  const pieces = hit.covered ?? [];
  const rows = useMemo(
    () =>
      genes.map((g) => {
        const c = geneCover(g, hit);
        return { g, ...c, share: shareOf(c.pct) };
      }),
    [genes, hit],
  );
  const count = (k: GeneShare) => rows.filter((r) => r.share === k).length;
  const shown = filter === "all" ? rows : rows.filter((r) => r.share === filter);
  const len = Math.max(1, report.element_len);
  const x = (bp: number) => (1000 * (bp - 1)) / len;

  return (
    <div className="space-y-3">
      <svg
        viewBox="0 0 1000 40"
        preserveAspectRatio="none"
        className="w-full h-10"
        role="img"
        aria-label={`Map of the element: genes above, the stretches ${hit.query_name} holds below`}
      >
        <rect x="0" y="27" width="1000" height="8" rx="2" className="fill-zinc-200 dark:fill-zinc-700" />
        {pieces.map((p) => (
          <rect
            key={`${p.start}-${p.end}`}
            x={x(p.start)}
            y="27"
            width={Math.max(1, x(p.end + 1) - x(p.start))}
            height="8"
            className="fill-sky-600 dark:fill-sky-400"
            opacity={p.identity >= 99 ? 1 : p.identity >= 95 ? 0.7 : 0.45}
          >
            <title>
              {`${p.start.toLocaleString("en-US")}-${p.end.toLocaleString("en-US")} (${fmtBp(p.end - p.start + 1)}), ${p.identity.toFixed(1)} % identity`}
            </title>
          </rect>
        ))}
        {rows.map(({ g, pct, share }) => (
          <rect
            key={`${g.start}-${g.end}-${g.locus_tag}`}
            x={x(g.start)}
            y={g.strand < 0 ? 13 : 2}
            width={Math.max(1.5, x(g.end + 1) - x(g.start))}
            height="9"
            rx="1"
            className={`${SHARE_FILL[share]} ${g.name.toLowerCase() === report.gene_id.toLowerCase() ? "stroke-zinc-900 dark:stroke-zinc-100" : ""}`}
            strokeWidth="1.5"
          >
            <title>{`${g.name || g.locus_tag}: ${g.product || "no product"} (${pct.toFixed(0)} % held)`}</title>
          </rect>
        ))}
      </svg>
      <p className="flex flex-wrap gap-x-4 gap-y-1 text-xs text-zinc-500 dark:text-zinc-400">
        {genes.length > 0 && (
          <span>Genes above (forward strand on top), coloured by how much of each {hit.query_name} holds.</span>
        )}
        <span>
          <span className="inline-block w-3 h-2 align-middle rounded-sm bg-sky-600 dark:bg-sky-400" /> stretches
          found (paler = less alike)
        </span>
      </p>

      {genes.length > 0 ? (
        <>
          <div className="flex flex-wrap items-center gap-2 text-xs">
            {(
              [
                ["all", `All ${rows.length}`],
                ["present", `Present ${count("present")}`],
                ["partly", `Partly ${count("partly")}`],
                ["missing", `Missing ${count("missing")}`],
              ] as const
            ).map(([k, label]) => (
              <button
                key={k}
                onClick={() => setFilter(k)}
                className={`h-7 px-2.5 rounded-md border ${filter === k ? "border-zinc-900 bg-zinc-900 text-white dark:border-zinc-100 dark:bg-zinc-100 dark:text-zinc-900" : "border-zinc-300 text-zinc-700 hover:bg-zinc-100 dark:border-zinc-700 dark:text-zinc-300 dark:hover:bg-zinc-800"}`}
              >
                {label}
              </button>
            ))}
          </div>
          <div className="max-h-[28rem] overflow-y-auto border border-zinc-200 rounded-lg bg-white dark:border-zinc-800 dark:bg-zinc-900">
            <table className="w-full text-sm">
              <thead className="sticky top-0 bg-zinc-50 text-zinc-500 text-xs dark:bg-zinc-800 dark:text-zinc-400">
                <tr>
                  <th className="text-left font-medium px-3 py-2">Gene</th>
                  <th className="text-left font-medium px-3 py-2">Product</th>
                  <th className="text-right font-medium px-3 py-2">Place</th>
                  <th className="text-left font-medium px-3 py-2">In {hit.query_name}</th>
                </tr>
              </thead>
              <tbody>
                {shown.map(({ g, pct, identity, share }) => (
                  <tr
                    key={`${g.start}-${g.end}-${g.locus_tag}`}
                    className="border-t border-zinc-100 align-top dark:border-zinc-800"
                  >
                    <td className="px-3 py-1.5">
                      <span className={g.name ? "italic" : "font-mono text-xs"}>
                        {g.name || g.locus_tag || "-"}
                      </span>
                      {g.name && g.locus_tag && (
                        <span className="block font-mono text-xs text-zinc-500 dark:text-zinc-400">
                          {g.locus_tag}
                        </span>
                      )}
                    </td>
                    <td className="px-3 py-1.5 text-zinc-700 dark:text-zinc-300">{g.product || "-"}</td>
                    <td className="px-3 py-1.5 text-right tabular-nums text-xs text-zinc-500 dark:text-zinc-400">
                      {g.start.toLocaleString("en-US")}-{g.end.toLocaleString("en-US")}
                      <span className="block">{g.strand < 0 ? "reverse" : "forward"}</span>
                    </td>
                    <td className={`px-3 py-1.5 tabular-nums ${SHARE_TEXT[share]}`}>
                      {share === "present" ? "present" : share === "partly" ? "partly" : "missing"}
                      <span className="block text-xs text-zinc-500 dark:text-zinc-400">
                        {pct.toFixed(0)} % of it
                        {identity != null && <>, {identity.toFixed(1)} % identity</>}
                      </span>
                    </td>
                  </tr>
                ))}
              </tbody>
            </table>
          </div>
        </>
      ) : (
        <>
          <p className="text-sm text-zinc-600 dark:text-zinc-400">
            {report.element_kind === "accession"
              ? "This record has no annotated genes to list, so only the stretches are shown."
              : "A strain's own contig carries no gene annotation, so only the stretches are shown. To see them gene by gene, compare a complete record, e.g. a library plasmid from the section above."}
          </p>
          {pieces.length > 0 && (
            <ul className="text-xs tabular-nums text-zinc-600 space-y-0.5 dark:text-zinc-400">
              {pieces.map((p) => (
                <li key={`${p.start}-${p.end}`}>
                  {p.start.toLocaleString("en-US")}-{p.end.toLocaleString("en-US")} ({fmtBp(p.end - p.start + 1)}),{" "}
                  {p.identity.toFixed(1)} % identity
                </li>
              ))}
            </ul>
          )}
        </>
      )}
    </div>
  );
}

/** NCBI asks for at most one poll a minute per search; the server
 * enforces that, so polling a little faster only catches the answer
 * sooner after it lands. */
const ORIGIN_POLL_MS = 20_000;

function GeneOriginSection({
  run,
  geneId,
  onCompare,
}: {
  run: Run;
  geneId: string;
  onCompare: (accession: string) => void;
}) {
  const [origin, setOrigin] = useState<GeneOrigin | null>(null);
  const [started, setStarted] = useState(false);
  const [wide, setWide] = useState(false);
  const [attempt, setAttempt] = useState(0);
  const [error, setError] = useState<string | null>(null);

  useEffect(() => {
    setOrigin(null);
    setStarted(false);
    setWide(false);
    setError(null);
  }, [run.id, geneId]);

  useEffect(() => {
    if (!started) return;
    let cancelled = false;
    let timer: ReturnType<typeof setTimeout> | undefined;
    setError(null);
    const tick = () => {
      api
        .panelOrigin(run.id, geneId, wide)
        .then((o) => {
          if (cancelled) return;
          setOrigin(o);
          if (o.state === "running") timer = setTimeout(tick, ORIGIN_POLL_MS);
        })
        .catch((e) => {
          if (!cancelled) setError((e as Error).message);
        });
    };
    tick();
    return () => {
      cancelled = true;
      if (timer) clearTimeout(timer);
    };
  }, [started, wide, attempt, run.id, geneId]);

  const done = origin?.state === "done";
  const scope = origin?.scope === "all bacteria" ? "bacterial" : origin?.scope;

  return (
    <section className="space-y-3">
      <SectionTitle>Where does this gene usually occur?</SectionTitle>
      <p className="text-sm text-zinc-600 dark:text-zinc-400">
        Searches NCBI for other genomes that carry {geneId} over its whole length and counts how
        many of them have it on a plasmid or on the chromosome. Only the gene{"’"}s sequence is
        sent to NCBI, never your strains. This describes other genomes, not your strains.
      </p>
      {!started && (
        <button
          onClick={() => setStarted(true)}
          className="h-9 px-4 rounded-lg bg-zinc-900 text-white text-sm font-medium hover:bg-zinc-700 dark:bg-zinc-100 dark:text-zinc-900 dark:hover:bg-zinc-300"
        >
          Look it up in NCBI
        </button>
      )}
      {error && (
        <div className="space-y-2">
          <ErrorBox message={error} />
          <button
            onClick={() => setAttempt((n) => n + 1)}
            className="h-9 px-4 rounded-lg border border-zinc-300 bg-white text-sm font-medium hover:bg-zinc-100 dark:border-zinc-700 dark:bg-zinc-900 dark:hover:bg-zinc-800"
          >
            Try again
          </button>
        </div>
      )}
      {started && !done && !error && (
        <p className="flex items-center gap-2 text-sm text-zinc-500 dark:text-zinc-400">
          <Spinner className="text-zinc-400" />
          {origin?.message || "Asking NCBI BLAST..."}{" "}
          {wide
            ? "A search of all bacteria often waits 10 to 30 minutes in NCBI\u2019s queue."
            : "This usually takes one to three minutes."}{" "}
          You can close the dialog; reopen it and look the gene up again to pick up the same search.
        </p>
      )}
      {origin?.searched_with && (
        <p className="text-sm text-zinc-600 dark:text-zinc-400">
          Searched with the {geneId} sequence the strains matched best:{" "}
          <span className="font-medium text-zinc-800 dark:text-zinc-200">{origin.searched_with}</span>.
        </p>
      )}
      {done && origin && (
        <>
          {origin.n_matches === 0 ? (
            <p className="text-[15px] text-zinc-700 dark:text-zinc-300">{origin.message}</p>
          ) : (
            <div className="rounded-lg border border-zinc-200 bg-zinc-50 px-4 py-3 text-[15px] dark:border-zinc-800 dark:bg-zinc-800/60">
              <p className="text-zinc-900 dark:text-zinc-100">
                Found in <b>{origin.n_matches}</b> {scope} record{origin.n_matches === 1 ? "" : "s"}:{" "}
                <b>{origin.n_plasmid}</b> on a plasmid, <b>{origin.n_chromosome}</b> on a chromosome
                {origin.n_contig > 0 && (
                  <>
                    , <b>{origin.n_contig}</b> on draft-assembly contigs that NCBI does not place
                  </>
                )}
                .
              </p>
              <p className="mt-1 text-sm text-zinc-600 dark:text-zinc-400">
                {originReading(origin)} {origin.scope_note}
              </p>
            </div>
          )}
          {origin.n_matches === 0 && (
            <p className="text-sm text-zinc-600 dark:text-zinc-400">{origin.scope_note}</p>
          )}
          {origin.can_widen && (
            <button
              onClick={() => {
                setOrigin(null);
                setWide(true);
              }}
              className="h-9 px-4 rounded-lg border border-zinc-300 bg-white text-sm font-medium hover:bg-zinc-100 dark:border-zinc-700 dark:bg-zinc-900 dark:hover:bg-zinc-800"
            >
              Search all bacteria (slow: often 10 to 30 minutes)
            </button>
          )}
          {origin.plasmids.length > 0 && (
            <div>
              <p className="text-sm font-medium text-zinc-700 dark:text-zinc-300 mb-2">
                Complete plasmids carrying {geneId}
              </p>
              <OriginRecords records={origin.plasmids} onCompare={onCompare} />
            </div>
          )}
          {origin.chromosomes.length > 0 && (
            <div>
              <p className="text-sm font-medium text-zinc-700 dark:text-zinc-300 mb-2">
                Chromosomes carrying {geneId}
              </p>
              <OriginRecords records={origin.chromosomes} />
            </div>
          )}
        </>
      )}
    </section>
  );
}

/** A plain-language reading of the counts, never stronger than they are. */
function originReading(o: GeneOrigin): string {
  const placed = o.n_plasmid + o.n_chromosome;
  if (placed === 0) return "None of the records says whether it is a plasmid or a chromosome.";
  const share = o.n_plasmid / placed;
  if (share >= 0.8) return "Mostly plasmid-borne.";
  if (share <= 0.2) return "Mostly chromosomal.";
  return "Found on both plasmids and chromosomes.";
}

function OriginRecords({
  records,
  onCompare,
}: {
  records: OriginRecord[];
  onCompare?: (accession: string) => void;
}) {
  return (
    <div className="border border-zinc-200 rounded-lg overflow-x-auto dark:border-zinc-800">
      <table className="w-full text-sm">
        <thead className="bg-zinc-50 text-zinc-500 text-xs dark:bg-zinc-800/60 dark:text-zinc-400">
          <tr>
            <th className="text-left font-medium px-3 py-2">Accession</th>
            <th className="text-left font-medium px-3 py-2">Description</th>
            <th className="text-right font-medium px-3 py-2">Length</th>
            <th
              className="text-right font-medium px-3 py-2 whitespace-nowrap"
              title="DNA identity of the best match to the searched sequence; listed records match at least 90 % of it at 90 % identity or more"
            >
              DNA identity
            </th>
            {onCompare && <th />}
          </tr>
        </thead>
        <tbody>
          {records.map((r) => (
            <tr key={r.accession} className="border-t first:border-t-0 border-zinc-100 dark:border-zinc-800">
              <td className="px-3 py-1.5 whitespace-nowrap">
                <a
                  href={`https://www.ncbi.nlm.nih.gov/nuccore/${r.accession}`}
                  target="_blank"
                  rel="noreferrer"
                  className="font-mono underline hover:text-zinc-900 dark:hover:text-zinc-100"
                >
                  {r.accession}
                </a>
              </td>
              <td className="px-3 py-1.5 text-zinc-700 dark:text-zinc-300">{r.title}</td>
              <td className="px-3 py-1.5 text-right tabular-nums whitespace-nowrap">{fmtBp(r.length)}</td>
              <td className="px-3 py-1.5 text-right tabular-nums whitespace-nowrap">
                {r.identity.toFixed(1)} %
              </td>
              {onCompare && (
                <td className="px-3 py-1.5 text-right">
                  <button
                    onClick={() => onCompare(r.accession)}
                    className="h-8 px-3 rounded-lg border border-zinc-300 bg-white text-xs font-medium whitespace-nowrap hover:bg-zinc-100 dark:border-zinc-700 dark:bg-zinc-900 dark:hover:bg-zinc-800"
                    title="Compare this plasmid with every strain of the run (section below)"
                  >
                    Compare with my strains
                  </button>
                </td>
              )}
            </tr>
          ))}
        </tbody>
      </table>
    </div>
  );
}
