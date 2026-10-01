import { useEffect, useRef, useState } from "react";
import { api, type ParamSpecLike } from "../api";
import type { ProjectFile, RunParams } from "../types";
import { Button, DropZone, ErrorBox, InfoIcon, Spinner } from "./ui";

export interface WizardResult {
  ran: boolean;
}

const DEFAULT_PARAMS: RunParams = {
  min_gap: 200,
  present_cov: 95,
  partial_cov: 0,
  blast_cov: 90,
  blast_pid: 90,
  blast_evalue: 1e-10,
  nucmer_minmatch: null,
  nucmer_breaklen: null,
  dnadiff: true,
};

/**
 * The guided input dialog: 1 reference, 2 queries, 3 optional panel,
 * 4 review + compare. Big buttons, plain language, immediate validation.
 */
export function InputWizard({
  open,
  onClose,
  projectId,
  files,
  onFilesChanged,
  onCompare,
}: {
  open: boolean;
  onClose: () => void;
  projectId: number;
  files: ProjectFile[];
  onFilesChanged: () => void;
  onCompare: (queryIds: number[], params: RunParams) => void;
}) {
  const [step, setStep] = useState(1);
  const [error, setError] = useState<string | null>(null);
  const [notice, setNotice] = useState<PanelNotice | null>(null);
  const [geneList, setGeneList] = useState("");
  const [busy, setBusy] = useState<string | null>(null);
  const [accession, setAccession] = useState("");
  const [schema, setSchema] = useState<ParamSpecLike[]>([]);
  const [showAdvanced, setShowAdvanced] = useState(false);
  const [params, setParams] = useState<RunParams>(DEFAULT_PARAMS);
  const [selectedQueries, setSelectedQueries] = useState<number[]>([]);
  const ranRef = useRef(false);

  useEffect(() => {
    if (open) {
      setStep(1);
      setError(null);
      setNotice(null);
      ranRef.current = false;
      api
        .getPresets()
        .then((p) => setSchema(p.schema))
        .catch(() => setSchema([]));
    }
  }, [open]);

  useEffect(() => {
    const queries = files.filter((f) => f.role === "query");
    setSelectedQueries((prev) => {
      const ids = queries.map((q) => q.id);
      if (prev.every((p) => ids.includes(p)) && prev.length) return prev;
      return ids;
    });
  }, [files]);

  if (!open) return null;

  const refFasta = files.find((f) => f.role === "reference_fasta");
  const refGff = files.find((f) => f.role === "reference_gff");
  const queries = files.filter((f) => f.role === "query");
  const panel = files.find((f) => f.role === "panel");

  async function uploadReference(fasta: File, gff: File) {
    setBusy("Uploading the reference genome...");
    setError(null);
    try {
      await api.uploadReference(projectId, fasta, gff);
      onFilesChanged();
    } catch (e) {
      setError((e as Error).message);
    } finally {
      setBusy(null);
    }
  }

  async function fetchNcbi() {
    if (!accession.trim()) {
      setError("Please type an NCBI assembly accession, for example GCF_000196035.1.");
      return;
    }
    setBusy("Downloading the reference from NCBI (this can take a minute)...");
    setError(null);
    try {
      await api.fetchReferenceFromNcbi(projectId, accession.trim());
      onFilesChanged();
    } catch (e) {
      setError((e as Error).message);
    } finally {
      setBusy(null);
    }
  }

  async function uploadQueries(fs: File[]) {
    setBusy("Adding the query genomes...");
    setError(null);
    try {
      await api.uploadQueries(projectId, fs);
      onFilesChanged();
    } catch (e) {
      setError((e as Error).message);
    } finally {
      setBusy(null);
    }
  }

  async function uploadPanel(f: File) {
    setBusy("Adding the gene panel...");
    setError(null);
    setNotice(null);
    try {
      let r: PanelResult;
      if (/\.(csv|tsv|txt)$/i.test(f.name)) {
        r = await api.uploadPanelIds(projectId, f);
      } else {
        await api.uploadPanel(projectId, f);
        r = { found: [], from_ncbi: [], from_catalog: [], missing: [] };
      }
      onFilesChanged();
      showPanelNotice(r);
    } catch (e) {
      setError((e as Error).message);
    } finally {
      setBusy(null);
    }
  }

  async function buildPanelFromText(append: boolean) {
    if (!geneList.trim()) {
      setError("Please type the gene names first.");
      return;
    }
    setBusy(
      `${append ? "Adding to" : "Building"} the gene panel (genes not in the reference are taken from the curated AMRFinderPlus and VFDB databases, or else fetched from NCBI)...`,
    );
    setError(null);
    setNotice(null);
    try {
      const r = await api.buildPanelFromText(projectId, geneList, append);
      onFilesChanged();
      showPanelNotice(r, append);
    } catch (e) {
      setError((e as Error).message);
    } finally {
      setBusy(null);
    }
  }

  function showPanelNotice(r: PanelResult, append = false) {
    const catalog = r.from_catalog ?? [];
    const n = r.found.length + r.from_ncbi.length + catalog.length;
    if (n === 0) return;
    const parts = [`${r.found.length} from the reference`];
    if (catalog.length > 0) parts.push(`${catalog.length} from the curated databases`);
    if (r.from_ncbi.length > 0) parts.push(`${r.from_ncbi.length} from NCBI`);
    setNotice({
      summary: append
        ? `Added ${n} gene${n === 1 ? "" : "s"} to the panel: ${parts.join(", ")}.`
        : `The panel was built with ${n} genes: ${parts.join(", ")}.`,
      catalog,
      hints: r.hints ?? [],
      fetched: r.from_ncbi,
      missing: r.missing,
    });
  }

  const steps = [
    { n: 1, label: "Reference" },
    { n: 2, label: "Query genomes" },
    { n: 3, label: "Gene panel (optional)" },
    { n: 4, label: "Compare" },
  ];
  const stepReady = [
    Boolean(refFasta && refGff),
    selectedQueries.length > 0,
    true,
    true,
  ];

  function compare() {
    if (!refFasta || !refGff || selectedQueries.length === 0) return;
    ranRef.current = true;
    onCompare(selectedQueries, params);
    onClose();
  }

  // The step buttons sit both above and below the step: with many query
  // FASTA files the footer is a long scroll away.
  const nav = (big: boolean) => (
    <>
      <Button variant="ghost" onClick={onClose}>
        Cancel
      </Button>
      <div className="flex gap-3">
        {step > 1 && (
          <Button variant="secondary" onClick={() => setStep(step - 1)}>
            Back
          </Button>
        )}
        {step < 4 ? (
          <Button disabled={!stepReady[step - 1]} onClick={() => setStep(step + 1)}>
            Continue
          </Button>
        ) : (
          <Button
            disabled={!stepReady[3] || selectedQueries.length === 0}
            onClick={compare}
            size={big ? "lg" : undefined}
          >
            Compare genomes
          </Button>
        )}
      </div>
    </>
  );

  return (
    <div className="fixed inset-0 z-50 flex items-start justify-center bg-zinc-900/40 p-4 sm:p-8 overflow-y-auto">
      <div className="bg-white rounded-xl shadow-xl border border-zinc-200 w-full max-w-3xl my-auto dark:bg-zinc-900 dark:border-zinc-800">
        {/* header */}
        <div className="px-8 pt-6 pb-4">
          <div className="flex flex-wrap items-center justify-between gap-3">
            <h2 className="text-xl font-semibold">Set up the comparison</h2>
            <div className="flex items-center gap-3">{nav(false)}</div>
          </div>
          <div className="flex gap-2 mt-4" aria-hidden>
            {steps.map((s, i) => (
              <div key={s.n} className="flex items-center flex-1 last:flex-none">
                <div className="flex items-center gap-2">
                  <span
                    className={`w-8 h-8 rounded-full grid place-items-center text-sm font-semibold ${
                      step > s.n
                        ? "bg-zinc-900 text-white dark:bg-zinc-100 dark:text-zinc-900"
                        : step === s.n
                          ? "bg-zinc-900 text-white ring-4 ring-zinc-200 dark:bg-zinc-100 dark:text-zinc-900 dark:ring-zinc-700"
                          : "bg-zinc-100 text-zinc-400 dark:bg-zinc-800 dark:text-zinc-600"
                    }`}
                  >
                    {step > s.n ? "\u2713" : s.n}
                  </span>
                  <span
                    className={`text-sm font-medium hidden sm:block ${step >= s.n ? "text-zinc-900 dark:text-zinc-100" : "text-zinc-400 dark:text-zinc-600"}`}
                  >
                    {s.label}
                  </span>
                </div>
                {i < steps.length - 1 && <div className="flex-1 h-px bg-zinc-200 mx-3 dark:bg-zinc-800" />}
              </div>
            ))}
          </div>
        </div>

        {error && (
          <div className="px-8 pb-2">
            <ErrorBox message={error} />
          </div>
        )}

        <div className="px-8 pb-2 min-h-[260px]">
          {busy && (
            <div className="flex items-center gap-3 text-zinc-600 py-8 dark:text-zinc-400">
              <Spinner /> <span>{busy}</span>
            </div>
          )}

          {!busy && step === 1 && (
            <div className="space-y-4">
              {refFasta && refGff ? (
                <div className="rounded-lg border border-emerald-200 bg-emerald-50 px-4 py-3 flex items-start gap-3 dark:border-emerald-900 dark:bg-emerald-950/40">
                  <span className="text-emerald-700 mt-0.5 dark:text-emerald-400">{"\u2713"}</span>
                  <div>
                    <p className="font-medium text-emerald-900 dark:text-emerald-300">Reference is ready</p>
                    <p className="text-sm text-emerald-800 mt-0.5 dark:text-emerald-400">
                      {refFasta.display_name} and {refGff.display_name}
                    </p>
                  </div>
                </div>
              ) : null}
              <p className="text-[15px] text-zinc-600 dark:text-zinc-400">
                The reference is the annotated genome the others are compared
                against. Provide the genome file (FASTA) and its annotation
                (GFF), or fetch both from NCBI by accession.
              </p>
              <div className="grid sm:grid-cols-2 gap-3">
                <DropZone
                  compact
                  hint="Reference genome (FASTA), e.g. LM259.fasta"
                  onFiles={(fs) => {
                    const gff = fs.find((f) => /\.(gff|gff3|gtf)$/i.test(f.name));
                    const fasta = fs.find(
                      (f) => /\.(fasta|fa|fna|fsa)$/i.test(f.name),
                    );
                    if (fasta && gff) uploadReference(fasta, gff);
                    else
                      setError(
                        "Please provide both files at once: the genome (FASTA) and its annotation (GFF). Drop them together in the box.",
                      );
                  }}
                  accept=".fasta,.fa,.fna,.fsa,.gff,.gff3"
                />
                <div className="rounded-xl border border-zinc-200 p-4 flex flex-col dark:border-zinc-800">
                  <p className="text-sm font-medium text-zinc-700 dark:text-zinc-300">
                    or fetch from NCBI
                  </p>
                  <input
                    value={accession}
                    onChange={(e) => setAccession(e.target.value)}
                    placeholder="GCF_000196035.1"
                    className="mt-2 h-11 px-3 rounded-lg border border-zinc-300 text-[15px] focus:border-zinc-500 outline-none w-full dark:border-zinc-700 dark:bg-zinc-900"
                  />
                  <Button
                    className="mt-3 w-full"
                    variant="secondary"
                    onClick={fetchNcbi}
                  >
                    Search and download
                  </Button>
                  <p className="text-xs text-zinc-400 mt-2 dark:text-zinc-500">
                    Downloads the genome and its annotation together.
                  </p>
                </div>
              </div>
            </div>
          )}

          {!busy && step === 2 && (
            <div className="space-y-4">
              <p className="text-[15px] text-zinc-600 dark:text-zinc-400">
                Add the genomes you want to compare against the reference
                (draft or complete). You can add several at once; they are
                compared in parallel.
              </p>
              <DropZone
                multiple
                hint="Query genomes (FASTA), one or many, e.g. strain_A.fasta"
                onFiles={uploadQueries}
              />
              {queries.length > 0 && (
                <ul className="divide-y divide-zinc-100 border border-zinc-200 rounded-lg dark:divide-zinc-800 dark:border-zinc-800">
                  {queries.map((q) => (
                    <li
                      key={q.id}
                      className="flex items-center justify-between px-4 py-2.5"
                    >
                      <span className="text-[15px]">{q.display_name}</span>
                      <button
                        className="text-sm text-zinc-400 hover:text-red-600 px-2 h-9 rounded dark:text-zinc-500 dark:hover:text-red-400"
                        onClick={async () => {
                          await api.deleteFile(projectId, q.id);
                          onFilesChanged();
                        }}
                      >
                        Remove
                      </button>
                    </li>
                  ))}
                </ul>
              )}
            </div>
          )}

          {!busy && step === 3 && (
            <div className="space-y-4">
              <div className="flex items-center justify-between gap-3">
                <p className="text-[15px] text-zinc-600 dark:text-zinc-400">
                  A gene panel is a set of genes checked extra strictly, in
                  addition to the genome comparison. Optional.
                </p>
                <PanelHelp />
              </div>
              <div className="flex gap-2">
                <textarea
                  className="flex-1 min-h-24 border border-zinc-300 rounded-lg px-3 py-2 text-[15px] font-mono text-sm focus:outline-none focus:border-zinc-500 dark:border-zinc-700 dark:bg-zinc-900"
                  placeholder="inlA, inlB, qacH, lmo0444, ..."
                  value={geneList}
                  onChange={(e) => setGeneList(e.target.value)}
                />
                {panel ? (
                  <div className="self-start flex flex-col gap-2">
                    <button
                      className="h-11 px-4 rounded-lg bg-zinc-900 text-white text-[15px] hover:bg-zinc-700 whitespace-nowrap dark:bg-zinc-100 dark:text-zinc-900 dark:hover:bg-zinc-300"
                      title="Keep the current panel and add these genes to it"
                      onClick={() => buildPanelFromText(true)}
                    >
                      Add to panel
                    </button>
                    <button
                      className="h-11 px-4 rounded-lg border border-zinc-300 text-[15px] text-zinc-700 hover:bg-zinc-100 whitespace-nowrap dark:border-zinc-700 dark:text-zinc-300 dark:hover:bg-zinc-800"
                      title="Discard the current panel and build a new one from these genes"
                      onClick={() => buildPanelFromText(false)}
                    >
                      Replace panel
                    </button>
                  </div>
                ) : (
                  <button
                    className="self-start h-11 px-4 rounded-lg bg-zinc-900 text-white text-[15px] hover:bg-zinc-700 whitespace-nowrap dark:bg-zinc-100 dark:text-zinc-900 dark:hover:bg-zinc-300"
                    onClick={() => buildPanelFromText(false)}
                  >
                    Build panel
                  </button>
                )}
              </div>
              {notice && (
                <div className="rounded-lg border border-zinc-200 bg-zinc-50 px-4 py-3 text-[15px] text-zinc-700 dark:border-zinc-800 dark:bg-zinc-800/60 dark:text-zinc-300">
                  <p>{notice.summary}</p>
                  {notice.hints.length > 0 && (
                    <div className="mt-2 rounded-md border border-amber-200 bg-amber-50 px-3 py-2 text-sm text-amber-900 dark:border-amber-900 dark:bg-amber-950/40 dark:text-amber-200">
                      {notice.hints.map((h) => (
                        <p key={h}>{h}</p>
                      ))}
                    </div>
                  )}
                  {notice.catalog.length > 0 && (
                    <NoticeList
                      title="From the curated databases (AMRFinderPlus, VFDB), as entry: product. Check that each is the gene you meant:"
                      items={notice.catalog}
                    />
                  )}
                  {notice.fetched.length > 0 && (
                    <NoticeList title="Fetched from NCBI:" items={notice.fetched} />
                  )}
                  {notice.missing.length > 0 && (
                    <p className="mt-2">
                      Not found anywhere: {notice.missing.join(", ")}. Pin a GenBank record, e.g.{" "}
                      <span className="font-mono">gene (ACCESSION:start-end rev)</span>, or add a
                      FASTA file.
                    </p>
                  )}
                </div>
              )}
              {panel ? (
                <div className="rounded-lg border border-emerald-200 bg-emerald-50 px-4 py-3 flex items-center justify-between dark:border-emerald-900 dark:bg-emerald-950/40">
                  <div>
                    <p className="font-medium text-emerald-900 dark:text-emerald-300">
                      Panel is ready
                    </p>
                    <p className="text-sm text-emerald-800 mt-0.5 dark:text-emerald-400">
                      {panel.display_name}
                    </p>
                  </div>
                  <button
                    className="text-sm text-zinc-400 hover:text-red-600 px-2 h-9 rounded dark:text-zinc-500 dark:hover:text-red-400"
                    onClick={async () => {
                      await api.deleteFile(projectId, panel.id);
                      onFilesChanged();
                      setNotice(null);
                    }}
                  >
                    Remove
                  </button>
                </div>
              ) : (
                <DropZone
                  compact
                  hint="Gene list (CSV/TSV) or gene panel (FASTA), optional"
                  accept=".fasta,.fa,.fna,.fsa,.csv,.tsv,.txt"
                  onFiles={(fs) => uploadPanel(fs[0])}
                />
              )}
            </div>
          )}

          {!busy && step === 4 && (
            <div className="space-y-4">
              <ul className="text-[15px] divide-y divide-zinc-100 border border-zinc-200 rounded-lg dark:divide-zinc-800 dark:border-zinc-800">
                <ReviewRow
                  ok={Boolean(refFasta && refGff)}
                  label="Reference genome"
                  value={
                    refFasta && refGff
                      ? `${refFasta.display_name} + ${refGff.display_name}`
                      : "missing"
                  }
                />
                <ReviewRow
                  ok={selectedQueries.length > 0}
                  label="Query genomes"
                  value={
                    selectedQueries.length === 1
                      ? queries.find((q) => q.id === selectedQueries[0])
                          ?.display_name ?? "1 genome"
                      : `${selectedQueries.length} genomes`
                  }
                />
                <ReviewRow
                  ok
                  label="Gene panel"
                  value={panel ? panel.display_name : "not used"}
                />
              </ul>

              <div className="border border-zinc-200 rounded-lg dark:border-zinc-800">
                <button
                  className="w-full flex items-center justify-between px-4 h-12 text-[15px] font-medium text-zinc-600 hover:text-zinc-900 dark:text-zinc-400 dark:hover:text-zinc-100"
                  onClick={() => setShowAdvanced(!showAdvanced)}
                >
                  <span className="flex items-center gap-2">
                    Advanced settings
                    <InfoIcon text="Optional thresholds for the comparison. Defaults are sensible for most bacteria; you never need to open this." />
                  </span>
                  <span>{showAdvanced ? "\u25B2" : "\u25BC"}</span>
                </button>
                {showAdvanced && (
                  <div className="px-4 pb-4 pt-1 border-t border-zinc-100 dark:border-zinc-800">
                    <ParamsForm
                      schema={schema}
                      params={params}
                      onChange={setParams}
                    />
                  </div>
                )}
              </div>
            </div>
          )}
        </div>

        {/* footer */}
        <div className="flex items-center justify-between px-8 py-4 border-t border-zinc-200 dark:border-zinc-800">
          {nav(true)}
        </div>
      </div>
    </div>
  );
}

interface PanelResult {
  found: string[];
  from_ncbi: string[];
  from_catalog?: string[];
  hints?: string[];
  missing: string[];
}

interface PanelNotice {
  summary: string;
  hints: string[];
  catalog: string[];
  fetched: string[];
  missing: string[];
}

function NoticeList({ title, items }: { title: string; items: string[] }) {
  return (
    <div className="mt-2">
      <p className="text-sm text-zinc-500 dark:text-zinc-400">{title}</p>
      <ul className="mt-1 space-y-0.5 text-sm">
        {items.map((i) => (
          <li key={i}>{i}</li>
        ))}
      </ul>
    </div>
  );
}

/** Click-to-open tutorial for the gene panel step, kept out of the way until asked for. */
function PanelHelp() {
  const [open, setOpen] = useState(false);
  return (
    <div className="relative shrink-0">
      <button
        type="button"
        onClick={() => setOpen((o) => !o)}
        className="text-sm font-medium text-zinc-500 hover:text-zinc-900 underline underline-offset-2 whitespace-nowrap dark:text-zinc-400 dark:hover:text-zinc-100"
      >
        How does this work?
      </button>
      {open && (
        <div className="absolute right-0 top-8 z-30 w-96 max-w-[80vw] bg-white border border-zinc-200 rounded-lg shadow-lg p-4 space-y-3 text-[15px] text-zinc-600 dark:bg-zinc-900 dark:border-zinc-800 dark:text-zinc-400">
          <p>
            Paste a list of genes (separated by commas or new lines: symbols
            like inlA, locus tags like lmo0444) or drop a CSV file, and the
            sequences are collected automatically, in this order: from your
            reference genome; then from the curated AMRFinderPlus
            (resistance, disinfectant, metal, stress) and VFDB (virulence)
            databases, picking the entry for your project{"\u2019"}s organism;
            and only then from NCBI by name. A ready-made FASTA panel also
            works.
          </p>
          <p>
            Gene names are not unique, so check the list of what was taken
            from where after building. A database entry can be asked for by
            its own name (e.g. cadA_Lm). To take a gene from one exact
            GenBank record instead, pin it:{" "}
            {"\"qacH (HF565366.1)\""} uses the record's own annotation, and
            {" "}{"\"emrC (CP038643.1:1496-1882 rev)\""} pinpoints the exact
            spot when the record does not name the gene.
          </p>
          <button
            type="button"
            onClick={() => setOpen(false)}
            className="text-sm font-medium text-zinc-500 hover:text-zinc-900 dark:text-zinc-400 dark:hover:text-zinc-100"
          >
            Got it
          </button>
        </div>
      )}
    </div>
  );
}

function ReviewRow({
  ok,
  label,
  value,
}: {
  ok: boolean;
  label: string;
  value: string;
}) {
  return (
    <li className="flex items-center justify-between px-4 py-3">
      <span className="text-zinc-500 dark:text-zinc-400">{label}</span>
      <span
        className={`font-medium ${ok ? "text-zinc-900 dark:text-zinc-100" : "text-red-600 dark:text-red-400"}`}
      >
        {value}
      </span>
    </li>
  );
}

export function ParamsForm({
  schema,
  params,
  onChange,
}: {
  schema: ParamSpecLike[];
  params: RunParams;
  onChange: (p: RunParams) => void;
}) {
  return (
    <div className="space-y-4 pt-3">
      {schema.map((spec) => (
        <div key={spec.name} className="grid sm:grid-cols-[1fr_170px] gap-2 sm:gap-4 items-center">
          <div>
            <label className="text-sm font-medium text-zinc-800 dark:text-zinc-200">
              {spec.label}
              {spec.layer === "align" && (
                <span className="ml-2 text-xs text-zinc-400 dark:text-zinc-500">
                  (slower to change)
                </span>
              )}
            </label>
            <p className="text-xs text-zinc-500 mt-0.5 dark:text-zinc-400">{spec.help}</p>
          </div>
          {spec.kind.kind === "bool" ? (
            <label className="flex items-center gap-2 h-11 cursor-pointer">
              <input
                type="checkbox"
                className="w-5 h-5 accent-zinc-900 dark:accent-zinc-100"
                checked={params[spec.name as keyof RunParams] as boolean}
                onChange={(e) =>
                  onChange({ ...params, [spec.name]: e.target.checked })
                }
              />
              <span className="text-sm text-zinc-600 dark:text-zinc-400">
                {params[spec.name as keyof RunParams] ? "On" : "Off"}
              </span>
            </label>
          ) : spec.kind.kind === "optional_int" ? (
            <input
              type="number"
              className="w-full h-11 px-3 rounded-lg border border-zinc-300 text-[15px] dark:border-zinc-700 dark:bg-zinc-900"
              placeholder="tool default"
              value={
                spec.kind.kind === "optional_int" && params[spec.name as keyof RunParams] === null
                  ? ""
                  : String(params[spec.name as keyof RunParams] ?? "")
              }
              onChange={(e) =>
                onChange({
                  ...params,
                  [spec.name]: e.target.value === "" ? null : Number(e.target.value),
                })
              }
            />
          ) : (
            <input
              type="number"
              className="w-full h-11 px-3 rounded-lg border border-zinc-300 text-[15px] dark:border-zinc-700 dark:bg-zinc-900"
              value={params[spec.name as keyof RunParams] as number}
              step={spec.kind.kind === "float" ? "any" : "1"}
              onChange={(e) =>
                onChange({ ...params, [spec.name]: Number(e.target.value) })
              }
            />
          )}
        </div>
      ))}
    </div>
  );
}
