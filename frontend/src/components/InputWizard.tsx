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
  partial_cov: 1,
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
  const [notice, setNotice] = useState<string | null>(null);
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
      if (/\.(csv|tsv|txt)$/i.test(f.name)) {
        const r = await api.uploadPanelIds(projectId, f);
        onFilesChanged();
        if (r.missing.length > 0) {
          setNotice(
            `The panel was built with ${r.found.length} of ${
              r.found.length + r.missing.length
            } genes. Not found in the reference: ${r.missing.join(", ")}.`,
          );
        } else {
          setNotice(
            `The panel was built from all ${r.found.length} genes of the list.`,
          );
        }
      } else {
        await api.uploadPanel(projectId, f);
        onFilesChanged();
      }
    } catch (e) {
      setError((e as Error).message);
    } finally {
      setBusy(null);
    }
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

  return (
    <div className="fixed inset-0 z-50 flex items-start justify-center bg-zinc-900/40 p-4 sm:p-8 overflow-y-auto">
      <div className="bg-white rounded-xl shadow-xl border border-zinc-200 w-full max-w-3xl my-auto">
        {/* header */}
        <div className="px-8 pt-6 pb-4">
          <h2 className="text-xl font-semibold">Set up the comparison</h2>
          <div className="flex gap-2 mt-4" aria-hidden>
            {steps.map((s, i) => (
              <div key={s.n} className="flex items-center flex-1 last:flex-none">
                <div className="flex items-center gap-2">
                  <span
                    className={`w-8 h-8 rounded-full grid place-items-center text-sm font-semibold ${
                      step > s.n
                        ? "bg-zinc-900 text-white"
                        : step === s.n
                          ? "bg-zinc-900 text-white ring-4 ring-zinc-200"
                          : "bg-zinc-100 text-zinc-400"
                    }`}
                  >
                    {step > s.n ? "\u2713" : s.n}
                  </span>
                  <span
                    className={`text-sm font-medium hidden sm:block ${step >= s.n ? "text-zinc-900" : "text-zinc-400"}`}
                  >
                    {s.label}
                  </span>
                </div>
                {i < steps.length - 1 && <div className="flex-1 h-px bg-zinc-200 mx-3" />}
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
            <div className="flex items-center gap-3 text-zinc-600 py-8">
              <Spinner /> <span>{busy}</span>
            </div>
          )}

          {!busy && step === 1 && (
            <div className="space-y-4">
              {refFasta && refGff ? (
                <div className="rounded-lg border border-emerald-200 bg-emerald-50 px-4 py-3 flex items-start gap-3">
                  <span className="text-emerald-700 mt-0.5">{"\u2713"}</span>
                  <div>
                    <p className="font-medium text-emerald-900">Reference is ready</p>
                    <p className="text-sm text-emerald-800 mt-0.5">
                      {refFasta.display_name} and {refGff.display_name}
                    </p>
                  </div>
                </div>
              ) : null}
              <p className="text-[15px] text-zinc-600">
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
                <div className="rounded-xl border border-zinc-200 p-4 flex flex-col">
                  <p className="text-sm font-medium text-zinc-700">
                    or fetch from NCBI
                  </p>
                  <input
                    value={accession}
                    onChange={(e) => setAccession(e.target.value)}
                    placeholder="GCF_000196035.1"
                    className="mt-2 h-11 px-3 rounded-lg border border-zinc-300 text-[15px] focus:border-zinc-500 outline-none w-full"
                  />
                  <Button
                    className="mt-3 w-full"
                    variant="secondary"
                    onClick={fetchNcbi}
                  >
                    Search and download
                  </Button>
                  <p className="text-xs text-zinc-400 mt-2">
                    Downloads the genome and its annotation together.
                  </p>
                </div>
              </div>
            </div>
          )}

          {!busy && step === 2 && (
            <div className="space-y-4">
              <p className="text-[15px] text-zinc-600">
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
                <ul className="divide-y divide-zinc-100 border border-zinc-200 rounded-lg">
                  {queries.map((q) => (
                    <li
                      key={q.id}
                      className="flex items-center justify-between px-4 py-2.5"
                    >
                      <span className="text-[15px]">{q.display_name}</span>
                      <button
                        className="text-sm text-zinc-400 hover:text-red-600 px-2 h-9 rounded"
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
              <p className="text-[15px] text-zinc-600">
                A gene panel is a set of genes checked extra strictly (a
                precise search, in addition to the genome comparison). This is
                optional.
              </p>
              <p className="text-[15px] text-zinc-600">
                Drop a list of genes (CSV/TSV, one per line: locus tags like
                lmo0444, symbols like inlA, or {"\"pva (lmo0446)\""} ) and the
                gene sequences are collected from the reference automatically.
                A ready-made FASTA panel works too.
              </p>
              {notice && (
                <div className="rounded-lg border border-zinc-200 bg-zinc-50 px-4 py-3 text-[15px] text-zinc-700">
                  {notice}
                </div>
              )}
              {panel ? (
                <div className="rounded-lg border border-emerald-200 bg-emerald-50 px-4 py-3 flex items-center justify-between">
                  <div>
                    <p className="font-medium text-emerald-900">
                      Panel is ready
                    </p>
                    <p className="text-sm text-emerald-800 mt-0.5">
                      {panel.display_name}
                    </p>
                  </div>
                  <button
                    className="text-sm text-zinc-400 hover:text-red-600 px-2 h-9 rounded"
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
              <ul className="text-[15px] divide-y divide-zinc-100 border border-zinc-200 rounded-lg">
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

              <div className="border border-zinc-200 rounded-lg">
                <button
                  className="w-full flex items-center justify-between px-4 h-12 text-[15px] font-medium text-zinc-600 hover:text-zinc-900"
                  onClick={() => setShowAdvanced(!showAdvanced)}
                >
                  <span className="flex items-center gap-2">
                    Advanced settings
                    <InfoIcon text="Optional thresholds for the comparison. Defaults are sensible for most bacteria; you never need to open this." />
                  </span>
                  <span>{showAdvanced ? "\u25B2" : "\u25BC"}</span>
                </button>
                {showAdvanced && (
                  <div className="px-4 pb-4 pt-1 border-t border-zinc-100">
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
        <div className="flex items-center justify-between px-8 py-4 border-t border-zinc-200">
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
              <Button
                disabled={!stepReady[step - 1]}
                onClick={() => setStep(step + 1)}
              >
                Continue
              </Button>
            ) : (
              <Button
                disabled={!stepReady[3] || selectedQueries.length === 0}
                onClick={compare}
                size="lg"
              >
                Compare genomes
              </Button>
            )}
          </div>
        </div>
      </div>
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
      <span className="text-zinc-500">{label}</span>
      <span
        className={`font-medium ${ok ? "text-zinc-900" : "text-red-600"}`}
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
            <label className="text-sm font-medium text-zinc-800">
              {spec.label}
              {spec.layer === "align" && (
                <span className="ml-2 text-xs text-zinc-400">
                  (slower to change)
                </span>
              )}
            </label>
            <p className="text-xs text-zinc-500 mt-0.5">{spec.help}</p>
          </div>
          {spec.kind.kind === "bool" ? (
            <label className="flex items-center gap-2 h-11 cursor-pointer">
              <input
                type="checkbox"
                className="w-5 h-5 accent-zinc-900"
                checked={params[spec.name as keyof RunParams] as boolean}
                onChange={(e) =>
                  onChange({ ...params, [spec.name]: e.target.checked })
                }
              />
              <span className="text-sm text-zinc-600">
                {params[spec.name as keyof RunParams] ? "On" : "Off"}
              </span>
            </label>
          ) : spec.kind.kind === "optional_int" ? (
            <input
              type="number"
              className="w-full h-11 px-3 rounded-lg border border-zinc-300 text-[15px]"
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
              className="w-full h-11 px-3 rounded-lg border border-zinc-300 text-[15px]"
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
