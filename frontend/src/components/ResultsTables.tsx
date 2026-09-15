import { useEffect, useMemo, useRef, useState } from "react";
import { useVirtualizer } from "@tanstack/react-virtual";
import { api } from "../api";
import type {
  Call,
  GapRow,
  GeneDetail,
  GeneCoverageRow,
  MatrixRow,
  Page,
  PanelRow,
  Run,
  TableQuery,
} from "../types";
import { CallBadge, Spinner } from "./ui";

export type TableKind = "genes_coverage" | "unaligned_gaps" | "panel_recheck" | "matrix";

interface Column {
  key: string;
  label: string;
  numeric?: boolean;
  width?: number;
}

const COVERAGE_COLUMNS: Column[] = [
  { key: "locus_tag", label: "Locus tag", width: 150 },
  { key: "symbol", label: "Symbol", width: 120 },
  { key: "biotype", label: "Biotype", width: 130 },
  { key: "protein_id", label: "Protein", width: 140 },
  { key: "seqid", label: "Sequence", width: 120 },
  { key: "start", label: "Start", numeric: true, width: 90 },
  { key: "end", label: "End", numeric: true, width: 90 },
  { key: "length", label: "Length", numeric: true, width: 90 },
  { key: "cov_bp", label: "Covered bases", numeric: true, width: 110 },
  { key: "cov_pct", label: "Coverage %", numeric: true, width: 100 },
  { key: "call", label: "Call", width: 120 },
  { key: "best_identity", label: "Best identity %", numeric: true, width: 110 },
  { key: "mismatches", label: "Mismatches", numeric: true, width: 100 },
  { key: "indels", label: "Indels", numeric: true, width: 80 },
];

const GAP_COLUMNS: Column[] = [
  { key: "seqid", label: "Sequence", width: 140 },
  { key: "start", label: "Start", numeric: true, width: 100 },
  { key: "end", label: "End", numeric: true, width: 100 },
  { key: "length", label: "Length", numeric: true, width: 100 },
  { key: "n_genes", label: "Genes inside", numeric: true, width: 110 },
  { key: "genes", label: "Locus tags of genes inside", width: 420 },
];

const PANEL_COLUMNS: Column[] = [
  { key: "gene_id", label: "Gene", width: 180 },
  { key: "qlen", label: "Length", numeric: true, width: 100 },
  { key: "cov_pct", label: "Coverage %", numeric: true, width: 110 },
  { key: "identity", label: "Identity %", numeric: true, width: 110 },
  { key: "best_evalue", label: "Best match significance", numeric: true, width: 140 },
  { key: "call", label: "Call", width: 120 },
];

const PAGE_SIZE = 200;
const COL_WIDTHS_KEY = "bactiment-col-widths";
const HIDDEN_COLS_KEY = "bactiment-hidden-cols";

// Columns hidden by default per table, until the user changes it via the
// Columns picker (then their choice is remembered instead).
const DEFAULT_HIDDEN_COLS: Partial<Record<TableKind, string[]>> = {
  genes_coverage: ["start", "end", "length", "cov_bp"],
};

export function ResultsTables({
  run,
  initialTable,
  initialQueryId,
  onStateChange,
  onOpenGene,
}: {
  run: Run;
  initialTable: TableKind;
  initialQueryId: number | undefined;
  onStateChange: (table: TableKind, queryId: number | undefined) => void;
  onOpenGene: (locus: string) => void;
}) {
  const safeInitialTable: TableKind =
    initialTable === "panel_recheck" && !run.has_panel
      ? "genes_coverage"
      : initialTable;
  const [table, setTable] = useState<TableKind>(safeInitialTable);
  const [queryId, setQueryId] = useState<number | undefined>(initialQueryId);
  const [search, setSearch] = useState("");
  const [debouncedSearch, setDebouncedSearch] = useState("");
  const [call, setCall] = useState<string>("");
  const [sortBy, setSortBy] = useState<string | null>(null);
  const [sortDir, setSortDir] = useState<"asc" | "desc">("asc");
  const [page, setPage] = useState(0);
  const [hiddenCols, setHiddenCols] = useState<Record<string, string[]>>(() => {
    try {
      return JSON.parse(localStorage.getItem(HIDDEN_COLS_KEY) ?? "{}");
    } catch {
      return {};
    }
  });
  const [data, setData] = useState<Page<GapRow | GeneCoverageRow | PanelRow | MatrixRow> | null>(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [pinnedGene, setPinnedGene] = useState<string | null>(null);
  const [colWidths, setColWidths] = useState<Record<string, number>>(() => {
    try {
      return JSON.parse(localStorage.getItem(COL_WIDTHS_KEY) ?? "{}");
    } catch {
      return {};
    }
  });

  useEffect(() => {
    try {
      localStorage.setItem(COL_WIDTHS_KEY, JSON.stringify(colWidths));
    } catch {}
  }, [colWidths]);

  useEffect(() => {
    try {
      localStorage.setItem(HIDDEN_COLS_KEY, JSON.stringify(hiddenCols));
    } catch {}
  }, [hiddenCols]);

  function resizeColumn(key: string, width: number) {
    setColWidths((prev) => ({ ...prev, [`${table}:${key}`]: width }));
  }

  const hidden = hiddenCols[table] ?? DEFAULT_HIDDEN_COLS[table] ?? [];
  function setHiddenForTable(next: string[]) {
    setHiddenCols((prev) => ({ ...prev, [table]: next }));
  }

  useEffect(() => {
    setTable(safeInitialTable);
    setQueryId(initialQueryId);
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [initialTable, initialQueryId, run.id]);

  useEffect(() => {
    const t = setTimeout(() => setDebouncedSearch(search), 350);
    return () => clearTimeout(t);
  }, [search]);

  useEffect(() => {
    onStateChange(table, queryId);
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [table, queryId]);

  const columns: Column[] = useMemo(() => {
    if (table === "genes_coverage") return COVERAGE_COLUMNS;
    if (table === "unaligned_gaps") return GAP_COLUMNS;
    if (table === "panel_recheck") return PANEL_COLUMNS;
    // matrix: fixed columns + one per query
    return [
      { key: "locus_tag", label: "Locus tag", width: 150 },
      { key: "symbol", label: "Symbol", width: 120 },
      { key: "biotype", label: "Biotype", width: 130 },
      ...run.queries.map((q) => ({
        key: `q_${q.file_id}`,
        label: q.name,
        width: 140,
      })),
    ];
  }, [table, run]);

  const visibleColumns = useMemo(
    () => columns.filter((c) => !hidden.includes(c.key)),
    [columns, hidden],
  );

  const query: TableQuery = useMemo(
    () => ({
      query_id: queryId,
      page,
      page_size: PAGE_SIZE,
      sort_by: sortBy ?? undefined,
      sort_dir: sortBy ? sortDir : undefined,
      search: debouncedSearch || undefined,
      call: call || undefined,
    }),
    [queryId, page, sortBy, sortDir, debouncedSearch, call],
  );

  useEffect(() => {
    setPage(0);
    setSortBy(null);
    setCall("");
  }, [table]);

  useEffect(() => {
    let cancelled = false;
    setLoading(true);
    setError(null);
    const fetcher =
      table === "genes_coverage"
        ? api.genesCoverage(run.id, query)
        : table === "unaligned_gaps"
          ? api.unalignedGaps(run.id, query)
          : table === "panel_recheck"
            ? api.panelRecheck(run.id, query)
            : api.matrix(run.id, query);
    fetcher
      .then((d) => {
        if (!cancelled) setData(d);
      })
      .catch((e) => {
        if (!cancelled) setError((e as Error).message);
      })
      .finally(() => {
        if (!cancelled) setLoading(false);
      });
    return () => {
      cancelled = true;
    };
  }, [run.id, table, JSON.stringify(query)]);

  const total = data?.total ?? 0;
  const from = total === 0 ? 0 : page * PAGE_SIZE + 1;
  const to = Math.min(total, (page + 1) * PAGE_SIZE);

  function toggleSort(key: string) {
    if (sortBy === key) {
      setSortDir(sortDir === "asc" ? "desc" : "asc");
    } else {
      setSortBy(key);
      setSortDir("asc");
    }
  }

  return (
    <div>
      {/* toolbar */}
      <div className="flex flex-wrap items-center gap-2 py-3">
        <div className="flex rounded-lg border border-zinc-300 overflow-hidden h-11 dark:border-zinc-700">
          {(
            [
              ["genes_coverage", "Genes coverage"],
              ["unaligned_gaps", "Unaligned gaps"],
              ["panel_recheck", "Panel recheck"],
              ["matrix", "Presence / absence"],
            ] as [TableKind, string][]
          )
            .filter(([k]) => k !== "panel_recheck" || run.has_panel)
            .map(([k, label]) => (
            <button
              key={k}
              onClick={() => setTable(k)}
              className={`px-4 text-[15px] font-medium transition-colors ${
                table === k
                  ? "bg-zinc-900 text-white dark:bg-zinc-100 dark:text-zinc-900"
                  : "bg-white text-zinc-600 hover:bg-zinc-100 dark:bg-zinc-900 dark:text-zinc-400 dark:hover:bg-zinc-800"
              }`}
            >
              {label}
            </button>
          ))}
        </div>

        {(table === "genes_coverage" || table === "panel_recheck") &&
          run.queries.length > 1 && (
            <select
              value={queryId ?? ""}
              onChange={(e) => setQueryId(Number(e.target.value) || undefined)}
              className="h-11 px-3 rounded-lg border border-zinc-300 bg-white text-[15px] dark:border-zinc-700 dark:bg-zinc-900"
            >
              {run.queries.map((q) => (
                <option key={q.file_id} value={q.file_id}>
                  {q.name}
                </option>
              ))}
            </select>
          )}

        <input
          value={search}
          onChange={(e) => setSearch(e.target.value)}
          placeholder="Search genes..."
          className="h-11 px-3 rounded-lg border border-zinc-300 text-[15px] w-56 dark:border-zinc-700 dark:bg-zinc-900"
        />

        {table === "genes_coverage" && (
          <div className="flex rounded-lg border border-zinc-300 overflow-hidden h-11 dark:border-zinc-700">
            {[
              ["", "All"],
              ["present", "Present"],
              ["partial", "Partial"],
              ["absent", "Absent"],
            ].map(([v, label]) => (
              <button
                key={v}
                onClick={() => setCall(v)}
                className={`px-3 text-sm font-medium transition-colors ${
                  call === v
                    ? "bg-zinc-900 text-white dark:bg-zinc-100 dark:text-zinc-900"
                    : "bg-white text-zinc-600 hover:bg-zinc-100 dark:bg-zinc-900 dark:text-zinc-400 dark:hover:bg-zinc-800"
                }`}
              >
                {label}
              </button>
            ))}
          </div>
        )}

        {table === "matrix" && run.queries.length > 1 && (
          <div className="flex rounded-lg border border-zinc-300 overflow-hidden h-11 dark:border-zinc-700">
            {[
              ["", "All genes"],
              ["not_present", "Not present everywhere"],
            ].map(([v, label]) => (
              <button
                key={v}
                onClick={() => setCall(v)}
                className={`px-3 text-sm font-medium transition-colors ${
                  call === v
                    ? "bg-zinc-900 text-white dark:bg-zinc-100 dark:text-zinc-900"
                    : "bg-white text-zinc-600 hover:bg-zinc-100 dark:bg-zinc-900 dark:text-zinc-400 dark:hover:bg-zinc-800"
                }`}
              >
                {label}
              </button>
            ))}
          </div>
        )}

        <ColumnPicker
          columns={columns}
          hidden={hidden}
          onChange={setHiddenForTable}
        />

        <ExportButton
          run={run}
          table={table}
          query={{ ...query, page: undefined, page_size: undefined }}
        />

        {loading && <Spinner className="text-zinc-400" />}
      </div>

      {error && (
        <p className="text-red-700 text-[15px] py-2 dark:text-red-400">{error}</p>
      )}

      {/* table */}
      <div className="border border-zinc-200 rounded-xl bg-white overflow-hidden dark:border-zinc-800 dark:bg-zinc-900">
        <div className="overflow-x-auto">
          <VirtualTable
            columns={visibleColumns}
            rows={data?.rows ?? []}
            table={table}
            sortBy={sortBy}
            sortDir={sortDir}
            onSort={toggleSort}
            onPinGene={setPinnedGene}
            onOpenGene={onOpenGene}
            run={run}
            colWidths={colWidths}
            onResizeColumn={resizeColumn}
          />
        </div>
        {/* pagination */}
        <div className="flex items-center justify-between px-4 h-12 border-t border-zinc-200 text-sm text-zinc-500 dark:border-zinc-800 dark:text-zinc-400">
          <span>
            {total === 0
              ? "No rows match the current filters"
              : `Showing ${from}-${to} of ${total}`}
          </span>
          <div className="flex gap-2">
            <button
              disabled={page === 0}
              onClick={() => setPage(page - 1)}
              className="h-9 px-3 rounded-md border border-zinc-300 disabled:opacity-40 hover:bg-zinc-100 dark:border-zinc-700 dark:hover:bg-zinc-800"
            >
              Previous
            </button>
            <button
              disabled={to >= total}
              onClick={() => setPage(page + 1)}
              className="h-9 px-3 rounded-md border border-zinc-300 disabled:opacity-40 hover:bg-zinc-100 dark:border-zinc-700 dark:hover:bg-zinc-800"
            >
              Next
            </button>
          </div>
        </div>
      </div>

      {/* right-click preview, pinned until another row is right-clicked or this is closed */}
      {pinnedGene && table === "genes_coverage" && (
        <GenePreview
          runId={run.id}
          locus={pinnedGene}
          onOpen={() => onOpenGene(pinnedGene)}
          onClose={() => setPinnedGene(null)}
        />
      )}
    </div>
  );
}

function VirtualTable({
  columns,
  rows,
  table,
  sortBy,
  sortDir,
  onSort,
  onPinGene,
  onOpenGene,
  run,
  colWidths,
  onResizeColumn,
}: {
  columns: Column[];
  rows: (GapRow | GeneCoverageRow | PanelRow | MatrixRow)[];
  table: TableKind;
  sortBy: string | null;
  sortDir: "asc" | "desc";
  onSort: (key: string) => void;
  onPinGene: (locus: string | null) => void;
  onOpenGene: (locus: string) => void;
  run: Run;
  colWidths: Record<string, number>;
  onResizeColumn: (key: string, width: number) => void;
}) {
  const parentRef = useRef<HTMLDivElement>(null);
  const virtualizer = useVirtualizer({
    count: rows.length,
    getScrollElement: () => parentRef.current,
    estimateSize: () => 44,
    overscan: 12,
  });

  // A manually resized column locks to its exact width (flex-shrink/grow 0);
  // an untouched one keeps stretching to fill the table, as before.
  function flexStyleFor(c: Column): React.CSSProperties {
    const baseWidth = c.width ?? 120;
    const override = colWidths[`${table}:${c.key}`];
    const width = override ?? baseWidth;
    return {
      flex: override != null ? `0 0 ${width}px` : `${baseWidth} 1 ${baseWidth}px`,
      minWidth: width,
    };
  }

  function startResize(e: React.MouseEvent<HTMLDivElement>, key: string) {
    e.preventDefault();
    e.stopPropagation();
    // Start from the column's actual on-screen width (it may currently be
    // stretched by flex-grow beyond its configured base width), not the
    // logical base width, so the column doesn't jump under the cursor.
    const startWidth = e.currentTarget.parentElement!.getBoundingClientRect().width;
    const startX = e.clientX;
    function onMove(ev: MouseEvent) {
      onResizeColumn(key, Math.max(50, Math.round(startWidth + (ev.clientX - startX))));
    }
    function onUp() {
      window.removeEventListener("mousemove", onMove);
      window.removeEventListener("mouseup", onUp);
    }
    window.addEventListener("mousemove", onMove);
    window.addEventListener("mouseup", onUp);
  }

  return (
    <div className="min-w-full">
      {/* header */}
      <div className="flex bg-zinc-50 border-b border-zinc-200 sticky top-0 z-10 dark:bg-zinc-900 dark:border-zinc-800">
        {columns.map((c) => (
          <div
            key={c.key}
            className="relative flex items-stretch border-r border-zinc-200 last:border-r-0 dark:border-zinc-800"
            style={flexStyleFor(c)}
          >
            <button
              onClick={() => onSort(c.key)}
              className={`flex-1 min-w-0 flex items-center gap-1 px-3 h-11 text-left text-xs font-semibold uppercase tracking-wide text-zinc-500 hover:text-zinc-900 dark:text-zinc-400 dark:hover:text-zinc-100 ${c.numeric ? "justify-end" : ""}`}
            >
              <span className="truncate">{c.label}</span>
              {sortBy === c.key && (
                <span className="text-zinc-900 dark:text-zinc-100">{sortDir === "asc" ? "\u2191" : "\u2193"}</span>
              )}
            </button>
            <div
              onMouseDown={(e) => startResize(e, c.key)}
              title="Drag to resize"
              className="w-1.5 shrink-0 cursor-col-resize hover:bg-blue-400/60 active:bg-blue-500"
            />
          </div>
        ))}
      </div>
      {/* rows */}
      <div ref={parentRef} className="overflow-auto thin-scroll" style={{ height: "min(62vh, 700px)" }}>
        <div style={{ height: virtualizer.getTotalSize(), position: "relative" }}>
          {virtualizer.getVirtualItems().map((v) => {
            const row = rows[v.index] as unknown as Record<string, unknown>;
            return (
              <div
                key={v.key}
                className={`flex items-center border-b border-zinc-100 text-[15px] dark:border-zinc-800 ${
                  v.index % 2 ? "bg-zinc-50/60 dark:bg-zinc-800/30" : "bg-white dark:bg-zinc-900"
                } hover:bg-blue-50/50 dark:hover:bg-blue-950/30`}
                style={{
                  position: "absolute",
                  top: v.start,
                  left: 0,
                  width: "100%",
                  height: v.size,
                }}
                onContextMenu={(e) => {
                  if (table === "genes_coverage") {
                    e.preventDefault();
                    const locus = row["locus_tag"] as string;
                    onPinGene(locus);
                  }
                }}
                onClick={() => {
                  if (table === "genes_coverage" || table === "matrix") {
                    const locus = row["locus_tag"] as string;
                    if (locus) onOpenGene(locus);
                  }
                }}
              >
                {columns.map((c) => (
                  <div
                    key={c.key}
                    className={`px-3 flex items-center truncate border-r border-zinc-100 last:border-r-0 dark:border-zinc-800 ${
                      c.numeric ? "justify-end font-mono text-sm tabular-nums" : ""
                    }`}
                    style={flexStyleFor(c)}
                  >
                    <Cell
                      col={c.key}
                      row={row}
                      table={table}
                      run={run}
                    />
                  </div>
                ))}
              </div>
            );
          })}
        </div>
      </div>
    </div>
  );
}

function Cell({
  col,
  row,
  table,
  run,
}: {
  col: string;
  row: Record<string, unknown>;
  table: TableKind;
  run: Run;
}) {
  const v = row[col];
  if (col === "call") return <CallBadge call={v as Call} />;
  if (col === "protein_id" && v) {
    return (
      <a
        href={`https://www.ncbi.nlm.nih.gov/protein/${v}`}
        target="_blank"
        rel="noreferrer"
        title={v as string}
        className="truncate underline text-zinc-700 hover:text-zinc-900 dark:text-zinc-300 dark:hover:text-zinc-100"
        onClick={(e) => e.stopPropagation()}
      >
        {v as string}
      </a>
    );
  }
  if (col === "genes") {
    const genes = v as string[];
    if (!genes || genes.length === 0) return <span className="text-zinc-300 dark:text-zinc-700">-</span>;
    return (
      <span className="truncate">
        {genes.slice(0, 6).join(", ")}
        {genes.length > 6 && (
          <span className="text-zinc-400 dark:text-zinc-500"> +{genes.length - 6} more</span>
        )}
      </span>
    );
  }
  if (col.startsWith("q_") && table === "matrix") {
    const idx = run.queries.findIndex((q) => `q_${q.file_id}` === col);
    const calls = row["calls"] as Call[];
    const covs = row["cov_pcts"] as number[];
    if (idx >= 0 && calls) {
      return (
        <span title={`Coverage ${covs[idx]?.toFixed(1) ?? 0}%`}>
          <CallBadge call={calls[idx]} />
        </span>
      );
    }
  }
  if (typeof v === "number") {
    if (col === "cov_pct" || col === "identity" || col === "best_identity" || col === "cov_pcts")
      return <span>{v.toFixed(2)}</span>;
    return <span>{v.toLocaleString("en-US")}</span>;
  }
  if (v === null || v === undefined || v === "") return <span className="text-zinc-300 dark:text-zinc-700">-</span>;
  return <span className="truncate">{String(v)}</span>;
}

/** Closes an open popover on Escape or on any click outside `ref`'s subtree. */
function usePopoverDismiss(open: boolean, onClose: () => void) {
  const ref = useRef<HTMLDivElement>(null);
  useEffect(() => {
    if (!open) return;
    function onPointerDown(e: MouseEvent) {
      if (ref.current && !ref.current.contains(e.target as Node)) onClose();
    }
    function onKeyDown(e: KeyboardEvent) {
      if (e.key === "Escape") onClose();
    }
    document.addEventListener("mousedown", onPointerDown);
    document.addEventListener("keydown", onKeyDown);
    return () => {
      document.removeEventListener("mousedown", onPointerDown);
      document.removeEventListener("keydown", onKeyDown);
    };
  }, [open, onClose]);
  return ref;
}

function ColumnPicker({
  columns,
  hidden,
  onChange,
}: {
  columns: Column[];
  hidden: string[];
  onChange: (hidden: string[]) => void;
}) {
  const [open, setOpen] = useState(false);
  const ref = usePopoverDismiss(open, () => setOpen(false));
  return (
    <div className="relative" ref={ref}>
      <button
        onClick={() => setOpen(!open)}
        className="h-11 px-3 rounded-lg border border-zinc-300 bg-white text-[15px] hover:bg-zinc-100 dark:border-zinc-700 dark:bg-zinc-900 dark:hover:bg-zinc-800"
      >
        Columns
      </button>
      {open && (
        <div className="absolute right-0 top-12 z-30 bg-white border border-zinc-200 rounded-lg shadow-lg p-3 w-64 dark:bg-zinc-900 dark:border-zinc-800">
          <p className="text-xs font-semibold text-zinc-400 uppercase mb-2 dark:text-zinc-500">Visible columns</p>
          <div className="space-y-1 max-h-80 overflow-y-auto">
            {columns.map((c) => {
              const checked = !hidden.includes(c.key);
              return (
                <label key={c.key} className="flex items-center gap-2 h-9 px-2 rounded hover:bg-zinc-50 cursor-pointer dark:hover:bg-zinc-800">
                  <input
                    type="checkbox"
                    className="w-4 h-4 accent-zinc-900"
                    checked={checked}
                    onChange={() => {
                      onChange(
                        checked
                          ? [...hidden, c.key]
                          : hidden.filter((k) => k !== c.key),
                      );
                    }}
                  />
                  <span className="text-sm">{c.label}</span>
                </label>
              );
            })}
          </div>
          {hidden.length > 0 && (
            <button
              className="mt-2 w-full h-9 text-sm text-zinc-500 hover:text-zinc-900 dark:text-zinc-400 dark:hover:text-zinc-100"
              onClick={() => onChange([])}
            >
              Show all columns
            </button>
          )}
        </div>
      )}
    </div>
  );
}

function ExportButton({
  run,
  table,
  query,
}: {
  run: Run;
  table: TableKind;
  query: TableQuery;
}) {
  const [open, setOpen] = useState(false);
  const allCols = table === "genes_coverage" ? undefined : undefined;
  return (
    <div className="relative">
      <button
        onClick={() => setOpen(!open)}
        className="h-11 px-4 rounded-lg bg-zinc-900 text-white text-[15px] font-medium hover:bg-zinc-700 inline-flex items-center gap-2 dark:bg-zinc-100 dark:text-zinc-900 dark:hover:bg-zinc-300"
      >
        <svg width="14" height="14" viewBox="0 0 16 16" fill="none">
          <path d="M8 2v8m0 0l-3-3m3 3l3-3M3 13h10" stroke="currentColor" strokeWidth="1.8" strokeLinecap="round" strokeLinejoin="round" />
        </svg>
        Export
      </button>
      {open && (
        <div className="absolute right-0 top-12 z-30 bg-white border border-zinc-200 rounded-lg shadow-lg p-2 w-72 dark:bg-zinc-900 dark:border-zinc-800">
          {["tsv", "csv"].map((format) => (
            <a
              key={format}
              href={api.exportUrl(run.id, table, query, format)}
              onClick={() => setOpen(false)}
              className="flex items-center justify-between px-3 h-11 rounded-md hover:bg-zinc-100 text-[15px] dark:hover:bg-zinc-800"
            >
              <span>{format.toUpperCase()}</span>
              <span className="text-xs text-zinc-400 dark:text-zinc-500">current view</span>
            </a>
          ))}
          {table === "genes_coverage" &&
            ["tsv", "csv"].map((format) => (
              <a
                key={`all-${format}`}
                href={api.exportUrl(run.id, table, { ...query, cols: allCols }, format)}
                onClick={() => setOpen(false)}
                className="flex items-center justify-between px-3 h-11 rounded-md hover:bg-zinc-100 text-[15px] dark:hover:bg-zinc-800"
              >
                <span>{format.toUpperCase()}</span>
                <span className="text-xs text-zinc-400 dark:text-zinc-500">all columns</span>
              </a>
            ))}
        </div>
      )}
    </div>
  );
}

/** Right-click preview: quick stats for each query, fetched once per gene. */
function GenePreview({
  runId,
  locus,
  onOpen,
  onClose,
}: {
  runId: number;
  locus: string;
  onOpen: () => void;
  onClose: () => void;
}) {
  const cache = GenePreviewCache.get(runId, locus);
  const [detail, setDetail] = useState<GeneDetail | null>(cache);
  useEffect(() => {
    if (cache) return;
    let cancelled = false;
    api
      .geneDetail(runId, locus)
      .then((d) => {
        if (!cancelled) {
          GenePreviewCache.set(runId, locus, d);
          setDetail(d);
        }
      })
      .catch(() => {});
    return () => {
      cancelled = true;
    };
  }, [runId, locus, cache]);
  if (!detail) return null;
  return (
    <div className="fixed right-4 top-56 bottom-6 z-40 w-72 flex flex-col bg-white border border-zinc-200 rounded-xl shadow-xl p-4 pointer-events-auto overflow-y-auto thin-scroll dark:bg-zinc-900 dark:border-zinc-800">
      <div className="flex items-start justify-between gap-2">
        <p className="font-semibold truncate" title={detail.locus_tag}>
          {detail.locus_tag}
        </p>
        <button
          onClick={onClose}
          aria-label="Close"
          className="shrink-0 w-6 h-6 grid place-items-center rounded text-zinc-400 hover:text-zinc-900 hover:bg-zinc-100 dark:text-zinc-500 dark:hover:text-zinc-100 dark:hover:bg-zinc-800"
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
      <p className="text-xs text-zinc-400 mt-0.5 dark:text-zinc-500">
        {detail.symbol ? `${detail.symbol} - ` : ""}
        {detail.biotype}
        <br />
        {detail.seqid}:{detail.start.toLocaleString("en-US")}-
        {detail.end.toLocaleString("en-US")}
        {detail.protein_id && (
          <>
            <br />
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
      </p>
      <button
        onClick={onOpen}
        className="mt-3 h-9 w-full rounded-md bg-zinc-900 text-white text-sm hover:bg-zinc-700 dark:bg-zinc-100 dark:text-zinc-900 dark:hover:bg-zinc-300"
      >
        Show alignment
      </button>
      <div className="mt-3 space-y-3">
        {detail.queries.map((q) => (
          <div key={q.query_id} className="pt-2 border-t border-zinc-100 text-sm dark:border-zinc-800">
            <p className="font-medium truncate" title={q.query_name}>
              {q.query_name}
            </p>
            <div className="mt-1 grid grid-cols-2 gap-x-2 gap-y-0.5 text-xs">
              <span className="text-zinc-400 dark:text-zinc-500">Coverage</span>
              <span className="text-right font-mono">{q.cov_pct.toFixed(1)}%</span>
              <span className="text-zinc-400 dark:text-zinc-500">Identity</span>
              <span className="text-right font-mono">
                {q.best_identity > 0 ? `${q.best_identity.toFixed(1)}%` : "-"}
              </span>
              <span className="text-zinc-400 dark:text-zinc-500">Mism. / indels</span>
              <span className="text-right font-mono">
                {q.mismatches} / {q.indels}
              </span>
            </div>
          </div>
        ))}
      </div>
    </div>
  );
}

const genePreviewCache = new Map<string, GeneDetail>();
const GenePreviewCache = {
  get(runId: number, locus: string) {
    return genePreviewCache.get(`${runId}:${locus}`) ?? null;
  },
  set(runId: number, locus: string, d: GeneDetail) {
    genePreviewCache.set(`${runId}:${locus}`, d);
  },
};
