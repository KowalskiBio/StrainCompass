import { useEffect, useMemo, useRef, useState } from "react";
import { useVirtualizer } from "@tanstack/react-virtual";
import { api } from "../api";
import type {
  Call,
  GainedBlastHit,
  GainedIdentify,
  GainedRow,
  GainedVerify,
  GapRow,
  GeneDetail,
  GeneCoverageRow,
  IdentifiedOrf,
  MatrixRow,
  Page,
  PanelRow,
  Run,
  TableQuery,
} from "../types";
import { nuccoreRangeUrl } from "../types";
import { CallBadge, Spinner, usePopoverDismiss } from "./ui";

export type TableKind =
  | "genes_coverage"
  | "unaligned_gaps"
  | "gained"
  | "panel_recheck"
  | "matrix";

/** Every row shape the table can render. */
type AnyRow = GapRow | GeneCoverageRow | PanelRow | MatrixRow | GainedRow;

/** Stable identity of a gained region within one query's table. */
function gainedKey(row: GainedRow): string {
  return `${row.qry_seqid}:${row.start}`;
}

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

const GAINED_COLUMNS: Column[] = [
  { key: "qry_seqid", label: "Query contig", width: 150 },
  { key: "start", label: "Start", numeric: true, width: 100 },
  { key: "end", label: "End", numeric: true, width: 100 },
  { key: "length", label: "Length", numeric: true, width: 100 },
  { key: "gc_pct", label: "GC %", numeric: true, width: 80 },
  { key: "anchor", label: "Placed", width: 120 },
  { key: "anchor_seqid", label: "Reference sequence", width: 160 },
  { key: "anchor_start", label: "Reference position", numeric: true, width: 140 },
  { key: "anchor_end", label: "Reference end", numeric: true, width: 130 },
  { key: "region_seq", label: "Sequence", width: 230 },
  { key: "gene_names", label: "Genes inside (named)", width: 280 },
  { key: "n_orfs_complete", label: "Genes predicted", numeric: true, width: 130 },
  { key: "n_orfs", label: "Genes incl. partial", numeric: true, width: 150 },
];

const PANEL_COLUMNS: Column[] = [
  { key: "gene_id", label: "Gene", width: 180 },
  { key: "qlen", label: "Length", numeric: true, width: 100 },
  { key: "cov_pct", label: "Coverage %", numeric: true, width: 110 },
  { key: "identity", label: "Identity %", numeric: true, width: 110 },
  { key: "best_evalue", label: "Best match significance", numeric: true, width: 140 },
  { key: "call", label: "Call", width: 120 },
];

// The backend caps page_size at 1000 per request; to show the whole table
// (no pagination UI) we fetch every page at this size and concatenate them.
const FETCH_PAGE_SIZE = 1000;
// Header row height; the virtualizer offsets the first row by it because the
// header shares the rows' scroll container (so both scroll sideways together).
const HEADER_HEIGHT = 44;
const COL_WIDTHS_KEY = "straincompass-col-widths";
const HIDDEN_COLS_KEY = "straincompass-hidden-cols";

// Columns hidden by default per table, until the user changes it via the
// Columns picker (then their choice is remembered instead).
const DEFAULT_HIDDEN_COLS: Partial<Record<TableKind, string[]>> = {
  genes_coverage: ["start", "end", "length", "cov_bp"],
  gained: ["end", "anchor_end", "n_orfs"],
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
    (initialTable === "panel_recheck" && !run.has_panel) ||
    (initialTable === "gained" && !run.has_gained)
      ? "genes_coverage"
      : initialTable;
  const [table, setTable] = useState<TableKind>(safeInitialTable);
  const [queryId, setQueryId] = useState<number | undefined>(initialQueryId);
  const [search, setSearch] = useState("");
  const [debouncedSearch, setDebouncedSearch] = useState("");
  const [call, setCall] = useState<string>("");
  const [sortBy, setSortBy] = useState<string | null>(null);
  const [sortDir, setSortDir] = useState<"asc" | "desc">("asc");
  const [hiddenCols, setHiddenCols] = useState<Record<string, string[]>>(() => {
    try {
      return JSON.parse(localStorage.getItem(HIDDEN_COLS_KEY) ?? "{}");
    } catch {
      return {};
    }
  });
  /** The fetched page, tagged with the table it belongs to: right after a
   * tab switch the columns are already the new table's while the rows are
   * still the old one's, and rendering e.g. gc_pct over a genes row
   * crashed the page (the error boundary caught `toFixed of undefined`).
   * Rows that do not belong to the shown table are treated as absent
   * until their own fetch lands. */
  const [data, setData] = useState<(Page<AnyRow> & { table: TableKind }) | null>(null);
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
    if (table === "gained") return GAINED_COLUMNS;
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
      sort_by: sortBy ?? undefined,
      sort_dir: sortBy ? sortDir : undefined,
      search: debouncedSearch || undefined,
      call: call || undefined,
    }),
    [queryId, sortBy, sortDir, debouncedSearch, call],
  );

  useEffect(() => {
    setSortBy(null);
    setCall("");
    setPinnedGene(null);
  }, [table]);

  useEffect(() => {
    let cancelled = false;
    setLoading(true);
    setError(null);

    function fetchPage(page: number) {
      const q: TableQuery = { ...query, page, page_size: FETCH_PAGE_SIZE };
      if (table === "genes_coverage") return api.genesCoverage(run.id, q);
      if (table === "unaligned_gaps") return api.unalignedGaps(run.id, q);
      if (table === "gained") return api.gained(run.id, q);
      if (table === "panel_recheck") return api.panelRecheck(run.id, q);
      return api.matrix(run.id, q);
    }

    // The backend paginates at up to FETCH_PAGE_SIZE rows per request; fetch
    // every page and concatenate so the whole filtered/sorted table renders
    // in one virtualized scroll instead of behind Previous/Next.
    async function fetchAll() {
      const rows: AnyRow[] = [];
      let total = 0;
      let page = 0;
      for (;;) {
        const d = await fetchPage(page);
        if (cancelled) return;
        rows.push(...d.rows);
        total = d.total;
        if (d.rows.length === 0 || rows.length >= total) break;
        page++;
      }
      if (!cancelled) setData({ table, rows, total });
    }

    fetchAll()
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

  const total = data && data.table === table ? data.total : 0;

  // The gained table's copy-to-clipboard column: the sequences of all
  // of this query's regions, fetched once per view and kept client-side,
  // so a click copies synchronously from memory instead of racing the
  // clipboard's user-gesture window against the network.
  const [gainedSeqs, setGainedSeqs] = useState<Map<string, string> | null>(null);
  useEffect(() => {
    if (table !== "gained" || queryId === undefined) {
      setGainedSeqs(null);
      return;
    }
    let cancelled = false;
    setGainedSeqs(null);
    api
      .gainedSequences(run.id, queryId)
      .then((m) => {
        if (!cancelled) setGainedSeqs(m);
      })
      .catch(() => {});
    return () => {
      cancelled = true;
    };
  }, [run.id, table, queryId]);

  // The pin slot holds a locus tag for the gene tables and a synthetic
  // "contig:start" for gained regions, which have no locus tag of their own.
  const pinnedGainedRow = useMemo(() => {
    if (table !== "gained" || !pinnedGene || data?.table !== "gained") {
      return null;
    }
    return (
      (data.rows.find(
        (r) => gainedKey(r as GainedRow) === pinnedGene,
      ) as GainedRow | undefined) ?? null
    );
  }, [table, pinnedGene, data]);

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
              ["gained", "Gained"],
              ["panel_recheck", "Panel recheck"],
              ["matrix", "Presence / absence"],
            ] as [TableKind, string][]
          )
            .filter(
              ([k]) =>
                (k !== "panel_recheck" || run.has_panel) &&
                (k !== "gained" || run.has_gained),
            )
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

        {(table === "genes_coverage" ||
          table === "panel_recheck" ||
          table === "gained") &&
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

        {table === "gained" && (
          <div className="flex rounded-lg border border-zinc-300 overflow-hidden h-11 dark:border-zinc-700">
            {[
              ["", "All"],
              ["anchored", "Placed on the reference"],
              ["unanchored", "Not placed"],
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

        <ExportButton run={run} table={table} query={query} />

        {loading && <Spinner className="text-zinc-400" />}
      </div>

      {error && (
        <p className="text-red-700 text-[15px] py-2 dark:text-red-400">{error}</p>
      )}

      {/* table */}
      <div className="border border-zinc-200 rounded-xl bg-white overflow-hidden dark:border-zinc-800 dark:bg-zinc-900">
        <VirtualTable
          columns={visibleColumns}
          rows={data && data.table === table ? data.rows : []}
          table={table}
          sortBy={sortBy}
          sortDir={sortDir}
          onSort={toggleSort}
          onPinGene={setPinnedGene}
          onOpenGene={onOpenGene}
          run={run}
          colWidths={colWidths}
          onResizeColumn={resizeColumn}
          gainedSeqs={table === "gained" ? gainedSeqs : null}
        />
        {/* row count */}
        <div className="flex items-center px-4 h-12 border-t border-zinc-200 text-sm text-zinc-500 dark:border-zinc-800 dark:text-zinc-400">
          <span>
            {total === 0
              ? "No rows match the current filters"
              : `Showing all ${total.toLocaleString("en-US")} row${total === 1 ? "" : "s"}`}
          </span>
        </div>
      </div>

      {/* right-click preview, pinned until another row is right-clicked or this is closed */}
      {pinnedGene && (table === "genes_coverage" || table === "matrix") && (
        <GenePreview
          runId={run.id}
          locus={pinnedGene}
          onOpen={() => onOpenGene(pinnedGene)}
          onClose={() => setPinnedGene(null)}
        />
      )}
      {pinnedGene && table === "gained" && pinnedGainedRow && (
        <GainedOrfsCard
          row={pinnedGainedRow}
          runId={run.id}
          queryId={queryId ?? run.queries[0].file_id}
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
  gainedSeqs,
}: {
  columns: Column[];
  rows: AnyRow[];
  table: TableKind;
  sortBy: string | null;
  sortDir: "asc" | "desc";
  onSort: (key: string) => void;
  onPinGene: React.Dispatch<React.SetStateAction<string | null>>;
  onOpenGene: (locus: string) => void;
  run: Run;
  colWidths: Record<string, number>;
  onResizeColumn: (key: string, width: number) => void;
  gainedSeqs: Map<string, string> | null;
}) {
  const parentRef = useRef<HTMLDivElement>(null);
  const virtualizer = useVirtualizer({
    count: rows.length,
    getScrollElement: () => parentRef.current,
    estimateSize: () => 44,
    overscan: 12,
    // The header row shares the scroll container with the rows, so the first
    // row starts one header height down.
    paddingStart: HEADER_HEIGHT,
  });

  function widthFor(c: Column) {
    return colWidths[`${table}:${c.key}`] ?? c.width ?? 120;
  }

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

  // The header and the rows live in one scroll container and are sized by the
  // same track, so they can never drift apart horizontally: the track is at
  // least as wide as the container (columns stretch to fill) and grows to the
  // summed column widths when those no longer fit (one shared scrollbar).
  const trackMinWidth = columns.reduce((sum, c) => sum + widthFor(c), 0);

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
    // One scroll container for the header and the rows: a single horizontal
    // scrollbar moves both together, and the header stays pinned vertically.
    <div
      ref={parentRef}
      className="overflow-auto thin-scroll"
      style={{ height: "min(62vh, 700px)" }}
    >
      <div
        style={{
          minWidth: trackMinWidth,
          height: virtualizer.getTotalSize(),
          position: "relative",
        }}
      >
        {/* header */}
        <div
          className="flex bg-zinc-50 border-b border-zinc-200 sticky top-0 z-10 dark:bg-zinc-900 dark:border-zinc-800"
          style={{ height: HEADER_HEIGHT, boxSizing: "border-box" }}
        >
          {columns.map((c) => (
            <div
              key={c.key}
              className="relative flex items-stretch border-r border-zinc-200 last:border-r-0 dark:border-zinc-800"
              style={flexStyleFor(c)}
            >
              <button
                onClick={() => onSort(c.key)}
                className={`flex-1 min-w-0 flex items-center gap-1 px-3 h-full text-left text-xs font-semibold uppercase tracking-wide text-zinc-500 hover:text-zinc-900 dark:text-zinc-400 dark:hover:text-zinc-100 ${c.numeric ? "justify-end" : ""}`}
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
                if (table === "genes_coverage" || table === "matrix") {
                  e.preventDefault();
                  const locus = row["locus_tag"] as string;
                  if (!locus) return;
                  // Right-clicking the row whose preview is already open closes
                  // it; any other row moves the preview over to that row.
                  onPinGene((prev) => (prev === locus ? null : locus));
                } else if (table === "gained") {
                  e.preventDefault();
                  const key = gainedKey(row as unknown as GainedRow);
                  onPinGene((prev) => (prev === key ? null : key));
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
                    gainedSeqs={gainedSeqs}
                  />
                </div>
              ))}
            </div>
          );
        })}
      </div>
    </div>
  );
}

function Cell({
  col,
  row,
  table,
  run,
  gainedSeqs,
}: {
  col: string;
  row: Record<string, unknown>;
  table: TableKind;
  run: Run;
  gainedSeqs: Map<string, string> | null;
}) {
  const v = row[col];
  const [copied, setCopied] = useState(false);
  if (col === "call") return <CallBadge call={v as Call} />;
  if (table === "gained") {
    if (col === "anchor") return <AnchorBadge row={row as unknown as GainedRow} />;
    if (col === "gc_pct")
      return (
        <span>
          {typeof v === "number" ? v.toFixed(1) : "-"}
        </span>
      );
    if (col === "n_orfs" || col === "n_orfs_complete") {
      // null is "the gene finder did not run", which is not a count of
      // zero and must never be rendered as one.
      if (v === null || v === undefined)
        return (
          <span
            className="text-zinc-400 text-sm dark:text-zinc-500"
            title="The gene finder was not available for this run, so the genes inside this region were not predicted."
          >
            not available
          </span>
        );
      return <span>{(v as number).toLocaleString("en-US")}</span>;
    }
  }
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
  if (col === "seqid" && v && typeof row["start"] === "number" && typeof row["end"] === "number") {
    return (
      <a
        href={nuccoreRangeUrl(v as string, row["start"] as number, row["end"] as number)}
        target="_blank"
        rel="noreferrer"
        title={`View ${v} on NCBI, zoomed to this range`}
        className="truncate underline text-zinc-700 hover:text-zinc-900 dark:text-zinc-300 dark:hover:text-zinc-100"
        onClick={(e) => e.stopPropagation()}
      >
        {v as string}
      </a>
    );
  }
  if (col === "region_seq" && table === "gained") {
    const key = `${row["qry_seqid"]}:${row["start"]}-${row["end"]}`;
    const seq = gainedSeqs?.get(key);
    if (!seq) {
      // Not loaded (yet): the fetch is one request per table view.
      return (
        <span className="font-mono text-xs text-zinc-300 dark:text-zinc-600">loading...</span>
      );
    }
    return (
      <button
        title="Click to copy the whole sequence"
        onClick={(e) => {
          e.stopPropagation();
          navigator.clipboard?.writeText(seq).then(
            () => {
              setCopied(true);
              window.setTimeout(() => setCopied(false), 1200);
            },
            () => {},
          );
        }}
        className={`font-mono text-xs tracking-tight truncate cursor-pointer rounded px-1 -mx-1 ${
          copied
            ? "text-teal-700 dark:text-teal-400 bg-teal-50 dark:bg-teal-900/30"
            : "text-zinc-500 hover:text-zinc-900 hover:bg-zinc-100 dark:text-zinc-400 dark:hover:text-zinc-100 dark:hover:bg-zinc-800"
        }`}
      >
        {copied ? "copied ✓" : `${seq.slice(0, 28)}…`}
      </button>
    );
  }
  if (col === "gene_names") {
    const names = v as string[];
    if (!names || names.length === 0) {
      if (row["n_orfs"] === null || row["n_orfs"] === undefined)
        return <span className="text-zinc-300 dark:text-zinc-700">-</span>;
      if (!row["named"]) {
        // The naming pass never ran for this result: empty here means
        // unknown, not unmatched, and must not be read as "novel".
        return (
          <span
            className="text-zinc-400 dark:text-zinc-500"
            title="This run was computed before gene naming existed (or the naming search failed). Re-run the comparison to name the genes."
          >
            not available
          </span>
        );
      }
      return (
        <span
          className="text-zinc-400 dark:text-zinc-500"
          title="None of the predicted genes matched any protein of the reference."
        >
          all novel
        </span>
      );
    }
    return (
      <span className="truncate" title={names.join(", ")}>
        {names.slice(0, 6).join(", ")}
        {names.length > 6 && (
          <span className="text-zinc-400 dark:text-zinc-500"> +{names.length - 6} more</span>
        )}
        {(() => {
          // Named is a property of some genes, not of the region: the
          // ones left out are novel to the reference, and a column that
          // shows only the names invites reading them as the content.
          const nOrfs = (row["n_orfs"] as number | null) ?? 0;
          const novel = nOrfs - names.length;
          if (novel > 0)
            return (
              <span
                className="text-zinc-400 dark:text-zinc-500"
                title={`${novel} of the predicted genes match no protein of the reference - they are novel to it. Use the region card to search them at NCBI BLAST.`}
              >
                {" "}
                +{novel} novel
              </span>
            );
          return null;
        })()}
      </span>
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

/**
 * How confidently a gained region is placed. A region flanked on one side
 * only, or flanked inconsistently, must not read the same as one pinned
 * between two agreeing blocks - on a fragmented assembly most of them are.
 */
function AnchorBadge({ row }: { row: GainedRow }) {
  if (row.anchor === "unanchored")
    return (
      <span
        className="px-2 py-0.5 rounded text-xs bg-zinc-100 text-zinc-500 dark:bg-zinc-800 dark:text-zinc-400"
        title="This region is on a query contig with no alignment to the reference at all - usually a plasmid or a phage - so it cannot be placed."
      >
        not placed
      </span>
    );
  if (row.anchor === "flank" || row.flanks_disagree)
    return (
      <span
        className="px-2 py-0.5 rounded text-xs bg-amber-100 text-amber-800 dark:bg-amber-900/40 dark:text-amber-300"
        title={
          row.flanks_disagree
            ? "The alignments on either side point at different places, so this position is approximate."
            : "Only one side of this region aligns to the reference, so it is placed at that junction alone."
        }
      >
        junction
      </span>
    );
  return (
    <span
      className="px-2 py-0.5 rounded text-xs bg-teal-100 text-teal-800 dark:bg-teal-900/40 dark:text-teal-300"
      title="Both flanking alignments agree on where this region sits."
    >
      between
    </span>
  );
}

/** The genes predicted inside one gained region, pinned by right-click,
 * plus the reference back-check of the region itself. */
function GainedOrfsCard({
  row,
  runId,
  queryId,
  onClose,
}: {
  row: GainedRow;
  runId: number;
  queryId: number;
  onClose: () => void;
}) {
  const ref = usePopoverDismiss(true, onClose);
  const [verify, setVerify] = useState<GainedVerify | null>(null);
  const [checking, setChecking] = useState(false);
  const [verifyError, setVerifyError] = useState<string | null>(null);
  const [identify, setIdentify] = useState<GainedIdentify | null>(null);
  const [identifying, setIdentifying] = useState(false);
  const [identifyError, setIdentifyError] = useState<string | null>(null);

  function check() {
    setChecking(true);
    setVerifyError(null);
    api
      .gainedVerify(runId, queryId, row.qry_seqid, row.start, row.end)
      .then(setVerify)
      .catch((e) => setVerifyError((e as Error).message))
      .finally(() => setChecking(false));
  }

  function identifyGenes() {
    setIdentifying(true);
    setIdentifyError(null);
    api
      .gainedIdentify(runId, queryId, row.qry_seqid, row.start, row.end)
      .then(setIdentify)
      .catch((e) => setIdentifyError((e as Error).message))
      .finally(() => setIdentifying(false));
  }

  return (
    <div
      ref={ref}
      className="mt-3 rounded-xl border border-zinc-200 bg-white p-4 dark:border-zinc-800 dark:bg-zinc-900"
    >
      <div className="flex items-start justify-between gap-4">
        <div>
          <p className="font-semibold text-[15px]">
            {row.qry_seqid}:{row.start.toLocaleString("en-US")}-
            {row.end.toLocaleString("en-US")}
          </p>
          <p className="text-sm text-zinc-500 dark:text-zinc-400">
            {row.length.toLocaleString("en-US")} bp with no alignment to the
            reference, GC {row.gc_pct.toFixed(1)}%
            {row.at_contig_end && " - at a contig end"}
          </p>
        </div>
        <button
          onClick={onClose}
          className="text-sm text-zinc-500 hover:text-zinc-900 dark:hover:text-zinc-100"
        >
          Close
        </button>
      </div>
      {row.n_orfs === null ? (
        <p className="mt-3 text-sm text-zinc-500 dark:text-zinc-400">
          The gene finder was not available for this run, so the genes inside
          this region were not predicted.
        </p>
      ) : row.orfs.length === 0 ? (
        <p className="mt-3 text-sm text-zinc-500 dark:text-zinc-400">
          No genes were predicted inside this region.
        </p>
      ) : (
        <ul className="mt-3 space-y-1">
          {row.orfs.map((o, i) => {
            const id = identify?.orfs[i];
            // The run's naming pass already named the ORF: show the name
            // immediately. The on-demand search adds the sequence (for
            // the NCBI links) and covers results from before the pass.
            const m = id?.match ?? o.best ?? null;
            return (
              <li key={i} className="text-sm">
                <div className="font-mono tabular-nums">
                  {o.start.toLocaleString("en-US")}-
                  {o.end.toLocaleString("en-US")} ({o.strand < 0 ? "-" : "+"}){" "}
                  <span
                    className={
                      o.partial
                        ? "text-amber-700 dark:text-amber-400"
                        : "text-zinc-500 dark:text-zinc-400"
                    }
                  >
                    {o.partial ? "partial" : "complete"}
                  </span>{" "}
                  <span className="text-zinc-400 dark:text-zinc-500">
                    confidence {o.confidence.toFixed(1)}
                  </span>
                </div>
                {id && <IdentifiedLine o={id} />}
                {!id && m && <IdentifiedLine o={{ start: o.start, end: o.end, strand: o.strand, match: m, seq: "" }} />}
              </li>
            );
          })}
        </ul>
      )}
      {row.orfs.length > 0 && (
        <div className="mt-2">
          {identify ? (
            <p className="text-xs text-zinc-400 dark:text-zinc-500">
              {identify.orfs.filter((o) => o.match).length} of{" "}
              {identify.orfs.length} gene
              {identify.orfs.length === 1 ? "" : "s"} named by similarity to the
              reference's proteins.{" "}
              {identify.region_seq.length <= 20000 ? (
                <>
                  Unnamed ones are novel to this reference -{" "}
                  <a
                    href={ncbiBlastUrl("blastn", identify.region_seq) ?? undefined}
                    target="_blank"
                    rel="noreferrer"
                    className="underline hover:text-zinc-700 dark:hover:text-zinc-300"
                  >
                    search the whole region at NCBI BLAST
                  </a>{" "}
                  to identify them.
                </>
              ) : (
                <>Unnamed ones are novel to this reference.</>
              )}
            </p>
          ) : identifyError ? (
            <div>
              <p className="text-sm text-red-700 dark:text-red-400">
                {identifyError}
              </p>
              <button
                onClick={identifyGenes}
                className="mt-1 text-sm text-zinc-500 underline hover:text-zinc-900 dark:hover:text-zinc-100"
              >
                Try again
              </button>
            </div>
          ) : (
            <button
              onClick={identifyGenes}
              disabled={identifying}
              className="text-sm text-zinc-500 underline hover:text-zinc-900 disabled:opacity-50 dark:hover:text-zinc-100"
            >
              {identifying
                ? "Searching the reference's proteins..."
                : "Identify the genes against the reference's proteins"}
            </button>
          )}
        </div>
      )}
      <div className="mt-4 pt-3 border-t border-zinc-100 dark:border-zinc-800">
        {verify ? (
          <GainedVerifyResult v={verify} regionLength={row.length} />
        ) : verifyError ? (
          <div>
            <p className="text-sm text-red-700 dark:text-red-400">{verifyError}</p>
            <button
              onClick={check}
              className="mt-2 text-sm text-zinc-500 underline hover:text-zinc-900 dark:hover:text-zinc-100"
            >
              Try again
            </button>
          </div>
        ) : (
          <div className="flex flex-wrap items-center gap-3">
            <button
              onClick={check}
              disabled={checking}
              className="h-9 px-3 rounded-lg bg-zinc-900 text-white text-sm font-medium hover:bg-zinc-700 disabled:opacity-50 dark:bg-zinc-100 dark:text-zinc-900 dark:hover:bg-zinc-300"
            >
              {checking ? "Searching the reference..." : "Check against the reference"}
            </button>
            {checking && <Spinner className="text-zinc-400" />}
            {!checking && (
              <p className="text-xs text-zinc-400 dark:text-zinc-500">
                Searches the region's sequence against the whole reference
                genome with a sensitive blastn, the back-check of "no
                alignment".
              </p>
            )}
          </div>
        )}
      </div>
    </div>
  );
}

/** An NCBI BLAST submission URL with the sequence embedded, or null when
 * the sequence is too long for a URL (the site would truncate it
 * silently). */
function ncbiBlastUrl(
  program: "blastn" | "blastx",
  seq: string,
): string | null {
  if (seq.length > 20000) return null;
  return (
    `https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=${program}` +
    `&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome&BLAST_DATABASE=nr` +
    `&QUERY=${encodeURIComponent(seq)}`
  );
}

/** The name (or novelness) of one predicted gene of a gained region,
 * from the identification search. */
function IdentifiedLine({ o }: { o: IdentifiedOrf }) {
  const m = o.match;
  if (!m) {
    if (!o.seq) return null;
    const blast = ncbiBlastUrl("blastx", o.seq);
    return (
      <p className="ml-4 text-xs text-zinc-500 dark:text-zinc-400">
        no similar gene in the reference
        {blast && (
          <>
            {" - "}
            <a
              href={blast}
              target="_blank"
              rel="noreferrer"
              className="underline hover:text-zinc-700 dark:hover:text-zinc-300"
            >
              search this gene at NCBI BLAST
            </a>
          </>
        )}
      </p>
    );
  }
  return (
    <p className="ml-4 text-xs">
      <span className="text-zinc-700 dark:text-zinc-300">
        {m.protein_id ? (
          <a
            href={`https://www.ncbi.nlm.nih.gov/protein/${encodeURIComponent(m.protein_id)}`}
            target="_blank"
            rel="noreferrer"
            className="underline hover:text-zinc-900 dark:hover:text-zinc-100"
          >
            {m.label || m.locus_tag}
          </a>
        ) : (
          m.label || m.locus_tag
        )}
      </span>{" "}
      <span className="text-zinc-400 dark:text-zinc-500">
        {m.identity.toFixed(0)}% aa identity over {m.coverage.toFixed(0)}% of
        the gene (E {fmtE(m.evalue)})
      </span>
    </p>
  );
}

/** Share of a gained region covered by any hit, in percent. Hit intervals
 * are 1-based inclusive on the region. */
function unionCoverage(hits: GainedBlastHit[], regionLength: number): number {
  const ivs = hits
    .map((h) => [h.qry_start, h.qry_end] as [number, number])
    .sort((a, b) => a[0] - b[0]);
  let covered = 0;
  let last = 0;
  for (const [s, e] of ivs) {
    if (e > last) {
      covered += e - Math.max(s - 1, last);
      last = e;
    }
  }
  return (100 * covered) / Math.max(regionLength, 1);
}

function fmtE(e: number): string {
  if (e === 0) return "0.0";
  if (e >= 0.001) return e.toFixed(2);
  return e.toExponential(1);
}

/** The hit rows shared by every tier: where on the reference, how
 * similar, how long, how significant. */
function HitGrid({ hits, unit }: { hits: GainedBlastHit[]; unit: "bp" | "aa" }) {
  return (
    <div className="mt-2 grid grid-cols-[minmax(0,1fr)_auto_auto_auto] gap-x-4 text-sm">
      {hits.slice(0, 8).map((h, i) => (
        <div key={i} className="contents">
          <span className="font-mono tabular-nums truncate">
            {h.ref_seqid}:{h.ref_start.toLocaleString("en-US")}-
            {h.ref_end.toLocaleString("en-US")}
          </span>
          <span
            className="text-right font-mono tabular-nums"
            title={unit === "aa" ? "amino-acid identity" : "nucleotide identity"}
          >
            {h.identity.toFixed(1)}%
          </span>
          <span className="text-right font-mono tabular-nums">
            {h.length.toLocaleString("en-US")} {unit}
          </span>
          <span
            className="text-right font-mono tabular-nums text-zinc-400 dark:text-zinc-500"
            title={`bitscore ${h.bitscore}`}
          >
            {fmtE(h.evalue)}
          </span>
          {h.genes.length > 0 && (
            <span className="col-span-4 -mt-1 text-xs text-zinc-500 dark:text-zinc-400">
              {h.genes.join(", ")}
            </span>
          )}
        </div>
      ))}
      {hits.length > 8 && (
        <span className="col-span-4 text-xs text-zinc-400 dark:text-zinc-500">
          +{hits.length - 8} more hit{hits.length - 8 === 1 ? "" : "s"}
        </span>
      )}
    </div>
  );
}

/** The weak tier, collapsed: what the loose search saw below the
 * threshold the verdicts are built on. */
function WeakTier({ hits }: { hits: GainedBlastHit[] }) {
  return (
    <details className="mt-2">
      <summary className="text-xs text-zinc-400 dark:text-zinc-500 cursor-pointer">
        {hits.length} weak match{hits.length === 1 ? "" : "es"} (E &le; 10,
        mostly noise)
      </summary>
      <HitGrid hits={hits} unit="bp" />
    </details>
  );
}

/** The cutoff-free companion fact: the longest run of bases the region
 * shares verbatim with the reference, wherever it is. */
function LongestExactLine({ v }: { v: GainedVerify }) {
  if (v.longest_exact_bp === 0) {
    return (
      <p className="mt-2 text-xs text-zinc-400 dark:text-zinc-500">
        Longest exact match to the reference: the region shares no run of
        bases with the reference at all.
      </p>
    );
  }
  return (
    <p
      className="mt-2 text-xs text-zinc-400 dark:text-zinc-500"
      title="The longest run of bases, anywhere in the reference, that appears verbatim in this region. Unrelated DNA of these sizes shares about log4(region length x reference length) bases by chance alone - on a bacterial genome pair, roughly 15-20. A much longer run means shared sequence."
    >
      Longest exact match to the reference: {v.longest_exact_bp} bp at{" "}
      {v.longest_exact_seqid}:
      {v.longest_exact_start.toLocaleString("en-US")}-
      {v.longest_exact_end.toLocaleString("en-US")} (region{" "}
      {v.longest_exact_qry_start.toLocaleString("en-US")}-
      {v.longest_exact_qry_end.toLocaleString("en-US")})
    </p>
  );
}

/**
 * The verdict of the reference back-check, in descending order of what
 * it means: a strong nucleotide hit means the aligner could not use a
 * match the reference does carry (the caveat the gained table lives
 * with); a translated hit means a relative too diverged for nucleotide
 * comparison; weak-only is inconclusive; and nothing anywhere - the
 * confirmation the word "gained" wants, with the honest limit that no
 * sequence search proves absence outright.
 */
function GainedVerifyResult({
  v,
  regionLength,
}: {
  v: GainedVerify;
  regionLength: number;
}) {
  const covered = unionCoverage(v.hits, regionLength);
  const best = v.hits.reduce((m, h) => Math.max(m, h.identity), 0);
  const probablyPresent = covered >= 80 && best >= 95;
  const tx = v.tx_hits ?? [];
  const txBest = tx.length > 0 ? tx.reduce((m, h) => Math.max(m, h.identity), 0) : 0;

  let badge: React.ReactNode;
  let body: React.ReactNode;

  if (v.hits.length > 0) {
    badge = (
      <span
        className="px-2 py-0.5 rounded text-xs bg-amber-100 text-amber-800 dark:bg-amber-900/40 dark:text-amber-300"
        title={`${covered.toFixed(0)}% of the region matched the reference at up to ${best.toFixed(1)}% identity.`}
      >
        {probablyPresent
          ? "similar sequence in the reference"
          : "partial similarity in the reference"}
      </span>
    );
    body = (
      <>
        <p className="mt-2 text-sm text-zinc-500 dark:text-zinc-400">
          {probablyPresent
            ? `About ${covered.toFixed(0)}% of this region matched the reference at up to ${best.toFixed(1)}% identity. The whole-genome aligner anchors on matches unique to the reference side, so a copy of something the reference carries several times over can fail to align - this region is probably not a true gain.`
            : `About ${covered.toFixed(0)}% of the region matched at up to ${best.toFixed(1)}% identity. Short or divergent matches can be shared repeats, conserved domains or the remains of a longer gain - read the hits below before concluding.`}
        </p>
        <HitGrid hits={v.hits} unit="bp" />
        {v.weak_hits.length > 0 && <WeakTier hits={v.weak_hits} />}
      </>
    );
  } else if (tx.length > 0) {
    badge = (
      <span
        className="px-2 py-0.5 rounded text-xs bg-amber-100 text-amber-800 dark:bg-amber-900/40 dark:text-amber-300"
        title={`The translated search found amino-acid similarity up to ${txBest.toFixed(1)}% identity.`}
      >
        divergent coding homolog in the reference
      </span>
    );
    body = (
      <>
        <p className="mt-2 text-sm text-zinc-500 dark:text-zinc-400">
          The nucleotide search found nothing, but the translated search
          (tblastx) found amino-acid similarity up to {txBest.toFixed(1)}%
          identity: this region probably codes for a relative of something
          the reference carries, too diverged for nucleotide comparison to
          see. Probably not a true gain.
        </p>
        <HitGrid hits={tx} unit="aa" />
        {v.weak_hits.length > 0 && <WeakTier hits={v.weak_hits} />}
      </>
    );
  } else {
    // No strong nucleotide hit and no translated hit. What remains is
    // the weak tier, and it decides the tone: a genome-sized search
    // produces dozens of fragmentary matches below the threshold by
    // chance, so only a substantial one (80+ bases aligned, longer
    // than the noise floor) is worth an amber "no clear similarity" -
    // otherwise the fragments stay collapsed under the green verdict.
    const substantial = v.weak_hits.some((h) => h.length >= 80);
    const wCovered = unionCoverage(v.weak_hits, regionLength);
    const wBest = v.weak_hits.reduce((m, h) => Math.max(m, h.identity), 0);
    const txRan = v.tx_hits !== null;
    if (substantial) {
      badge = (
        <span
          className="px-2 py-0.5 rounded text-xs bg-amber-100 text-amber-800 dark:bg-amber-900/40 dark:text-amber-300"
          title="Only matches below the detection threshold were found."
        >
          no clear similarity in the reference
        </span>
      );
      body = (
        <>
          <p className="mt-2 text-sm text-zinc-500 dark:text-zinc-400">
            The only matches found are below the detection threshold (E
            between 1e-5 and 10): about {wCovered.toFixed(0)}% of the
            region matched at up to {wBest.toFixed(1)}% identity, including
            a stretch long enough that it may mean something.
            {txRan
              ? " The translated search found nothing."
              : ""}
          </p>
          <WeakTier hits={v.weak_hits} />
        </>
      );
    } else {
      badge = (
        <span
          className="px-2 py-0.5 rounded text-xs bg-teal-100 text-teal-800 dark:bg-teal-900/40 dark:text-teal-300"
          title={
            txRan
              ? "A sensitive nucleotide search and a translated (amino-acid) search of this region both found nothing similar anywhere in the reference genome."
              : "A sensitive nucleotide search of this region found nothing similar anywhere in the reference genome."
          }
        >
          {txRan
            ? "not found in the reference, even translated"
            : "not found in the reference"}
        </span>
      );
      body = (
        <>
          <p className="mt-2 text-sm text-zinc-500 dark:text-zinc-400">
            {txRan
              ? "Sensitive nucleotide and translated searches of this region found no similar sequence anywhere in the reference genome, at either the nucleotide or the amino-acid level. This is the closest available evidence of a true gain - with the standing limit that no sequence search can prove absence outright."
              : "A sensitive nucleotide search of this region found no similar sequence anywhere in the reference genome. This is the back-check the alignment alone cannot give, and it is consistent with the region being truly gained."}
          </p>
          {v.tx_note && (
            <p className="mt-1 text-xs text-zinc-400 dark:text-zinc-500">
              {v.tx_note}
            </p>
          )}
          {v.weak_hits.length > 0 && <WeakTier hits={v.weak_hits} />}
        </>
      );
    }
  }

  return (
    <div>
      {badge}
      {body}
      <LongestExactLine v={v} />
    </div>
  );
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
  const [detail, setDetail] = useState<GeneDetail | null>(null);
  useEffect(() => {
    let cancelled = false;
    api
      .geneDetail(runId, locus)
      .then((d) => {
        if (!cancelled) setDetail(d);
      })
      .catch(() => {});
    return () => {
      cancelled = true;
    };
  }, [runId, locus]);
  const ref = usePopoverDismiss(true, onClose);
  if (!detail) return null;
  return (
    <div
      ref={ref}
      className="fixed right-4 top-56 z-40 w-72 max-h-[calc(100vh-15.5rem)] flex flex-col bg-white border border-zinc-200 rounded-xl shadow-xl p-4 pointer-events-auto overflow-y-auto thin-scroll dark:bg-zinc-900 dark:border-zinc-800"
    >
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
        <a
          href={nuccoreRangeUrl(detail.seqid, detail.start, detail.end)}
          target="_blank"
          rel="noreferrer"
          className="underline hover:text-zinc-900 dark:hover:text-zinc-100"
        >
          {detail.seqid}:{detail.start.toLocaleString("en-US")}-
          {detail.end.toLocaleString("en-US")}
        </a>
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

