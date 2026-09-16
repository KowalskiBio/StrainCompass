import { useCallback, useEffect, useMemo, useRef, useState } from "react";
import { api } from "../api";
import type {
  AlignmentData,
  Call,
  RefseqWindow,
  Run,
  WgaData,
  WgaGene,
} from "../types";
import { Spinner } from "./ui";
import { useWheelGestures } from "./useWheelGestures";
import { clampRange, panRange, zoomRange } from "./genomeRange";
import type { Range } from "./genomeRange";

/** Show SNP/indel markers once the visible window is below this span. */
const VARIANT_SPAN = 20000;
/** A gene rectangle gets its label once it is this wide in pixels. */
const LABEL_MIN_W = 42;
/** Width used until the container has been measured. */
const FALLBACK_WIDTH = 1100;
/**
 * Draw the reference bases inside the genes once the window is this small,
 * with hysteresis so a zoom that hovers the boundary does not flicker.
 * Same thresholds the alignment view uses for its letters mode.
 */
const BASES_SPAN = 100;
const BASES_HYSTERESIS = 15;
/** The refseq endpoint caps a request at this many bases. */
const REFSEQ_MAX_WINDOW = 8192;
/**
 * Above this many bases per pixel, point variants are shaded by density rather
 * than drawn one mark each: below it every mark is at least a pixel wide on its
 * own, so individual marks are still close to faithful and stay hoverable.
 */
const DENSITY_BP_PER_PX = 2;
/** Floor opacity for a density column, so a lone mismatch is still visible. */
const MIN_DENSITY_INK = 0.15;
/** Variant marker colors, shared with the Alignment view (Oligool palette). */
const SNP_COLOR = "#dc2626";
const INS_COLOR = "#3b82f6";
const DEL_COLOR = "#9333ea";

/**
 * Strain map: the reference genome as one continuous line of gene
 * rectangles (beads on a string), colored by their presence call in
 * the selected query (green present, yellow partial, gray absent) or
 * by biotype when "reference annotation" is selected. Zooming in shows
 * gene names inside the boxes and, below a 20 kb window, SNP/indel
 * markers of the selected query. Hovering a gene shows its name,
 * position and function annotation.
 */
export function GenomeView({
  run,
  initialRange,
  initialGene,
  initialQuery,
  onOpenGene,
  onRangeChange,
  onQueryChange,
}: {
  run: Run;
  initialRange?: { seqid: string; start: number; end: number } | null;
  initialGene?: string | null;
  /** Query file id whose calls drive the gene colors; 0 = reference mode. */
  initialQuery?: number | null;
  onOpenGene: (locus: string) => void;
  onRangeChange: (r: { seqid: string; start: number; end: number }) => void;
  onQueryChange?: (queryId: number | null) => void;
}) {
  const [data, setData] = useState<WgaData | null>(null);
  const [error, setError] = useState<string | null>(null);
  const [seqid, setSeqid] = useState(initialRange?.seqid ?? "");
  const [range, setRange] = useState<Range | null>(
    initialRange ? { start: initialRange.start, end: initialRange.end } : null,
  );
  const svgRef = useRef<SVGSVGElement>(null);
  /** The map container, measured for width and owning the wheel gestures. */
  const [mapEl, setMapEl] = useState<HTMLDivElement | null>(null);
  const [measuredW, setMeasuredW] = useState(0);
  /**
   * Zoom re-projects every gene rather than transforming a layer, so a burst of
   * wheel events would otherwise re-render the whole map several times a frame.
   * Gestures accumulate here and are committed once per frame.
   */
  const pendingRangeRef = useRef<Range | null>(null);
  const rangeRafRef = useRef(0);
  /** The exact range object our last frame committed, to tell it apart from
   * ranges set elsewhere in the component. */
  const lastCommittedRef = useRef<Range | null>(null);
  const initialGeneRef = useRef(initialGene);
  const [popup, setPopup] = useState<{
    x: number;
    y: number;
    gene?: WgaGene;
    gi?: number;
    variant?: VariantPopupInfo;
  } | null>(null);
  const dragRef = useRef<{ x: number; start: number; moved: boolean } | null>(null);
  /** Set on mouseup after a drag, so the click that follows is ignored. */
  const wasDragRef = useRef(false);
  /** Query file id whose calls drive the gene colors; null = biotype. */
  const [colorBy, setColorBy] = useState<number | null>(
    initialQuery === undefined ? null : initialQuery === 0 ? null : initialQuery,
  );
  const [colorByInit, setColorByInit] = useState(initialQuery !== undefined);
  /** Index into data.genes of the hovered gene. */
  const [hover, setHover] = useState<number | null>(null);
  /** Index into `markers` of the hovered variant. */
  const [hoverVariant, setHoverVariant] = useState<number | null>(null);

  // Variant events of the whole run, fetched lazily on first deep zoom.
  const [alignment, setAlignment] = useState<AlignmentData | null>(null);
  const [alignmentPending, setAlignmentPending] = useState(false);
  const [variantError, setVariantError] = useState<string | null>(null);
  const alignmentFetched = useRef(false);
  /** Master switch for the variant layer. Off also skips fetching it at all. */
  const [variantsOn, setVariantsOn] = useState(true);
  /** Which kinds of variant are drawn, toggled from the legend. */
  const [variantKinds, setVariantKinds] = useState({
    snp: true,
    ins: true,
    del: true,
  });

  useEffect(() => {
    let cancelled = false;
    alignmentFetched.current = false;
    setAlignment(null);
    setVariantError(null);
    api
      .wga(run.id)
      .then((d) => {
        if (cancelled) return;
        setData(d);
        setSeqid((prev) => prev || (d.reference[0]?.[0] ?? ""));
        if (!colorByInit) {
          setColorBy(d.queries[0]?.query_id ?? null);
          setColorByInit(true);
        }
      })
      .catch((e) => !cancelled && setError((e as Error).message));
    return () => {
      cancelled = true;
    };
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [run.id]);

  const seqLength = useMemo(
    () => data?.reference.find((r) => r[0] === seqid)?.[1] ?? 0,
    [data, seqid],
  );

  /** The query whose presence calls drive the gene colors (null = biotype). */
  const selectedQuery = useMemo(
    () => data?.queries.find((q) => q.query_id === colorBy) ?? null,
    [data, colorBy],
  );

  /**
   * Genes of the shown contig, by ascending start, keeping each one's index into
   * data.genes (the per-query calls arrays are aligned with it). Rendering used
   * to walk all ~2900 genes on every zoom frame; this lets it bisect to the
   * visible slice instead.
   */
  const contigGenes = useMemo(() => {
    const out: { g: WgaGene; gi: number }[] = [];
    data?.genes.forEach((g, gi) => {
      if (g.seqid === seqid) out.push({ g, gi });
    });
    out.sort((a, b) => a.g.start - b.g.start);
    return out;
  }, [data, seqid]);

  /** Longest gene on the contig: how far back a gene can start and still overlap. */
  const maxGeneLen = useMemo(
    () => contigGenes.reduce((m, x) => Math.max(m, x.g.end - x.g.start), 0),
    [contigGenes],
  );

  const callCounts = useMemo(() => {
    if (!selectedQuery?.calls) return null;
    const counts = { PRESENT: 0, PARTIAL: 0, ABSENT: 0 } as Record<Call, number>;
    for (const c of selectedQuery.calls) counts[c]++;
    return counts;
  }, [selectedQuery]);

  useEffect(() => {
    if (seqLength > 0 && (!range || range.end > seqLength)) {
      setRange({ start: 1, end: seqLength });
    }
  }, [seqLength, range]);

  // publish the visible range for deep linking (debounced)
  useEffect(() => {
    if (!range || !seqid) return;
    const t = setTimeout(
      () => onRangeChange({ seqid, start: range.start, end: range.end }),
      600,
    );
    return () => clearTimeout(t);
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [range, seqid]);

  useEffect(() => {
    const gene = initialGeneRef.current;
    if (gene && data && !initialRange) {
      const g = data.genes.find((x) => x.locus_tag === gene);
      if (g) {
        setSeqid(g.seqid);
        const pad = Math.max(2000, (g.end - g.start) * 2);
        setRange({ start: Math.max(1, g.start - pad), end: g.end + pad });
      }
    }
  }, [initialRange, data]);

  // Fetch the run's variant events the first time the user zooms deep
  // enough to see them (computing them can take a moment on old runs).
  const span = range ? range.end - range.start : Infinity;
  // Deliberately a boolean, not the span itself: keying this effect on the
  // range re-ran it on every zoom and pan, and each re-run's cleanup cancelled
  // the in-flight request's state updates while the ref guard stopped it
  // starting a new one, so one gesture mid-fetch hung the spinner for good.
  const wantVariants = Boolean(data) && span < VARIANT_SPAN && variantsOn;
  useEffect(() => {
    if (!wantVariants || alignmentFetched.current) return;
    alignmentFetched.current = true;
    setVariantError(null);
    setAlignmentPending(true);
    api
      .alignment(run.id)
      .then(setAlignment)
      .catch((e: Error) => {
        // Let a later zoom try again rather than silently showing no markers.
        alignmentFetched.current = false;
        setVariantError(e.message);
      })
      .finally(() => setAlignmentPending(false));
  }, [wantVariants, run.id]);

  /** Variant events of the selected query on the visible seqid. */
  const variantEvents = useMemo(() => {
    if (!alignment || !colorBy) return null;
    const q = alignment.queries.find((x) => x.query_id === colorBy);
    const ev = q?.events[seqid];
    if (!ev) return null;
    return ev;
  }, [alignment, colorBy, seqid]);

  const showMarkers = variantsOn && Boolean(variantEvents) && span < VARIANT_SPAN;

  // Reference bases, drawn inside the genes at deep zoom.
  const [showBases, setShowBases] = useState(false);
  useEffect(() => {
    if (!showBases && span < BASES_SPAN - BASES_HYSTERESIS) setShowBases(true);
    else if (showBases && span > BASES_SPAN + BASES_HYSTERESIS) setShowBases(false);
  }, [span, showBases]);

  /** Fetched base windows per contig; panning reuses whatever already covers. */
  const refWindowsRef = useRef(new Map<string, RefseqWindow[]>());
  const [refseqVersion, setRefseqVersion] = useState(0);

  const baseFrom = range ? Math.max(1, Math.floor(range.start)) : 0;
  const baseTo = range ? Math.min(seqLength, Math.ceil(range.end)) : 0;

  useEffect(() => {
    if (!showBases || !seqid || !seqLength || baseTo < baseFrom) return;
    const covered = refWindowsRef.current
      .get(seqid)
      ?.some((w) => w.start <= baseFrom && w.end >= baseTo);
    if (covered) return;
    // Debounced: a fast zoom crosses many ranges before settling.
    const t = setTimeout(() => {
      // Snap to a 256 bp grid so panning keeps hitting the same window.
      const ws = Math.max(1, Math.floor((baseFrom - 256) / 256) * 256 + 1);
      const we = Math.min(seqLength, ws + REFSEQ_MAX_WINDOW - 1);
      api
        .refseq(run.id, seqid, ws, we)
        .then((w) => {
          const list = refWindowsRef.current.get(seqid) ?? [];
          list.push(w);
          refWindowsRef.current.set(seqid, list);
          setRefseqVersion((v) => v + 1);
        })
        .catch(() => {});
    }, 150);
    return () => clearTimeout(t);
  }, [showBases, seqid, seqLength, baseFrom, baseTo, run.id]);

  // Size the map to its container instead of a fixed width.
  useEffect(() => {
    if (!mapEl) return;
    const ro = new ResizeObserver(([entry]) => {
      setMeasuredW(Math.max(1, Math.round(entry.contentRect.width)));
    });
    ro.observe(mapEl);
    return () => ro.disconnect();
  }, [mapEl]);

  // Adopt ranges that came from the buttons, the range input or a contig switch,
  // so the next gesture builds on what is actually on screen. Deliberately NOT
  // our own commits: a wide view can take over 100 ms to render, and any gesture
  // landing in that window would be overwritten here by the range we just
  // committed, silently throwing away part of the zoom.
  useEffect(() => {
    if (range && range !== lastCommittedRef.current) pendingRangeRef.current = range;
  }, [range]);

  useEffect(() => () => cancelAnimationFrame(rangeRafRef.current), []);

  const applyRange = useCallback((r: Range) => {
    pendingRangeRef.current = r;
    if (rangeRafRef.current) return;
    rangeRafRef.current = requestAnimationFrame(() => {
      rangeRafRef.current = 0;
      lastCommittedRef.current = pendingRangeRef.current;
      if (pendingRangeRef.current) setRange(pendingRangeRef.current);
    });
  }, []);

  const setWheelEl = useWheelGestures<HTMLDivElement>((g) => {
    const base = pendingRangeRef.current;
    if (!base || !seqLength || !mapEl) return;
    const rect = mapEl.getBoundingClientRect();
    if (!rect.width) return;
    if (g.kind === "pan") {
      applyRange(
        panRange(base, (g.dx / rect.width) * (base.end - base.start), seqLength),
      );
      return;
    }
    const frac = (g.clientX - rect.left) / rect.width;
    const anchorBp = base.start + frac * (base.end - base.start);
    applyRange(zoomRange(base, g.factor, anchorBp, seqLength));
  });

  const setMapRef = useCallback(
    (n: HTMLDivElement | null) => {
      setMapEl(n);
      setWheelEl(n);
    },
    [setWheelEl],
  );

  if (error) return <p className="text-red-700 py-4 dark:text-red-400">{error}</p>;
  if (!data || !range || !seqid)
    return (
      <div className="flex items-center gap-3 text-zinc-500 py-16 justify-center dark:text-zinc-400">
        <Spinner /> Loading the genome view...
      </div>
    );

  const width = measuredW || FALLBACK_WIDTH;
  const rulerH = 28;
  const mapH = 168;
  const height = rulerH + mapH;
  const baselineY = rulerH + mapH / 2;
  const geneH = 46;
  const bpToX = (bp: number) =>
    ((bp - range.start) / Math.max(1, range.end - range.start)) * width;

  const zoom = (factor: number, anchorBp: number) => {
    const base: Range = pendingRangeRef.current ?? range;
    applyRange(zoomRange(base, factor, anchorBp, seqLength));
  };

  const bpAt = (clientX: number) => {
    const rect = svgRef.current!.getBoundingClientRect();
    const frac = (clientX - rect.left) / rect.width;
    return range.start + frac * (range.end - range.start);
  };

  // ticks for the ruler
  const tickStep = niceStep(span);
  const firstTick = Math.ceil(range.start / tickStep) * tickStep;
  const ticks: number[] = [];
  for (let t = firstTick; t <= range.end; t += tickStep) ticks.push(t);

  const inRange = (s: number, e: number) => e >= range.start && s <= range.end;

  // Only the genes that can touch the visible window. Bisect for the first gene
  // that could still overlap, then walk forward until past the right edge.
  const visibleGenes: { g: WgaGene; gi: number }[] = [];
  {
    const from = range.start - maxGeneLen;
    let lo = 0;
    let hi = contigGenes.length;
    while (lo < hi) {
      const mid = (lo + hi) >> 1;
      if (contigGenes[mid].g.start < from) lo = mid + 1;
      else hi = mid;
    }
    for (let i = lo; i < contigGenes.length; i++) {
      const e = contigGenes[i];
      if (e.g.start > range.end) break;
      if (inRange(e.g.start, e.g.end)) visibleGenes.push(e);
    }
  }

  /** Index into `markers` for the variant under an event target, if any. */
  const markerAt = (target: EventTarget | null): number | null => {
    const el = (target as Element | null)?.closest?.("[data-mi]");
    const mi = el ? Number(el.getAttribute("data-mi")) : NaN;
    return Number.isFinite(mi) ? mi : null;
  };

  /** The reference base at a 1-based position, if a fetched window covers it. */
  const baseAt = (pos: number): string | null => {
    for (const w of refWindowsRef.current.get(seqid) ?? []) {
      if (pos >= w.start && pos <= w.end) return w.seq[pos - w.start] ?? null;
    }
    return null;
  };

  // The visible bases, at most ~115 of them by the time this is on.
  const basePositions: number[] = [];
  if (showBases) {
    for (let p = baseFrom; p <= baseTo; p++) basePositions.push(p);
  }
  // Referenced so the letters re-render when a window arrives.
  void refseqVersion;

  /**
   * One delegated listener for the whole gene layer instead of three closures per
   * gene: at whole-genome view that was ~8600 new functions on every zoom frame.
   */
  const geneAt = (target: EventTarget | null): { g: WgaGene; gi: number } | null => {
    const el = (target as Element | null)?.closest?.("[data-gi]");
    const gi = el ? Number(el.getAttribute("data-gi")) : NaN;
    return Number.isFinite(gi) && data.genes[gi] ? { g: data.genes[gi], gi } : null;
  };

  // Variant markers of the selected query inside the visible range.
  const markers = showMarkers && variantEvents
    ? [
        ...(variantKinds.snp
          ? variantEvents.snps
              .filter((s) => s.pos >= range.start && s.pos <= range.end)
              .map((s) => ({ kind: "snp" as const, pos: s.pos, r: s.r, q: s.q }))
          : []),
        ...(variantKinds.del
          ? variantEvents.dels
              .filter((d) => d.pos <= range.end && d.pos + d.len - 1 >= range.start)
              .map((d) => ({ kind: "del" as const, pos: d.pos, len: d.len }))
          : []),
        ...(variantKinds.ins
          ? variantEvents.ins
              .filter((i) => i.pos >= range.start - 1 && i.pos <= range.end)
              .map((i) => ({ kind: "ins" as const, pos: i.pos, seq: i.seq }))
          : []),
      ].sort((a, b) => a.pos - b.pos)
    : [];

  const bpPerPx = span / Math.max(1, width);
  const markerY1 = baselineY - geneH / 2 - 5;
  const markerY2 = baselineY + geneH / 2 + 5;

  /**
   * Point variants bucketed into pixel columns, with the share of that column's
   * bases that differ. Only built when the map is too coarse to draw them
   * individually.
   */
  const densityColumns: { px: number; kind: "snp" | "ins"; frac: number }[] = [];
  if (showMarkers && bpPerPx > DENSITY_BP_PER_PX) {
    const counts = new Map<string, number>();
    for (const m of markers) {
      if (m.kind === "del") continue;
      const px = Math.floor(bpToX(m.kind === "ins" ? m.pos + 1 : m.pos));
      if (px < 0 || px > width) continue;
      const key = `${m.kind}:${px}`;
      counts.set(key, (counts.get(key) ?? 0) + 1);
    }
    for (const [key, n] of counts) {
      const [kind, px] = key.split(":");
      densityColumns.push({
        px: Number(px),
        kind: kind as "snp" | "ins",
        frac: n / bpPerPx,
      });
    }
  }

  /** Click in the gene area: open the nearest variant marker's popup. */
  function openVariantPopup(clientX: number) {
    if (markers.length === 0) return;
    const bp = bpAt(clientX);
    const bpPerPx = (range!.end - range!.start) / width;
    const tol = Math.max(4 * bpPerPx, 1);
    let best: (typeof markers)[number] | null = null;
    let bestDist = Infinity;
    for (const m of markers) {
      const d = Math.abs(m.pos - bp);
      if (d < bestDist) {
        bestDist = d;
        best = m;
      }
    }
    if (!best || bestDist > Math.max(tol, 8)) return;
    const rect = svgRef.current!.getBoundingClientRect();
    setPopup({
      x: clientX - rect.left,
      y: baselineY + geneH / 2,
      variant: { ...best, queryName: selectedQuery?.query_name ?? "" },
    });
  }

  return (
    <div className="space-y-3">
      {/* controls */}
      <div className="flex flex-wrap items-center gap-2 py-2">
        <select
          value={seqid}
          onChange={(e) => {
            setSeqid(e.target.value);
            const len = data.reference.find((r) => r[0] === e.target.value)?.[1] ?? 0;
            setRange({ start: 1, end: len });
          }}
          className="h-11 px-3 rounded-lg border border-zinc-300 bg-white text-[15px] dark:border-zinc-700 dark:bg-zinc-900"
          disabled={data.reference.length <= 1}
        >
          {data.reference.map(([id, len]) => (
            <option key={id} value={id}>
              {id} ({len.toLocaleString("en-US")} bp)
            </option>
          ))}
        </select>
        <label className="flex items-center gap-2 text-sm text-zinc-500 dark:text-zinc-400">
          Genes colored by
          <select
            value={colorBy ?? 0}
            onChange={(e) => {
              const v = Number(e.target.value);
              const next = v === 0 ? null : v;
              setColorBy(next);
              onQueryChange?.(next);
            }}
            className="h-11 px-3 rounded-lg border border-zinc-300 bg-white text-[15px] max-w-56 truncate dark:border-zinc-700 dark:bg-zinc-900"
          >
            <option value={0}>reference annotation</option>
            {data.queries.map((q) => (
              <option key={q.query_id} value={q.query_id}>
                {q.query_name}
              </option>
            ))}
          </select>
        </label>
        {selectedQuery && (
          <label
            className="flex items-center gap-2 h-11 px-3 rounded-lg border border-zinc-300 bg-white text-sm text-zinc-600 cursor-pointer select-none dark:border-zinc-700 dark:bg-zinc-900 dark:text-zinc-300"
            title={`Draw SNP and indel markers of ${selectedQuery.query_name} once the window is below ${(VARIANT_SPAN / 1000).toFixed(0)} kb. Turning this off also skips downloading them.`}
          >
            <input
              type="checkbox"
              checked={variantsOn}
              onChange={(e) => setVariantsOn(e.target.checked)}
              className="accent-zinc-900 dark:accent-zinc-100"
            />
            Variants
          </label>
        )}
        <RangeInput
          range={range}
          seqLength={seqLength}
          onChange={(r) => setRange(r)}
        />
        <div className="flex rounded-lg border border-zinc-300 overflow-hidden h-11 dark:border-zinc-700">
          <button
            className="px-3 bg-white hover:bg-zinc-100 text-[15px] dark:bg-zinc-900 dark:hover:bg-zinc-800"
            onClick={() => zoom(0.5, (range.start + range.end) / 2)}
            title="Zoom in"
          >
            +
          </button>
          <button
            className="px-3 bg-white hover:bg-zinc-100 text-[15px] border-l border-zinc-300 dark:bg-zinc-900 dark:hover:bg-zinc-800 dark:border-zinc-700"
            onClick={() => zoom(2, (range.start + range.end) / 2)}
            title="Zoom out"
          >
            &minus;
          </button>
          <button
            className="px-3 bg-white hover:bg-zinc-100 text-[15px] border-l border-zinc-300 dark:bg-zinc-900 dark:hover:bg-zinc-800 dark:border-zinc-700"
            onClick={() => setRange({ start: 1, end: seqLength })}
            title="Whole sequence"
          >
            Whole
          </button>
        </div>
        <span className="text-sm text-zinc-400 dark:text-zinc-500">
          {Math.round(range.start).toLocaleString("en-US")} -{" "}
          {Math.round(range.end).toLocaleString("en-US")} bp ({(
            (range.end - range.start) /
            1000
          ).toFixed(1)} kb shown)
        </span>
        {alignmentPending && (
          <span className="inline-flex items-center gap-2 text-sm text-zinc-400 dark:text-zinc-500">
            <Spinner /> loading variant markers...
          </span>
        )}
        {variantError && !alignmentPending && (
          <span className="text-sm text-red-700 dark:text-red-400">
            Variant markers failed to load ({variantError}). Zoom out and back in
            to retry.
          </span>
        )}
      </div>

      <div
        ref={setMapRef}
        /* overflow-x-clip, not overflow-hidden: the map already fits its box, and
           hiding both axes would cut off the hover tooltips below the line. */
        className="relative border border-zinc-200 rounded-xl bg-white overflow-x-clip touch-none overscroll-contain dark:border-zinc-800 dark:bg-zinc-900"
        onMouseLeave={() => {
          setPopup(null);
          setHover(null);
          setHoverVariant(null);
        }}
      >
        <svg
          ref={svgRef}
          width={width}
          height={height}
          viewBox={`0 0 ${width} ${height}`}
          className="block select-none"
          style={{ cursor: dragRef.current ? "grabbing" : "grab" }}
          onMouseDown={(e) => {
            dragRef.current = { x: e.clientX, start: range.start, moved: false };
          }}
          onMouseMove={(e) => {
            if (dragRef.current && svgRef.current) {
              if (Math.abs(e.clientX - dragRef.current.x) > 3) {
                dragRef.current.moved = true;
              }
              const rect = svgRef.current.getBoundingClientRect();
              const span = range.end - range.start;
              const shiftBp = ((dragRef.current.x - e.clientX) / rect.width) * span;
              applyRange(
                clampRange(dragRef.current.start + shiftBp, span, seqLength),
              );
            }
          }}
          onMouseUp={() => {
            wasDragRef.current = dragRef.current?.moved ?? false;
            dragRef.current = null;
          }}
          onMouseLeave={() => (dragRef.current = null)}
          onClick={(e) => {
            if (e.defaultPrevented) return;
            if (wasDragRef.current) {
              wasDragRef.current = false;
              return;
            }
            openVariantPopup(e.clientX);
          }}
        >
          {/* ruler */}
          <g>
            <line
              x1={0}
              x2={width}
              y1={rulerH - 8}
              y2={rulerH - 8}
              style={{ stroke: "var(--gv-ruler-line)" }}
            />
            {ticks.map((t) => (
              <g key={t}>
                <line
                  x1={bpToX(t)}
                  x2={bpToX(t)}
                  y1={rulerH - 14}
                  y2={rulerH - 8}
                  style={{ stroke: "var(--gv-tick)" }}
                />
                <text
                  x={bpToX(t)}
                  y={rulerH - 17}
                  fontSize="10"
                  style={{ fill: "var(--gv-tick-label)" }}
                  textAnchor="middle"
                  className="font-mono"
                >
                  {t >= 1000 ? `${(t / 1000).toFixed(tickStep >= 1000 ? 0 : 1)}k` : t}
                </text>
              </g>
            ))}
          </g>

          {/* the string: a thin baseline under all the gene beads */}
          <line
            x1={0}
            x2={width}
            y1={baselineY}
            y2={baselineY}
            style={{ stroke: "var(--gv-ruler-line)" }}
          />

          {/* genes: one continuous line of beads on the string */}
          <g
            onMouseOver={(e) => {
              setHover(geneAt(e.target)?.gi ?? null);
              // Markers sit on top of this layer, so reaching it means the
              // cursor has left any variant it was over.
              setHoverVariant(null);
            }}
            onMouseOut={() => setHover(null)}
            onClick={(e) => {
              const hit = geneAt(e.target);
              if (!hit) return;
              e.preventDefault();
              const rect = svgRef.current!.getBoundingClientRect();
              setPopup({
                x: e.clientX - rect.left,
                y: baselineY + geneH / 2,
                gene: hit.g,
                gi: hit.gi,
              });
            }}
          >
            {visibleGenes.map(({ g, gi }) => {
              const x = bpToX(Math.max(g.start, range.start));
              const x2 = bpToX(Math.min(g.end, range.end));
              const w = Math.max(2, x2 - x);
              const y = baselineY - geneH / 2;
              const fill = selectedQuery
                ? callColor(selectedQuery.calls?.[gi])
                : geneColor(g);
              const label = g.symbol || g.locus_tag;
              return (
                <g key={g.locus_tag} className="cursor-pointer" data-gi={gi}>
                  <GeneShape x={x} w={w} y={y} h={geneH} strand={g.strand} fill={fill} />
                  {w >= LABEL_MIN_W && (
                    <text
                      x={x + w / 2}
                      /* the bases take the middle line, so the name moves up */
                      y={showBases ? baselineY - geneH / 2 + 12 : baselineY + 3.5}
                      fontSize="10"
                      textAnchor="middle"
                      className="font-mono pointer-events-none select-none"
                      style={{
                        fill: "var(--gv-gene-label)",
                        stroke: "var(--gv-map-bg)",
                        strokeWidth: 2.5,
                        paintOrder: "stroke",
                      }}
                    >
                      {truncateLabel(label, w - 10)}
                    </text>
                  )}
                </g>
              );
            })}
          </g>

          {/*
            Variant markers. Reference base p occupies [p, p+1), the same span its
            letter is centred in, so a mismatch covers exactly its own base, a
            deletion covers the bases it removes, and an insertion sits on the
            boundary after its position rather than over a base.

            Once a base is well under a pixel those marks cannot be drawn
            faithfully: a minimum-width tick for every mismatch merges into solid
            colour and made a 3% divergent region look like a third of the map.
            Past DENSITY_BP_PER_PX the point events become one column per pixel,
            shaded by the share of bases in that column that actually differ.
          */}
          {showMarkers &&
            (bpPerPx > DENSITY_BP_PER_PX ? (
              <g className="pointer-events-none">
                {densityColumns.map(({ px, kind, frac }) => (
                  <rect
                    key={`${kind}${px}`}
                    x={px}
                    y={markerY1}
                    width={1}
                    height={markerY2 - markerY1}
                    fill={kind === "snp" ? SNP_COLOR : INS_COLOR}
                    fillOpacity={Math.min(1, Math.max(MIN_DENSITY_INK, frac))}
                  />
                ))}
                {markers.map((m, i) =>
                  m.kind === "del" ? (
                    <rect
                      key={i}
                      x={bpToX(m.pos)}
                      y={markerY1}
                      width={Math.max(1, bpToX(m.pos + m.len) - bpToX(m.pos))}
                      height={markerY2 - markerY1}
                      fill={DEL_COLOR}
                      fillOpacity={0.5}
                    />
                  ) : null,
                )}
              </g>
            ) : (
              <g
                onMouseOver={(e) => setHoverVariant(markerAt(e.target))}
                onMouseOut={() => setHoverVariant(null)}
              >
                {markers.map((m, i) => {
                  const color =
                    m.kind === "snp" ? SNP_COLOR : m.kind === "del" ? DEL_COLOR : INS_COLOR;
                  if (m.kind === "ins") {
                    const x = bpToX(m.pos + 1);
                    return (
                      <line
                        key={i}
                        data-mi={i}
                        x1={x}
                        x2={x}
                        y1={markerY1}
                        y2={markerY2}
                        stroke={color}
                        strokeWidth={2}
                        strokeLinecap="round"
                      />
                    );
                  }
                  const x = bpToX(m.pos);
                  const to = m.kind === "del" ? m.pos + m.len : m.pos + 1;
                  return (
                    <rect
                      key={i}
                      data-mi={i}
                      x={x}
                      y={markerY1}
                      width={Math.max(1, bpToX(to) - x)}
                      height={markerY2 - markerY1}
                      fill={color}
                      fillOpacity={m.kind === "del" ? 0.5 : showBases ? 0.35 : 1}
                    />
                  );
                })}
              </g>
            ))}

          {/* reference sequence, once a base is wide enough to read */}
          {showBases && (
            <g className="pointer-events-none">
              {basePositions.map((pos) => {
                const ch = baseAt(pos);
                if (!ch) return null;
                const cx = bpToX(pos + 0.5);
                return (
                  <text
                    key={pos}
                    x={cx}
                    y={baselineY + 4}
                    textAnchor="middle"
                    fontSize={Math.min(14, Math.max(8, (width / span) * 0.8))}
                    className="font-mono select-none"
                    style={{
                      fill: "var(--gv-gene-label)",
                      stroke: "var(--gv-map-bg)",
                      strokeWidth: 2.5,
                      paintOrder: "stroke",
                    }}
                  >
                    {ch}
                  </text>
                );
              })}
            </g>
          )}
        </svg>

        {hoverVariant !== null && markers[hoverVariant] && (
          <VariantTooltip
            variant={markers[hoverVariant]}
            x={bpToX(
              markers[hoverVariant].kind === "ins"
                ? markers[hoverVariant].pos + 1
                : markers[hoverVariant].pos + 0.5,
            )}
            y={baselineY + geneH / 2 + 6}
            queryName={selectedQuery?.query_name ?? ""}
            svgWidth={width}
          />
        )}

        {hover !== null && hoverVariant === null && data.genes[hover] && data.genes[hover].seqid === seqid && (
          <GeneTooltip
            gene={data.genes[hover]}
            x={bpToX(
              (Math.max(data.genes[hover].start, range.start) +
                Math.min(data.genes[hover].end, range.end)) /
                2,
            )}
            y={baselineY + geneH / 2}
            queryName={selectedQuery?.query_name ?? null}
            call={selectedQuery?.calls?.[hover]}
            covPct={selectedQuery?.cov_pcts?.[hover]}
            identity={selectedQuery?.identities?.[hover]}
            svgWidth={width}
          />
        )}
      </div>

      {/* legend */}
      <div className="flex flex-wrap items-center gap-4 text-xs text-zinc-500 dark:text-zinc-400">
        {selectedQuery ? (
          <>
            <span className="inline-flex items-center gap-1.5">
              <span className="w-3 h-3 rounded-sm" style={{ background: "var(--gv-call-present)" }} />
              present
              {callCounts && ` (${callCounts.PRESENT.toLocaleString("en-US")})`}
            </span>
            <span className="inline-flex items-center gap-1.5">
              <span className="w-3 h-3 rounded-sm" style={{ background: "var(--gv-call-partial)" }} />
              partial
              {callCounts && ` (${callCounts.PARTIAL.toLocaleString("en-US")})`}
            </span>
            <span className="inline-flex items-center gap-1.5">
              <span className="w-3 h-3 rounded-sm" style={{ background: "var(--gv-call-absent)" }} />
              absent
              {callCounts && ` (${callCounts.ABSENT.toLocaleString("en-US")})`}
            </span>
            <span className="text-zinc-400 dark:text-zinc-500">
              in {selectedQuery.query_name}
            </span>
            {!variantsOn ? (
              <span className="text-zinc-400 dark:text-zinc-500">
                variant markers off
              </span>
            ) : showMarkers ? (
              <>
                <KindToggle
                  on={variantKinds.snp}
                  color={SNP_COLOR}
                  label="SNP"
                  onToggle={() => setVariantKinds((v) => ({ ...v, snp: !v.snp }))}
                />
                <KindToggle
                  on={variantKinds.ins}
                  color={INS_COLOR}
                  label="insertion"
                  onToggle={() => setVariantKinds((v) => ({ ...v, ins: !v.ins }))}
                />
                <KindToggle
                  on={variantKinds.del}
                  color={DEL_COLOR}
                  label="deletion"
                  wide
                  onToggle={() => setVariantKinds((v) => ({ ...v, del: !v.del }))}
                />
              </>
            ) : (
              <span className="text-zinc-400 dark:text-zinc-500">
                zoom in below {(VARIANT_SPAN / 1000).toFixed(0)} kb to see SNP/indel markers
              </span>
            )}
          </>
        ) : (
          <>
            <span className="inline-flex items-center gap-1.5">
              <span className="w-3 h-3 rounded-sm" style={{ background: "var(--gv-gene-cds)" }} />
              protein coding
            </span>
            <span className="inline-flex items-center gap-1.5">
              <span className="w-3 h-3 rounded-sm" style={{ background: "var(--gv-gene-trna)" }} />
              tRNA
            </span>
            <span className="inline-flex items-center gap-1.5">
              <span className="w-3 h-3 rounded-sm" style={{ background: "var(--gv-gene-rrna)" }} />
              rRNA
            </span>
            <span className="inline-flex items-center gap-1.5">
              <span className="w-3 h-3 rounded-sm" style={{ background: "var(--gv-gene-pseudo)" }} />
              pseudo / other
            </span>
          </>
        )}
        <span>
          Pinch or scroll-wheel to zoom, two-finger scroll or drag to pan. Click a
          gene for details.
        </span>
      </div>

      {popup && (
        <GenePopup
          popup={popup}
          onClose={() => setPopup(null)}
          onOpenGene={(locus) => {
            setPopup(null);
            onOpenGene(locus);
          }}
          callInfo={
            selectedQuery && popup.gi !== undefined
              ? {
                  name: selectedQuery.query_name,
                  call: selectedQuery.calls?.[popup.gi] ?? "ABSENT",
                  covPct: selectedQuery.cov_pcts?.[popup.gi] ?? 0,
                  identity: selectedQuery.identities?.[popup.gi] ?? 0,
                }
              : null
          }
        />
      )}
    </div>
  );
}

/** What a variant is, on hover: the substitution, or the length of the indel. */
function VariantTooltip({
  variant,
  x,
  y,
  queryName,
  svgWidth,
}: {
  variant: VariantInfo;
  x: number;
  y: number;
  queryName: string;
  svgWidth: number;
}) {
  const color =
    variant.kind === "snp"
      ? SNP_COLOR
      : variant.kind === "del"
        ? DEL_COLOR
        : INS_COLOR;
  const title =
    variant.kind === "snp"
      ? "Mismatch"
      : variant.kind === "del"
        ? "Deletion"
        : "Insertion";
  return (
    <div
      className="pointer-events-none absolute z-30 w-60 rounded-lg border border-zinc-200 bg-white p-2.5 text-xs shadow-xl dark:border-zinc-800 dark:bg-zinc-900"
      style={{ left: Math.max(4, Math.min(x - 120, svgWidth - 244)), top: y }}
    >
      <p className="font-semibold text-[13px]" style={{ color }}>
        {title}
      </p>
      <p className="mt-1 font-mono text-[11px] text-zinc-500 dark:text-zinc-400">
        {variant.kind === "ins"
          ? `after reference ${variant.pos.toLocaleString("en-US")}`
          : `reference ${variant.pos.toLocaleString("en-US")}`}
        {variant.kind === "del" &&
          ` – ${(variant.pos + variant.len - 1).toLocaleString("en-US")}`}
      </p>
      {variant.kind === "snp" && (
        <p className="mt-1.5 font-mono text-[13px]">
          {String.fromCharCode(variant.r)} <span className="text-zinc-400">&rarr;</span>{" "}
          <span style={{ color }}>{String.fromCharCode(variant.q)}</span>
          <span className="ml-1.5 text-[11px] text-zinc-500 dark:text-zinc-400">
            (reference &rarr; query)
          </span>
        </p>
      )}
      {variant.kind === "del" && (
        <p className="mt-1.5 text-[12px]">
          <span className="font-mono font-semibold" style={{ color }}>
            {variant.len.toLocaleString("en-US")} bp
          </span>{" "}
          missing from the query
        </p>
      )}
      {variant.kind === "ins" && (
        <p className="mt-1.5 text-[12px] break-all">
          <span className="font-mono font-semibold" style={{ color }}>
            {variant.seq.length.toLocaleString("en-US")} bp
          </span>{" "}
          inserted:{" "}
          <span className="font-mono">
            {variant.seq.slice(0, 60)}
            {variant.seq.length > 60 ? "…" : ""}
          </span>
        </p>
      )}
      <p className="mt-1.5 text-[11px] text-zinc-400 dark:text-zinc-500">{queryName}</p>
    </div>
  );
}

/** A legend entry that doubles as the on/off switch for that kind of variant. */
function KindToggle({
  on,
  color,
  label,
  wide,
  onToggle,
}: {
  on: boolean;
  color: string;
  label: string;
  wide?: boolean;
  onToggle: () => void;
}) {
  return (
    <button
      type="button"
      onClick={onToggle}
      aria-pressed={on}
      title={`${on ? "Hide" : "Show"} ${label} markers`}
      className={`inline-flex items-center gap-1.5 -mx-1 px-1 rounded cursor-pointer hover:bg-zinc-100 dark:hover:bg-zinc-800 ${
        on ? "" : "opacity-40"
      }`}
    >
      <span
        className={`${wide ? "w-3" : "w-1"} h-3 rounded-sm`}
        style={{ background: on ? color : "transparent", outline: `1px solid ${color}` }}
      />
      <span className={on ? "" : "line-through"}>{label}</span>
    </button>
  );
}

/** A variant as the map holds it, before it is attributed to a query. */
type VariantInfo =
  | { kind: "snp"; pos: number; r: number; q: number }
  | { kind: "del"; pos: number; len: number }
  | { kind: "ins"; pos: number; seq: string };

type VariantPopupInfo = VariantInfo & { queryName: string };

/** A gene bead: a rectangle with a strand arrow tip. */
function GeneShape({
  x,
  w,
  y,
  h,
  strand,
  fill,
}: {
  x: number;
  w: number;
  y: number;
  h: number;
  strand: number;
  fill: string;
}) {
  const tip = Math.min(9, Math.max(3, h / 4));
  if (w < 2 * tip) {
    return (
      <rect x={x} y={y} width={w} height={h} rx={1.5} style={{ fill }} />
    );
  }
  const x2 = x + w;
  const points =
    strand >= 0
      ? `${x},${y} ${x2 - tip},${y} ${x2},${y + h / 2} ${x2 - tip},${y + h} ${x},${y + h}`
      : `${x + tip},${y} ${x2},${y} ${x2},${y + h} ${x + tip},${y + h} ${x},${y + h / 2}`;
  return <polygon points={points} style={{ fill }} />;
}

function truncateLabel(label: string, maxPx: number): string {
  const charW = 6.2; // 10px monospace approx
  const maxChars = Math.floor(maxPx / charW);
  if (maxChars < 2) return "";
  if (label.length <= maxChars) return label;
  return `${label.slice(0, Math.max(1, maxChars - 1))}\u2026`;
}

function GenePopup({
  popup,
  onClose,
  onOpenGene,
  callInfo,
}: {
  popup: { x: number; y: number; gene?: WgaGene; variant?: VariantPopupInfo };
  onClose: () => void;
  onOpenGene: (locus: string) => void;
  callInfo?: { name: string; call: Call; covPct: number; identity: number } | null;
}) {
  const ref = useRef<HTMLDivElement>(null);
  useEffect(() => {
    function onDocClick(e: MouseEvent) {
      if (ref.current && !ref.current.contains(e.target as Node)) onClose();
    }
    document.addEventListener("mousedown", onDocClick);
    return () => document.removeEventListener("mousedown", onDocClick);
  }, [onClose]);
  return (
    <div
      ref={ref}
      className="absolute z-30 bg-white border border-zinc-200 rounded-lg shadow-lg p-3 w-72 dark:bg-zinc-900 dark:border-zinc-800"
      style={{ left: Math.min(popup.x, 800), top: popup.y + 16 }}
    >
      {popup.gene && (
        <>
          <p className="font-semibold text-[15px]">{popup.gene.locus_tag}</p>
          <p className="text-xs text-zinc-500 mt-1 dark:text-zinc-400">
            {popup.gene.symbol ? `${popup.gene.symbol} - ` : ""}
            {popup.gene.biotype} on {popup.gene.seqid}
          </p>
          <p className="text-xs text-zinc-500 mt-1 font-mono dark:text-zinc-400">
            {popup.gene.start.toLocaleString("en-US")} -{" "}
            {popup.gene.end.toLocaleString("en-US")} (
            {popup.gene.strand > 0 ? "+" : "-"} strand)
          </p>
          {popup.gene.product && (
            <p className="text-xs text-zinc-600 mt-2 leading-snug dark:text-zinc-300">
              {popup.gene.product}
            </p>
          )}
          <FunctionBadges product={popup.gene.product} />
          {callInfo && (
            <div className="mt-2 flex items-center gap-2 text-xs">
              <CallBadge call={callInfo.call} />
              <span className="font-mono text-zinc-500 dark:text-zinc-400">
                {callInfo.covPct.toFixed(1)}% coverage,{" "}
                {callInfo.identity.toFixed(1)}% identity
              </span>
            </div>
          )}
          <button
            className="mt-3 w-full h-10 rounded-md bg-zinc-900 text-white text-sm hover:bg-zinc-700 dark:bg-zinc-100 dark:text-zinc-900 dark:hover:bg-zinc-300"
            onClick={() => onOpenGene(popup.gene!.locus_tag)}
          >
            Show alignment
          </button>
        </>
      )}
      {popup.variant && (
        <>
          <p className="font-semibold text-[15px]">
            {popup.variant.kind === "snp"
              ? "SNP"
              : popup.variant.kind === "del"
                ? "Deletion"
                : "Insertion"}
          </p>
          <p className="text-xs text-zinc-500 mt-1 dark:text-zinc-400">
            {popup.variant.queryName}
          </p>
          <p className="text-xs text-zinc-500 mt-1 font-mono dark:text-zinc-400">
            {popup.variant.kind === "ins"
              ? `after reference position ${popup.variant.pos.toLocaleString("en-US")}`
              : `reference position ${popup.variant.pos.toLocaleString("en-US")}`}
          </p>
          {popup.variant.kind === "snp" && (
            <p className="text-xs mt-1 font-mono dark:text-zinc-300">
              {String.fromCharCode(popup.variant.r)} &rarr;{" "}
              <span style={{ color: SNP_COLOR }}>
                {String.fromCharCode(popup.variant.q)}
              </span>
            </p>
          )}
          {popup.variant.kind === "del" && (
            <p className="text-xs mt-1 font-mono dark:text-zinc-300">
              <span style={{ color: DEL_COLOR }}>
                {popup.variant.len.toLocaleString("en-US")} bp
              </span>{" "}
              missing from the query
            </p>
          )}
          {popup.variant.kind === "ins" && (
            <p className="text-xs mt-1 font-mono break-all dark:text-zinc-300">
              <span style={{ color: INS_COLOR }}>
                {popup.variant.seq.length.toLocaleString("en-US")} bp
              </span>{" "}
              inserted: {popup.variant.seq.slice(0, 200)}
              {popup.variant.seq.length > 200 ? "\u2026" : ""}
            </p>
          )}
        </>
      )}
    </div>
  );
}

/** Hover tooltip for a gene: name, position, function annotation and,
 * when a query is selected, its presence call in that query. */
function GeneTooltip({
  gene,
  x,
  y,
  queryName,
  call,
  covPct,
  identity,
  svgWidth,
}: {
  gene: WgaGene;
  x: number;
  y: number;
  queryName: string | null;
  call?: Call;
  covPct?: number;
  identity?: number;
  svgWidth: number;
}) {
  return (
    <div
      className="pointer-events-none absolute z-20 w-64 rounded-lg border border-zinc-200 bg-white p-2.5 text-xs shadow-xl dark:border-zinc-800 dark:bg-zinc-900"
      style={{
        left: Math.max(4, Math.min(x - 128, svgWidth - 264)),
        top: y + 6,
      }}
    >
      <p className="font-semibold text-[13px]">
        {gene.locus_tag}
        {gene.symbol && (
          <span className="font-normal text-zinc-500 dark:text-zinc-400">
            {" "}
            ({gene.symbol})
          </span>
        )}
      </p>
      <p className="mt-0.5 font-mono text-[11px] text-zinc-500 dark:text-zinc-400">
        {gene.seqid}:{gene.start.toLocaleString("en-US")}-
        {gene.end.toLocaleString("en-US")} ({gene.strand > 0 ? "+" : "-"}),
        {gene.biotype}
      </p>
      {gene.product && (
        <p className="mt-1 leading-snug text-zinc-700 dark:text-zinc-300">
          {gene.product}
        </p>
      )}
      <FunctionBadges product={gene.product} />
      {queryName && call && (
        <p className="mt-1.5 flex items-center gap-1.5">
          <CallBadge call={call} />
          <span className="font-mono text-[11px] text-zinc-500 dark:text-zinc-400">
            {covPct?.toFixed(1)}% cov, {identity?.toFixed(1)}% id
          </span>
        </p>
      )}
    </div>
  );
}

function CallBadge({ call }: { call: Call }) {
  const { label, cls } =
    call === "PRESENT"
      ? {
          label: "present",
          cls: "bg-emerald-100 text-emerald-800 dark:bg-emerald-900/40 dark:text-emerald-300",
        }
      : call === "PARTIAL"
        ? {
            label: "partial",
            cls: "bg-amber-100 text-amber-800 dark:bg-amber-900/40 dark:text-amber-300",
          }
        : {
            label: "absent",
            cls: "bg-zinc-100 text-zinc-600 dark:bg-zinc-800 dark:text-zinc-400",
          };
  return (
    <span className={`inline-flex h-5 shrink-0 items-center rounded px-1.5 text-[11px] font-medium ${cls}`}>
      {label}
    </span>
  );
}

const VIRULENCE_KEYWORDS = [
  "virulen",
  "toxin",
  "hemolys",
  "haemolys",
  "leukocidin",
  "cytolys",
  "adhesin",
  "adhesi",
  "invasin",
  "invasion",
  "hemagglutin",
  "haemagglutin",
  "siderophore",
  "aerobactin",
  "enterobactin",
  "enterochelin",
  "iga protease",
  "autotransporter",
  "rtx ",
  "secretion system",
  "immune evasion",
  "capsular polysaccharide",
];

const RESISTANCE_KEYWORDS = [
  "resistance",
  "beta-lactamase",
  "lactamase",
  "aminoglycoside",
  "tetracycline",
  "chloramphenicol",
  "macrolide",
  "lincosamide",
  "sulfonamide",
  "trimethoprim",
  "vancomycin",
  "fosfomycin",
  "rifampin",
  "rifampicin",
  "quinolone",
  "carbapenem",
  "multidrug",
  "efflux pump",
  "mdr ",
];

const MOBILE_ELEMENT_KEYWORDS = [
  "transposase",
  "integrase",
  "recombinase",
  "insertion sequence",
  "phage",
  "prophage",
  "plasmid",
  "integron",
  "resolvase",
  "relaxase",
  "mobilization",
];

const FUNCTION_BADGE_CLS: Record<string, string> = {
  "virulence-related": "bg-red-100 text-red-800 dark:bg-red-900/40 dark:text-red-300",
  "resistance-related":
    "bg-violet-100 text-violet-800 dark:bg-violet-900/40 dark:text-violet-300",
  "mobile element": "bg-sky-100 text-sky-800 dark:bg-sky-900/40 dark:text-sky-300",
};

/** Rough keyword based classification of a GFF product string. */
function functionCategories(product: string | undefined): string[] {
  if (!product) return [];
  const p = product.toLowerCase();
  const cats: string[] = [];
  if (VIRULENCE_KEYWORDS.some((k) => p.includes(k))) cats.push("virulence-related");
  if (RESISTANCE_KEYWORDS.some((k) => p.includes(k))) cats.push("resistance-related");
  if (MOBILE_ELEMENT_KEYWORDS.some((k) => p.includes(k))) cats.push("mobile element");
  return cats;
}

function FunctionBadges({ product }: { product: string | undefined }) {
  const cats = functionCategories(product);
  if (cats.length === 0) return null;
  return (
    <div className="mt-1.5 flex flex-wrap gap-1">
      {cats.map((c) => (
        <span
          key={c}
          className={`inline-flex h-5 items-center rounded px-1.5 text-[11px] font-medium ${FUNCTION_BADGE_CLS[c]}`}
        >
          {c}
        </span>
      ))}
    </div>
  );
}

function RangeInput({
  range,
  seqLength,
  onChange,
}: {
  range: { start: number; end: number };
  seqLength: number;
  onChange: (r: { start: number; end: number }) => void;
}) {
  const [start, setStart] = useState(String(Math.round(range.start)));
  const [end, setEnd] = useState(String(Math.round(range.end)));
  useEffect(() => {
    setStart(String(Math.round(range.start)));
    setEnd(String(Math.round(range.end)));
  }, [range.start, range.end]);
  return (
    <form
      className="flex items-center gap-1"
      onSubmit={(e) => {
        e.preventDefault();
        const s = Math.max(1, Number(start) || 1);
        const e2 = Math.min(seqLength, Number(end) || seqLength);
        if (e2 > s) onChange({ start: s, end: e2 });
      }}
    >
      <input
        className="h-11 w-24 px-2 rounded-lg border border-zinc-300 font-mono text-sm dark:border-zinc-700 dark:bg-zinc-900"
        value={start}
        onChange={(e) => setStart(e.target.value)}
        aria-label="Range start"
      />
      <span className="text-zinc-400 dark:text-zinc-600">-</span>
      <input
        className="h-11 w-24 px-2 rounded-lg border border-zinc-300 font-mono text-sm dark:border-zinc-700 dark:bg-zinc-900"
        value={end}
        onChange={(e) => setEnd(e.target.value)}
        aria-label="Range end"
      />
      <button
        type="submit"
        className="h-11 px-3 rounded-lg border border-zinc-300 bg-white hover:bg-zinc-100 text-sm dark:border-zinc-700 dark:bg-zinc-900 dark:hover:bg-zinc-800"
      >
        Go
      </button>
    </form>
  );
}

function niceStep(span: number): number {
  const target = span / 8;
  const pow = Math.pow(10, Math.floor(Math.log10(Math.max(1, target))));
  const candidates = [1, 2, 5, 10].map((m) => m * pow);
  return candidates.find((c) => c >= target) ?? 10 * pow;
}

function geneColor(g: WgaGene): string {
  if (g.biotype === "tRNA" || g.biotype === "tmRNA") return "var(--gv-gene-trna)";
  if (g.biotype === "rRNA") return "var(--gv-gene-rrna)";
  if (g.biotype.startsWith("pseudo")) return "var(--gv-gene-pseudo)";
  return "var(--gv-gene-cds)";
}

/** Gene fill by presence call in the selected query. */
function callColor(call: Call | undefined): string {
  if (call === "PRESENT") return "var(--gv-call-present)";
  if (call === "PARTIAL") return "var(--gv-call-partial)";
  return "var(--gv-call-absent)";
}
