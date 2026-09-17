import { useCallback, useEffect, useMemo, useRef, useState } from "react";
import { api } from "../api";
import type { AlignmentData, DelEvent, InsEvent, RefseqWindow, Run, SnpEvent } from "../types";
import { Spinner } from "./ui";
import { useWheelGestures } from "./useWheelGestures";

/**
 * Whole-genome alignment viewer in the style of Oligool's MSA viewer
 * (zoom logic and colors ported from it). Row 0 is the reference
 * genome; every row below is one query FASTA aligned against it.
 * Red = mismatch, blue = insertion, purple = deletion, gray = same as
 * the reference, blank = not aligned. Zooming in past ~100 bp switches
 * to base letters (fetched from the server for the visible window).
 * The reference contigs are concatenated into one coordinate space,
 * with separator lines and per-contig ruler numbers.
 */

/* ── layout constants (ported from Oligool's MSAViewer) ── */
const RIGHT_PADDING = 20;
const CONTIG_BAND_H = 16;
const RULER_HEIGHT = 24;
const ROW_HEIGHT = 18;
const MAX_VIEWER_HEIGHT = 560;
const BP_THRESHOLD = 100;
const HYSTERESIS = 15;
/** Auto-switch threshold between bars and letters mode. */
const REFSEQ_MAX_WINDOW = 8192;
/**
 * Deepest zoom: this many bases across the sequence area, so the view always
 * reaches base letters. A fixed fraction can't do it: 0.005 of a 5 Mb genome
 * is 25 kb, which kept letters mode mathematically unreachable.
 */
const MIN_VISIBLE_BASES = 40;
/**
 * Zoomed out, point variants are shaded by the share of bases in their pixel
 * column that actually differ (same contract the strain map uses), with a
 * floor so a lone mismatch is still visible.
 */
const MIN_DENSITY_INK = 0.15;

/* ── colors (Oligool palette; event colors shared with the StrainMap) ── */
const SNP_COLOR = "#dc2626";
const INS_COLOR = "#3b82f6";
const DEL_COLOR = "#9333ea";

const getPrettyStep = (minPixels: number, pixelsPerUnit: number) => {
  const minUnits = minPixels / pixelsPerUnit;
  const steps = [
    1, 2, 5, 10, 20, 50, 100, 250, 500, 1000, 2500, 5000, 10000, 25000, 50000,
    100000, 250000, 500000, 1000000, 2500000, 5000000,
  ];
  return steps.find((s) => s >= minUnits) || steps[steps.length - 1];
};

interface PreparedContig {
  seqid: string;
  /** 0-based start in the concatenated coordinate space. */
  start: number;
  len: number;
}

interface PreparedRow {
  queryId: number;
  name: string;
  /** Merged aligned intervals on the reference, concat coords, sorted. */
  spans: [number, number][];
  /** Raw alignment blocks (for the rev markers), concat coords. */
  blocks: { s: number; e: number; rev: boolean }[];
  /** All three event lists are sorted by pos for binary-search lookups
   * and visible-window culling (the draw loop never scans whole-genome
   * event lists: with divergent queries there are tens of thousands of
   * SNPs per row and a full scan per frame made panning laggy). */
  snpList: SnpEvent[];
  delList: DelEvent[];
  insList: InsEvent[];
}

function inSpans(spans: [number, number][], a: number): boolean {
  let lo = 0;
  let hi = spans.length - 1;
  while (lo <= hi) {
    const mid = (lo + hi) >>> 1;
    const [s, e] = spans[mid];
    if (a < s) hi = mid - 1;
    else if (a > e) lo = mid + 1;
    else return true;
  }
  return false;
}

/** First index whose pos >= target (list.length when all are before). */
function lowerBoundPos<T extends { pos: number }>(list: T[], target: number): number {
  let lo = 0;
  let hi = list.length;
  while (lo < hi) {
    const mid = (lo + hi) >>> 1;
    if (list[mid].pos < target) lo = mid + 1;
    else hi = mid;
  }
  return lo;
}

/** SNP at an anchor exactly ('snps' is a point list), or null. */
function snpAt<T extends { pos: number }>(list: T[], a: number): T | null {
  const i = lowerBoundPos(list, a);
  return i < list.length && list[i].pos === a ? list[i] : null;
}

/** Insertion anchored at an anchor exactly, or null. */
function insAt<T extends { pos: number }>(list: T[], a: number): T | null {
  const i = lowerBoundPos(list, a);
  return i < list.length && list[i].pos === a ? list[i] : null;
}

/** Deletion covering an anchor, or null. Deletions do not start twice at
 * the same position and nucmer blocks do not overlap, so the last one
 * starting at or before the anchor is the only candidate. */
function delCovering(list: DelEvent[], a: number): DelEvent | null {
  const i = lowerBoundPos(list, a);
  const d = i < list.length && list[i].pos === a ? list[i] : i > 0 ? list[i - 1] : null;
  return d && a >= d.pos && a < d.pos + d.len ? d : null;
}

interface Prepared {
  contigs: PreparedContig[];
  totalLen: number;
  rows: PreparedRow[];
  offsetOf: Map<string, number>;
}

/** Building the sorted per-row event lists takes real time on divergent
 * runs (hundreds of thousands of events), and the component unmounts
 * whenever the user toggles back to the StrainMap — a fresh memo per
 * visit made every toggle as slow as the first. The prepared form is
 * cached per fetched payload, which the promise cache keeps alive. */
const preparedCache = new WeakMap<AlignmentData, Prepared>();

function prepareData(data: AlignmentData): Prepared {
  let prepared = preparedCache.get(data);
  if (prepared) return prepared;
  const contigs: PreparedContig[] = [];
  let off = 0;
  for (const [seqid, len] of data.reference) {
    contigs.push({ seqid, start: off, len });
    off += len;
  }
  const totalLen = off;
  const offsetOf = new Map(contigs.map((c) => [c.seqid, c.start]));

  const rows: PreparedRow[] = data.queries.map((q) => {
    const raw = q.blocks
      .map((b) => ({
        s: (offsetOf.get(b.ref_seqid) ?? 0) + b.ref_start - 1,
        e: (offsetOf.get(b.ref_seqid) ?? 0) + b.ref_end - 1,
        rev: b.qry_rev,
      }))
      .sort((a, b) => a.s - b.s);
    // merged spans
    const spans: [number, number][] = [];
    for (const b of raw) {
      const last = spans[spans.length - 1];
      if (last && b.s <= last[1] + 1) last[1] = Math.max(last[1], b.e);
      else spans.push([b.s, b.e]);
    }
    // events, converted into the concatenated space, in concat order
    const snpList: SnpEvent[] = [];
    const delList: DelEvent[] = [];
    const insList: InsEvent[] = [];
    for (const c of contigs) {
      const ev = q.events[c.seqid];
      if (!ev) continue;
      for (const s of ev.snps) {
        snpList.push({ ...s, pos: c.start + s.pos - 1 });
      }
      for (const d of ev.dels) {
        delList.push({ ...d, pos: c.start + d.pos - 1 });
      }
      for (const i of ev.ins) {
        // insertion sits after anchor (c.start + i.pos - 1); the tick
        // goes on the boundary anchor that follows it
        insList.push({ ...i, pos: c.start + i.pos });
      }
    }
    snpList.sort((a, b) => a.pos - b.pos);
    delList.sort((a, b) => a.pos - b.pos);
    insList.sort((a, b) => a.pos - b.pos);
    return {
      queryId: q.query_id,
      name: q.query_name,
      spans,
      blocks: raw,
      snpList,
      delList,
      insList,
    };
  });

  prepared = { contigs, totalLen, rows, offsetOf };
  preparedCache.set(data, prepared);
  return prepared;
}

export function AlignmentView({
  run,
  initialRange,
  onLocusChange,
}: {
  run: Run;
  /** {seqid, start, end}: where to center the view on first load. */
  initialRange?: { seqid: string; start: number; end: number } | null;
  onLocusChange?: (r: { seqid: string; start: number; end: number }) => void;
}) {
  const [data, setData] = useState<AlignmentData | null>(null);
  const [error, setError] = useState<string | null>(null);

  const canvasRef = useRef<HTMLCanvasElement>(null);
  const hoverOverlayRef = useRef<HTMLCanvasElement>(null);
  const scrollRef = useRef<HTMLDivElement>(null);
  const [availableWidth, setAvailableWidth] = useState(900);
  /** Virtual horizontal offset in css px: no DOM node is ever this wide. */
  const [scrollLeft, setScrollLeft] = useState(0);
  const [viewFraction, setViewFraction] = useState(1);
  const [viewMode, setViewMode] = useState<"bars" | "letters">("bars");
  const [labelWidth, setLabelWidth] = useState(140);
  const [isResizingLabel, setIsResizingLabel] = useState(false);
  const hoverColRef = useRef<number | null>(null);
  const hoverRafRef = useRef(0);
  const redrawRef = useRef<() => void>(() => {});
  const scrollPressRef = useRef<ReturnType<typeof setInterval> | null>(null);
  const [isBarDragging, setIsBarDragging] = useState(false);
  const [zoomRect, setZoomRect] = useState<{ a0: number; a1: number } | null>(null);
  const [hoverInfo, setHoverInfo] = useState<{
    x: number;
    y: number;
    flipX?: boolean;
    title: string;
    lines: string[];
  } | null>(null);
  const [labelHover, setLabelHover] = useState<{ text: string; x: number; y: number } | null>(null);

  // live refs for use inside setInterval
  const viewFractionRef = useRef(viewFraction);
  const seqAreaWRef = useRef(0);
  const totalVirtualWRef = useRef(0);

  const [isDark, setIsDark] = useState(() =>
    document.documentElement.classList.contains("dark"),
  );
  useEffect(() => {
    const obs = new MutationObserver(() =>
      setIsDark(document.documentElement.classList.contains("dark")),
    );
    obs.observe(document.documentElement, {
      attributes: true,
      attributeFilter: ["class"],
    });
    return () => obs.disconnect();
  }, []);

  useEffect(() => {
    let cancelled = false;
    api
      .alignment(run.id)
      .then((d) => !cancelled && setData(d))
      .catch((e) => !cancelled && setError((e as Error).message));
    return () => {
      cancelled = true;
    };
  }, [run.id]);

  /* ── prepare the concatenated coordinate space and per-query indexes ── */
  const prepared = useMemo(() => (data ? prepareData(data) : null), [data]);

  /* ── auto label width from the row names ── */
  useEffect(() => {
    if (!prepared) return;
    const names = ["Reference", ...prepared.rows.map((r) => r.name)];
    const canvas = document.createElement("canvas");
    const ctx = canvas.getContext("2d");
    if (!ctx) return;
    ctx.font = "bold 10px ui-sans-serif, system-ui, sans-serif";
    let maxW = 0;
    for (const n of names) maxW = Math.max(maxW, ctx.measureText(n).width);
    setLabelWidth(Math.min(280, Math.max(80, Math.ceil(maxW + 16))));
  }, [prepared]);

  /* ── sizing (ported) ── */
  const seqAreaW = Math.max(1, availableWidth - labelWidth - RIGHT_PADDING);
  const totalVirtualW = seqAreaW / viewFraction;
  const anchorLen = prepared?.totalLen ?? 0;
  /**
   * The horizontal offset is fully virtual: totalVirtualW only ever lives in
   * maths, never in the DOM. Reaching base level on a 5 Mb genome needs a
   * ~175-megapixel strip, several times the browsers' maximum layout width —
   * a virtual offset has no such ceiling, so zoom depth is limited only by
   * MIN_VISIBLE_BASES.
   */
  const minVF = anchorLen > 0 ? Math.min(1, MIN_VISIBLE_BASES / anchorLen) : 1;
  const cellW = anchorLen > 0 ? totalVirtualW / anchorLen : 1;
  const visibleBases = anchorLen * viewFraction;
  const headerH = CONTIG_BAND_H + RULER_HEIGHT;
  const totalH = headerH + (prepared ? prepared.rows.length + 1 : 0) * ROW_HEIGHT + 4;

  viewFractionRef.current = viewFraction;
  seqAreaWRef.current = seqAreaW;
  totalVirtualWRef.current = totalVirtualW;

  /* ── auto-switch bars/letters with hysteresis (ported) ── */
  useEffect(() => {
    // Before the data arrives anchorLen is 0, so visibleBases is 0 and
    // read as "zoomed past base level" — the viewer flipped into
    // letters mode while still showing nothing, and once the data
    // landed the first draw letter-rendered the WHOLE genome (tens of
    // millions of canvas ops, a ~half-minute frozen page). Only ever
    // switch on real geometry.
    if (!prepared) return;
    if (viewMode === "bars" && visibleBases < BP_THRESHOLD - HYSTERESIS) {
      setViewMode("letters");
    } else if (viewMode === "letters" && visibleBases > BP_THRESHOLD + HYSTERESIS) {
      setViewMode("bars");
    }
  }, [visibleBases, viewMode, prepared]);

  /* ── container resize tracking ── */
  // The observer lives in the setScrollEl callback ref, not in a mount
  // effect: until the data arrives the component renders only the spinner,
  // so a [] effect found no scroll container to observe and the map stayed
  // at the 900px fallback width forever. ResizeObserver also fires once on
  // observe, so the very first layout already gets the real width.
  const resizeObsRef = useRef<ResizeObserver | null>(null);
  useEffect(() => () => resizeObsRef.current?.disconnect(), []);

  /* ── zoom helpers ── */
  const applyZoom = (factor: number, anchorPx: number = seqAreaW / 2) => {
    const newVF = Math.max(minVF, Math.min(1, viewFraction * factor));
    const currentTotalVirtualW = seqAreaW / viewFraction;
    const anchorFracGlobal = (scrollLeft + anchorPx) / currentTotalVirtualW;
    const newTotalVirtualW = seqAreaW / newVF;
    let newSL = anchorFracGlobal * newTotalVirtualW - anchorPx;
    newSL = Math.max(0, Math.min(newTotalVirtualW - seqAreaW, newSL));
    setViewFraction(newVF);
    setScrollLeft(newSL);
  };
  const zoomIn = () => applyZoom(0.75, 0);
  const zoomOut = () => applyZoom(1.33, 0);

  /* ── first load: jump to the initial range, mirroring the zoom level ── */
  const initialDone = useRef(false);
  useEffect(() => {
    if (!prepared || initialDone.current || availableWidth <= 0) return;
    initialDone.current = true;
    if (!initialRange) return;
    const contig = prepared.contigs.find((c) => c.seqid === initialRange.seqid);
    if (!contig) return;
    const a0 = contig.start + Math.max(1, initialRange.start) - 1;
    const a1 = contig.start + initialRange.end - 1;
    const span = Math.max(40, a1 - a0 + 1);
    const newVF = Math.max(minVF, Math.min(1, span / prepared.totalLen));
    const centerFrac = (a0 + a1) / (2 * prepared.totalLen);
    const newTotalW = seqAreaW / newVF;
    let newSL = centerFrac * newTotalW - seqAreaW / 2;
    newSL = Math.max(0, Math.min(newTotalW - seqAreaW, newSL));
    setViewFraction(newVF);
    setScrollLeft(newSL);
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [prepared, availableWidth, initialRange]);

  /* ── reference bases for letters mode: cached window fetch ── */
  const refWindows = useRef(new Map<string, RefseqWindow[]>());
  const [refSeqVersion, setRefSeqVersion] = useState(0);

  const firstCol = Math.max(0, Math.floor(scrollLeft / cellW));
  const lastCol = Math.min(anchorLen - 1, Math.ceil((scrollLeft + seqAreaW) / cellW));

  useEffect(() => {
    if (!prepared || viewMode !== "letters") return;
    const t = setTimeout(() => {
      for (const c of prepared.contigs) {
        if (c.start > lastCol || c.start + c.len - 1 < firstCol) continue;
        const localStart = Math.max(1, firstCol - c.start + 1);
        const localEnd = Math.min(c.len, lastCol - c.start + 1);
        if (localEnd < localStart) continue;
        // align window starts to a 256 bp grid so panning reuses windows
        const ws = Math.max(1, Math.floor((localStart - 25) / 256) * 256 + 1);
        const we = Math.min(c.len, ws + REFSEQ_MAX_WINDOW - 1, localEnd + 25);
        const have = refWindows.current.get(c.seqid)?.some(
          (w) => w.start <= localStart && w.end >= localEnd,
        );
        if (have) continue;
        api
          .refseq(run.id, c.seqid, ws, we)
          .then((w) => {
            const list = refWindows.current.get(c.seqid) ?? [];
            list.push(w);
            refWindows.current.set(c.seqid, list);
            setRefSeqVersion((v) => v + 1);
          })
          .catch(() => {});
      }
    }, 150);
    return () => clearTimeout(t);
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [run.id, prepared, viewMode, firstCol, lastCol]);

  /** Reference base at a concatenated anchor, or null when not fetched. */
  const refBaseAt = useCallback(
    (a: number): string | null => {
      if (!prepared) return null;
      for (const c of prepared.contigs) {
        if (a >= c.start && a < c.start + c.len) {
          const local = a - c.start + 1;
          for (const w of refWindows.current.get(c.seqid) ?? []) {
            if (local >= w.start && local <= w.end) {
              return w.seq[local - w.start] ?? null;
            }
          }
          return null;
        }
      }
      return null;
    },
    [prepared],
  );

  /* ── publish the visible locus for mode switching (debounced) ── */
  useEffect(() => {
    if (!prepared || !onLocusChange || anchorLen === 0) return;
    const t = setTimeout(() => {
      const a0 = Math.max(0, firstCol);
      const a1 = Math.min(anchorLen - 1, Math.max(lastCol - 1, a0));
      const center = Math.floor((a0 + a1) / 2);
      const c = prepared.contigs.find(
        (x) => center >= x.start && center < x.start + x.len,
      );
      if (!c) return;
      const s = Math.max(c.start + 1, a0 - c.start + 1);
      const e = Math.min(c.len, a1 - c.start + 1);
      if (e < s) return;
      onLocusChange({ seqid: c.seqid, start: s, end: e });
    }, 600);
    return () => clearTimeout(t);
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [firstCol, lastCol, anchorLen, prepared]);

  /* ── pan helpers (ported) ── */
  const panBy = (dir: -1 | 1) => {
    const vf = viewFractionRef.current;
    const saw = seqAreaWRef.current;
    const tvw = totalVirtualWRef.current;
    const step = Math.max(40, saw * vf * 0.25);
    setScrollLeft((v) => Math.max(0, Math.min(tvw - saw, v + dir * step)));
  };
  const startPan = (dir: -1 | 1) => {
    panBy(dir);
    scrollPressRef.current = setInterval(() => panBy(dir), 120);
  };
  const stopPan = () => {
    if (scrollPressRef.current) {
      clearInterval(scrollPressRef.current);
      scrollPressRef.current = null;
    }
  };

  /* ── wheel = zoom, two-finger scroll / wheel = pan ── */
  // Attached natively rather than through React's onWheel: React registers wheel
  // as a passive listener, so preventDefault() there never runs and a macOS
  // trackpad pinch falls through to the browser's own page zoom.
  const setWheelEl = useWheelGestures<HTMLDivElement>((g) => {
    if (g.kind === "pan") {
      const saw = seqAreaWRef.current;
      const tvw = totalVirtualWRef.current;
      setScrollLeft((v) => Math.max(0, Math.min(tvw - saw, v + (g.dx || g.dy))));
      return;
    }
    const rect = scrollRef.current?.getBoundingClientRect();
    if (!rect) return;
    const offsetX = Math.max(
      0,
      Math.min(seqAreaW, g.clientX - rect.left - labelWidth),
    );
    const currentTotalVirtualW = seqAreaW / viewFraction;
    const mouseFracGlobal = (scrollLeft + offsetX) / currentTotalVirtualW;
    const newVF = Math.max(minVF, Math.min(1, viewFraction * g.factor));
    const newTotalVirtualW = seqAreaW / newVF;
    let newSL = mouseFracGlobal * newTotalVirtualW - offsetX;
    newSL = Math.max(0, Math.min(newTotalVirtualW - seqAreaW, newSL));
    setViewFraction(newVF);
    setScrollLeft(newSL);
  });

  const setScrollEl = useCallback(
    (n: HTMLDivElement | null) => {
      scrollRef.current = n;
      setWheelEl(n);
      resizeObsRef.current?.disconnect();
      resizeObsRef.current = null;
      if (n) {
        const obs = new ResizeObserver((entries) => {
          for (const entry of entries) setAvailableWidth(entry.contentRect.width);
        });
        obs.observe(n);
        resizeObsRef.current = obs;
      }
    },
    [setWheelEl],
  );

  /* ── main canvas drawing ── */
  const startFrac = totalVirtualW > 0 ? scrollLeft / totalVirtualW : 0;

  const draw = useCallback(() => {
    const cvs = canvasRef.current;
    if (!cvs || !prepared) return;
    const ctx = cvs.getContext("2d");
    if (!ctx) return;

    const canvasH = Math.min(totalH, MAX_VIEWER_HEIGHT);
    const dpr = window.devicePixelRatio || 1;
    cvs.width = availableWidth * dpr;
    cvs.height = canvasH * dpr;
    cvs.style.width = `${availableWidth}px`;
    cvs.style.height = `${canvasH}px`;
    ctx.scale(dpr, dpr);

    const nRows = prepared.rows.length + 1;
    ctx.fillStyle = isDark ? "#0f172a" : "#ffffff";
    ctx.fillRect(0, 0, availableWidth, canvasH);

    const fCol = Math.max(0, Math.floor(scrollLeft / cellW));
    const lCol = Math.min(anchorLen - 1, Math.ceil((scrollLeft + seqAreaW) / cellW));
    const colStep = cellW > 0 ? Math.max(1, Math.floor(1 / cellW)) : 1;
    const barW = Math.max(1, Math.ceil(cellW * colStep));

    /* ── rows (clipped to the sequence area) ── */
    ctx.save();
    ctx.beginPath();
    ctx.rect(labelWidth, 0, availableWidth - labelWidth, canvasH);
    ctx.clip();

    for (let row = 0; row < nRows; row++) {
      const y = headerH + row * ROW_HEIGHT;
      const isRef = row === 0;
      const r = isRef ? null : prepared.rows[row - 1];

      // aligned span bars
      if (isRef) {
        ctx.fillStyle = isDark ? "#1e3a8a" : "#bfdbfe";
        ctx.fillRect(labelWidth, y + 3, seqAreaW, ROW_HEIGHT - 6);
      } else if (r) {
        for (const [s, e] of r.spans) {
          if (e < fCol || s > lCol + 1) continue;
          const barX1 = Math.floor(Math.max(labelWidth, labelWidth + s * cellW - scrollLeft));
          const barX2 = Math.ceil(
            Math.min(labelWidth + seqAreaW, labelWidth + (e + 1) * cellW - scrollLeft),
          );
          if (barX2 > barX1) {
            // the match fill is a true gray, not a slate: slate read as
            // "blue" next to the red/blue event colors
            ctx.fillStyle = isDark ? "#3f3f46" : "#e4e4e7";
            ctx.fillRect(barX1, y + 3, barX2 - barX1, ROW_HEIGHT - 6);
          }
        }
        // rev markers on wide blocks
        for (const b of r.blocks) {
          if (!b.rev) continue;
          const bx1 = labelWidth + b.s * cellW - scrollLeft;
          const bw = (b.e - b.s + 1) * cellW;
          if (bw < 24 || bx1 > labelWidth + seqAreaW || bx1 + bw < labelWidth) continue;
          ctx.fillStyle = isDark ? "#94a3b8" : "#64748b";
          ctx.font = "8px ui-monospace, SFMono-Regular, monospace";
          ctx.textAlign = "center";
          ctx.textBaseline = "middle";
          ctx.fillText("rev", bx1 + bw / 2, y + ROW_HEIGHT / 2);
        }
      }

      if (viewMode === "letters") {
        // step at most one column per pixel: if a letters draw ever
        // runs against a zoomed-out window (a pending mode/zoom state
        // change), it must not iterate every anchor of the genome
        for (let a = fCol; a <= lCol; a += colStep) {
          const x = labelWidth + a * cellW - scrollLeft;
          if (x + cellW < labelWidth || x > labelWidth + seqAreaW) continue;
          let bg = isDark ? "#27272a" : "#f4f4f5";
          let fg = isDark ? "#cbd5e1" : "#374151";
          let ch: string | null;

          if (isRef) {
            ch = refBaseAt(a);
            if (ch === null) ch = "";
          } else if (r) {
            if (!inSpans(r.spans, a)) {
              ch = "-";
              fg = isDark ? "#475569" : "#9ca3af";
            } else if (delCovering(r.delList, a)) {
              ch = "-";
              bg = isDark ? "#3b0764" : "#f3e8ff";
              fg = isDark ? "#d8b4fe" : "#7e22ce";
            } else {
              const snp = snpAt(r.snpList, a);
              const refB = refBaseAt(a);
              if (snp) {
                ch = String.fromCharCode(snp.q);
                bg = isDark ? "#7f1d1d" : "#fee2e2";
                fg = isDark ? "#fecaca" : "#b91c1c";
              } else {
                ch = refB ?? "";
              }
            }
          } else {
            ch = "";
          }

          ctx.fillStyle = bg;
          ctx.fillRect(x, y, cellW + 0.5, ROW_HEIGHT);
          if (ch) {
            ctx.fillStyle = fg;
            const fs = Math.min(13, Math.max(8, cellW * 0.8));
            ctx.font = `${fs}px ui-monospace, SFMono-Regular, monospace`;
            ctx.textAlign = "center";
            ctx.textBaseline = "middle";
            ctx.fillText(ch, x + cellW / 2, y + ROW_HEIGHT / 2);
          }
        }
        // insertion markers with extent lines (ported), windowed by anchor
        if (r) {
          const lo = lowerBoundPos(r.insList, fCol);
          const hi = lowerBoundPos(r.insList, lCol + 2);
          for (let i = lo; i < hi; i++) {
            const ins = r.insList[i];
            const bx = labelWidth + ins.pos * cellW - scrollLeft;
            if (bx < labelWidth || bx > labelWidth + seqAreaW) continue;
            const n = ins.seq.length;
            ctx.fillStyle = INS_COLOR;
            if (n > 1) {
              const leftX = Math.max(labelWidth, bx);
              const rightX = Math.min(labelWidth + seqAreaW, bx + Math.min(n, 60) * cellW);
              if (rightX > leftX) {
                ctx.fillRect(leftX, y + 1.5, rightX - leftX, 1.5);
                ctx.fillRect(leftX, y + ROW_HEIGHT - 3, rightX - leftX, 1.5);
              }
            }
            ctx.fillRect(Math.floor(bx) - 1, y + 1, 2, ROW_HEIGHT - 2);
          }
        }
      } else {
        // bars mode: variant markers on top of the span bars, windowed
        // to the visible columns (point events bisect; deletions stay a
        // full scan because their interval can start far before fCol
        // and still reach into the window, and they are few)
        if (r) {
          const snpLo = lowerBoundPos(r.snpList, fCol);
          const snpHi = lowerBoundPos(r.snpList, lCol + 1);
          // Zoomed out, many SNPs land on the same pixel: draw one column per
          // pixel, shaded by the share of the column's bases that differ, and
          // binary-search straight to the next affected pixel, so the loop is
          // bounded by the sequence area width rather than by the visible
          // event count (a divergent query shows >100k SNPs in a whole-genome
          // overview). Full-strength marks only once bases get pixel-wide:
          // drawn solid, a 3% divergence read as a solid red row.
          const mergePixels = cellW < 1;
          if (!mergePixels) {
            for (let i = snpLo; i < snpHi; i++) {
              const x = Math.floor(labelWidth + r.snpList[i].pos * cellW - scrollLeft);
              const w = Math.min(barW, labelWidth + seqAreaW - x);
              if (w <= 0) break; // sorted: the rest are past the right edge
              ctx.globalAlpha = 1;
              ctx.fillStyle = SNP_COLOR;
              ctx.fillRect(x, y + 2, w, ROW_HEIGHT - 4);
            }
          } else {
            let i = snpLo;
            while (i < snpHi) {
              const x = Math.floor(labelWidth + r.snpList[i].pos * cellW - scrollLeft);
              if (x > labelWidth + seqAreaW) break;
              const nextPos = Math.ceil((x + 1 - labelWidth + scrollLeft) / cellW);
              const next = lowerBoundPos(r.snpList, nextPos);
              if (x >= labelWidth) {
                ctx.globalAlpha = Math.min(
                  1,
                  Math.max(MIN_DENSITY_INK, (next - i) * cellW),
                );
                ctx.fillStyle = SNP_COLOR;
                ctx.fillRect(x, y + 2, 1, ROW_HEIGHT - 4);
              }
              i = Math.max(i + 1, next);
            }
            ctx.globalAlpha = 1;
          }
          if (!mergePixels) {
            for (const d of r.delList) {
              if (d.pos + d.len - 1 < fCol || d.pos > lCol) continue;
              const x1 = Math.max(
                labelWidth,
                Math.floor(labelWidth + d.pos * cellW - scrollLeft),
              );
              const x2 = Math.ceil(
                labelWidth + Math.min(d.pos + d.len - 1, lCol + 1) * cellW - scrollLeft,
              );
              const w = Math.max(1, Math.min(x2 - x1, labelWidth + seqAreaW - x1));
              if (w <= 0 || x1 > labelWidth + seqAreaW) continue;
              ctx.fillStyle = DEL_COLOR;
              ctx.fillRect(x1, y + 2, w, ROW_HEIGHT - 4);
            }
          } else {
            // Same density contract as the point events: a short del must not
            // paint a solid column for a pixel it barely touches — only a
            // column whose bases are mostly missing stays solid purple.
            ctx.fillStyle = DEL_COLOR;
            const cover = new Float32Array(Math.ceil(seqAreaW) + 1);
            for (const d of r.delList) {
              if (d.pos + d.len - 1 < fCol || d.pos > lCol) continue;
              const a = Math.max(d.pos, fCol);
              const b = Math.min(d.pos + d.len, lCol + 1); // half-open
              const pxA = Math.max(
                labelWidth,
                Math.floor(labelWidth + a * cellW - scrollLeft),
              );
              const pxB = Math.min(
                labelWidth + seqAreaW - 1,
                Math.floor(labelWidth + b * cellW - scrollLeft),
              );
              for (let px = pxA; px <= pxB; px++) {
                const colA = (px - labelWidth + scrollLeft) / cellW;
                const ov = Math.max(0, Math.min(b, colA + 1 / cellW) - Math.max(a, colA));
                if (ov > 0)
                  cover[px - labelWidth] = Math.min(1, cover[px - labelWidth] + ov * cellW);
              }
            }
            for (let px = 0; px < cover.length; px++) {
              if (cover[px] === 0) continue;
              ctx.globalAlpha = Math.max(MIN_DENSITY_INK, cover[px]);
              ctx.fillRect(labelWidth + px, y + 2, 1, ROW_HEIGHT - 4);
            }
            ctx.globalAlpha = 1;
          }
          const insLo = lowerBoundPos(r.insList, fCol);
          const insHi = lowerBoundPos(r.insList, lCol + 2);
          // same pixel merge for the 2px insertion ticks
          if (!mergePixels) {
            for (let j = insLo; j < insHi; j++) {
              const bx = Math.floor(labelWidth + r.insList[j].pos * cellW - scrollLeft);
              if (bx > labelWidth + seqAreaW) break;
              if (bx >= labelWidth) {
                ctx.fillStyle = INS_COLOR;
                ctx.fillRect(bx - 1, y + 2, 2, ROW_HEIGHT - 4);
              }
            }
          } else {
            let j = insLo;
            while (j < insHi) {
              const bx = Math.floor(labelWidth + r.insList[j].pos * cellW - scrollLeft);
              if (bx > labelWidth + seqAreaW) break;
              const nextPos = Math.ceil((bx + 2 - labelWidth + scrollLeft) / cellW);
              const next = lowerBoundPos(r.insList, nextPos);
              if (bx >= labelWidth) {
                ctx.globalAlpha = Math.min(
                  1,
                  Math.max(MIN_DENSITY_INK, (next - j) * cellW),
                );
                ctx.fillStyle = INS_COLOR;
                ctx.fillRect(bx - 1, y + 2, 2, ROW_HEIGHT - 4);
              }
              j = Math.max(j + 1, next);
            }
            ctx.globalAlpha = 1;
          }
        }
      }

      ctx.fillStyle = isDark ? "#27272a" : "#f4f4f5";
      ctx.fillRect(labelWidth, y + ROW_HEIGHT - 0.5, seqAreaW, 0.5);
    }

    /* ── sticky header: contig band + ruler with local numbers ── */
    ctx.fillStyle = isDark ? "#1e293b" : "#f8fafc";
    ctx.fillRect(labelWidth, 0, availableWidth - labelWidth, headerH);

    // contig separators spanning the full height
    for (const c of prepared.contigs) {
      if (c.start === 0) continue;
      const x = Math.floor(labelWidth + c.start * cellW - scrollLeft) + 0.5;
      if (x < labelWidth || x > labelWidth + seqAreaW) continue;
      ctx.strokeStyle = isDark ? "#334155" : "#cbd5e1";
      ctx.lineWidth = 1;
      ctx.beginPath();
      ctx.moveTo(x, 0);
      ctx.lineTo(x, canvasH);
      ctx.stroke();
    }

    // contig names in the band
    ctx.font = "9px ui-monospace, SFMono-Regular, monospace";
    ctx.textBaseline = "middle";
    ctx.textAlign = "center";
    for (const c of prepared.contigs) {
      const x1 = labelWidth + c.start * cellW - scrollLeft;
      const x2 = labelWidth + (c.start + c.len) * cellW - scrollLeft;
      if (x2 < labelWidth || x1 > labelWidth + seqAreaW) continue;
      const w = Math.min(x2, labelWidth + seqAreaW) - Math.max(x1, labelWidth);
      if (w < 30) continue;
      ctx.save();
      ctx.beginPath();
      ctx.rect(Math.max(x1, labelWidth), 0, w, CONTIG_BAND_H);
      ctx.clip();
      ctx.fillStyle = isDark ? "#94a3b8" : "#64748b";
      ctx.fillText(c.seqid, Math.max(x1, labelWidth) + w / 2, CONTIG_BAND_H / 2);
      ctx.restore();
    }

    // ruler with per-contig (local) numbers
    const rulerY = CONTIG_BAND_H;
    ctx.strokeStyle = isDark ? "#334155" : "#e2e8f0";
    ctx.lineWidth = 1;
    ctx.beginPath();
    ctx.moveTo(labelWidth, rulerY + RULER_HEIGHT - 0.5);
    ctx.lineTo(availableWidth, rulerY + RULER_HEIGHT - 0.5);
    ctx.stroke();

    const tickInterval = getPrettyStep(100, cellW);
    ctx.fillStyle = "#94a3b8";
    ctx.font = "9px ui-monospace, SFMono-Regular, monospace";
    ctx.textAlign = "center";
    ctx.textBaseline = "bottom";
    for (const c of prepared.contigs) {
      const cLo = Math.max(0, fCol - c.start, 0);
      const cHi = Math.min(c.len - 1, lCol - c.start);
      if (cHi < cLo) continue;
      const firstTick = Math.max(1, Math.ceil((cLo + 1) / tickInterval) * tickInterval);
      for (let p = firstTick; p <= cHi + 1; p += tickInterval) {
        const x = labelWidth + (c.start + p - 1) * cellW - scrollLeft;
        if (x < labelWidth || x > labelWidth + seqAreaW) continue;
        ctx.fillStyle = isDark ? "#334155" : "#cbd5e1";
        ctx.fillRect(x, rulerY + RULER_HEIGHT - 6, 1, 6);
        ctx.fillStyle = "#94a3b8";
        const label = p >= 1000 ? `${(p / 1000).toFixed(tickInterval >= 1000 ? 0 : 1)}k` : p;
        ctx.fillText(String(label), x, rulerY + RULER_HEIGHT - 7);
      }
    }

    /* ── zoom rectangle preview (right-drag) ── */
    if (zoomRect) {
      const s = Math.min(zoomRect.a0, zoomRect.a1);
      const e = Math.max(zoomRect.a0, zoomRect.a1);
      const selX1 = Math.max(labelWidth, labelWidth + s * cellW - scrollLeft);
      const selX2 = Math.min(labelWidth + seqAreaW, labelWidth + (e + 1) * cellW - scrollLeft);
      if (selX2 > selX1) {
        ctx.fillStyle = "rgba(74, 222, 128, 0.25)";
        ctx.fillRect(selX1, 0, selX2 - selX1, canvasH);
        ctx.strokeStyle = "#22c55e";
        ctx.lineWidth = 1.5;
        ctx.strokeRect(selX1 + 0.5, 0.5, selX2 - selX1 - 1, canvasH - 1);
      }
    }

    ctx.restore();

    /* ── labels (outside the clip) ── */
    ctx.fillStyle = isDark ? "#0f172a" : "#ffffff";
    ctx.fillRect(0, 0, labelWidth, canvasH);
    ctx.fillStyle = isDark ? "#334155" : "#e2e8f0";
    ctx.fillRect(labelWidth - 1, 0, 1, canvasH);

    ctx.textAlign = "right";
    ctx.textBaseline = "middle";
    ctx.font = "10px ui-sans-serif, system-ui, sans-serif";
    for (let row = 0; row < nRows; row++) {
      const y = headerH + row * ROW_HEIGHT + ROW_HEIGHT / 2;
      if (y > canvasH) break;
      const isRef = row === 0;
      // in letters mode the rev markers on the blocks are usually off
      // screen, so flag reverse-strand blocks on the row label instead
      const revVisible =
        !isRef &&
        viewMode === "letters" &&
        prepared.rows[row - 1].blocks.some(
          (b) => b.rev && b.e >= fCol && b.s <= lCol,
        );
      const text = isRef ? "Reference" : prepared.rows[row - 1].name;
      ctx.fillStyle = isRef
        ? isDark
          ? "#3b82f6"
          : "#2563eb"
        : isDark
          ? "#94a3b8"
          : "#64748b";
      ctx.font = `${isRef ? "bold " : ""}10px ui-sans-serif, system-ui, sans-serif`;
      const badge = revVisible ? " rev" : "";
      const badgeW = badge ? ctx.measureText(badge).width : 0;
      let display = text;
      if (ctx.measureText(display).width + badgeW > labelWidth - 20) {
        while (
          display.length > 5 &&
          ctx.measureText(`${display}\u2026`).width + badgeW > labelWidth - 20
        ) {
          display = display.slice(0, -1);
        }
        display += "\u2026";
      }
      if (badge) {
        ctx.fillStyle = isDark ? "#f59e0b" : "#b45309";
        ctx.fillText(badge, labelWidth - 8, y);
      }
      ctx.fillStyle = isRef
        ? isDark
          ? "#3b82f6"
          : "#2563eb"
        : isDark
          ? "#94a3b8"
          : "#64748b";
      ctx.fillText(display, labelWidth - 8 - badgeW, y);
    }

    /* size the hover overlay to match */
    const hoverCvs = hoverOverlayRef.current;
    if (hoverCvs) {
      hoverCvs.width = availableWidth * dpr;
      hoverCvs.height = canvasH * dpr;
      hoverCvs.style.width = `${availableWidth}px`;
      hoverCvs.style.height = `${canvasH}px`;
    }
  }, [
    prepared,
    availableWidth,
    totalH,
    scrollLeft,
    cellW,
    seqAreaW,
    anchorLen,
    viewMode,
    isDark,
    zoomRect,
    refBaseAt,
    refSeqVersion,
  ]);

  /* ── hover overlay: a thin vertical line (ported) ── */
  const drawHoverOverlay = useCallback(() => {
    const cvs = hoverOverlayRef.current;
    const mainCvs = canvasRef.current;
    if (!cvs || !mainCvs) return;
    const ctx = cvs.getContext("2d");
    if (!ctx) return;
    const dpr = window.devicePixelRatio || 1;
    if (cvs.width !== mainCvs.width || cvs.height !== mainCvs.height) {
      cvs.width = mainCvs.width;
      cvs.height = mainCvs.height;
    }
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    const w = cvs.width / dpr;
    const h = cvs.height / dpr;
    ctx.clearRect(0, 0, w, h);
    const hoverCol = hoverColRef.current;
    if (hoverCol === null || hoverCol < 0 || hoverCol >= anchorLen) return;
    // red when any query has a SNP at this anchor, blue otherwise
    let isSnp = false;
    if (prepared) {
      for (const r of prepared.rows) {
        if (snpAt(r.snpList, hoverCol)) {
          isSnp = true;
          break;
        }
      }
    }
    const hx = labelWidth + hoverCol * cellW - scrollLeft + cellW / 2;
    if (hx < labelWidth || hx > labelWidth + seqAreaW) return;
    ctx.strokeStyle = isSnp ? "#ef4444" : "#3b82f6";
    ctx.lineWidth = 1;
    ctx.beginPath();
    ctx.moveTo(Math.floor(hx) + 0.5, 0);
    ctx.lineTo(Math.floor(hx) + 0.5, h);
    ctx.stroke();
  }, [cellW, scrollLeft, seqAreaW, anchorLen, prepared, labelWidth]);

  useEffect(() => {
    draw();
  }, [draw]);
  useEffect(() => {
    drawHoverOverlay();
  }, [drawHoverOverlay]);
  useEffect(() => {
    redrawRef.current = () => drawHoverOverlay();
  }, [drawHoverOverlay]);

  /* ── hover tooltips ── */
  const handleCanvasMouseMove = useCallback(
    (e: React.MouseEvent<HTMLCanvasElement>) => {
      const cvs = canvasRef.current;
      if (!cvs || !prepared) return;
      if (panDragRef.current) return;
      const rect = cvs.getBoundingClientRect();
      const mouseXRaw = e.clientX - rect.left;
      const mouseYRaw = e.clientY - rect.top;

      if (mouseXRaw < labelWidth) {
        const row = Math.floor((mouseYRaw - headerH) / ROW_HEIGHT);
        if (row >= 0 && row <= prepared.rows.length && mouseYRaw >= headerH) {
          setLabelHover({
            text: row === 0 ? "Reference" : prepared.rows[row - 1].name,
            x: e.clientX,
            y: e.clientY,
          });
          cvs.style.cursor = "help";
        } else {
          setLabelHover(null);
          cvs.style.cursor = "default";
        }
        if (hoverColRef.current !== null) {
          hoverColRef.current = null;
          cancelAnimationFrame(hoverRafRef.current);
          hoverRafRef.current = requestAnimationFrame(() => redrawRef.current());
        }
        setHoverInfo(null);
        return;
      }
      setLabelHover(null);

      const col = Math.floor((scrollLeft + mouseXRaw - labelWidth) / cellW);
      const flipX = mouseXRaw > labelWidth + seqAreaW / 2;
      const row = Math.floor((mouseYRaw - headerH) / ROW_HEIGHT);
      const anchor = col >= 0 && col < anchorLen ? col : null;

      if (anchor !== null && row >= 0 && row <= prepared.rows.length && mouseYRaw >= headerH) {
        const contig = prepared.contigs.find(
          (c) => anchor >= c.start && anchor < c.start + c.len,
        );
        const local = contig ? anchor - contig.start + 1 : 0;
        const lines: string[] = [];
        let title = "";
        if (row === 0) {
          title = "Reference";
          const b = refBaseAt(anchor);
          lines.push(b ? `base ${b}` : "base (loading)");
        } else {
          const r = prepared.rows[row - 1];
          title = r.name;
          // insertion near the boundary?
          const tol = Math.max(4, Math.min(cellW, 10));
          const bxA = (labelWidth + anchor * cellW - scrollLeft) - mouseXRaw;
          const ins =
            Math.abs(bxA) <= tol
              ? insAt(r.insList, anchor)
              : insAt(r.insList, anchor + 1) !== null &&
                  Math.abs(labelWidth + (anchor + 1) * cellW - scrollLeft - mouseXRaw) <= tol
                ? insAt(r.insList, anchor + 1)
                : null;
          const del = delCovering(r.delList, anchor);
          const snp = ins || del ? null : snpAt(r.snpList, anchor);
          if (ins) {
            lines.push(
              `Insertion ${ins.seq.length} bp: ${ins.seq.length > 60 ? `${ins.seq.slice(0, 60)}…` : ins.seq}`,
            );
          } else if (del) {
            lines.push(`Deletion ${del.len} bp`);
          } else if (snp) {
            lines.push(
              `SNP ${String.fromCharCode(snp.r)} → ${String.fromCharCode(snp.q)}`,
            );
          } else if (inSpans(r.spans, anchor)) {
            lines.push("same as reference");
          } else {
            lines.push("not aligned");
          }
        }
        if (contig) {
          lines.push(`${contig.seqid}:${local.toLocaleString("en-US")}`);
        }
        setHoverInfo({ x: e.clientX, y: e.clientY, flipX, title, lines });
      } else {
        setHoverInfo(null);
      }

      cvs.style.cursor = viewMode === "letters" ? "pointer" : "crosshair";
      const newCol = anchor;
      if (newCol === hoverColRef.current) return;
      hoverColRef.current = newCol;
      cancelAnimationFrame(hoverRafRef.current);
      hoverRafRef.current = requestAnimationFrame(() => redrawRef.current());
    },
    [prepared, scrollLeft, cellW, seqAreaW, anchorLen, labelWidth, headerH, viewMode, refBaseAt],
  );

  const handleCanvasMouseLeave = useCallback(() => {
    setLabelHover(null);
    setHoverInfo(null);
    if (hoverColRef.current === null) return;
    hoverColRef.current = null;
    cancelAnimationFrame(hoverRafRef.current);
    hoverRafRef.current = requestAnimationFrame(() => redrawRef.current());
  }, []);

  /* ── left-drag = pan, left-click = center on the locus ── */
  const panDragRef = useRef<{ moved: boolean } | null>(null);
  const [isPanning, setIsPanning] = useState(false);

  const handleCanvasMouseDown = useCallback(
    (e: React.MouseEvent<HTMLCanvasElement>) => {
      if (!prepared) return;
      const cvs = canvasRef.current;
      if (!cvs) return;
      const rect = cvs.getBoundingClientRect();
      const mouseXCanvas = e.clientX - rect.left;
      // the label column carries hover and the resize handle, no gestures
      if (mouseXCanvas < labelWidth) return;

      if (e.button === 0) {
        const el = scrollRef.current;
        if (!el) return;
        const startX = e.clientX;
        const startY = e.clientY;
        const startSL = scrollLeft;
        const startST = el.scrollTop;
        panDragRef.current = { moved: false };
        setIsPanning(true);
        // drop the tooltip anchored to the pre-pan position
        setHoverInfo(null);
        setLabelHover(null);
        hoverColRef.current = null;
        redrawRef.current();
        document.body.style.userSelect = "none";

        const onMove = (ev: MouseEvent) => {
          if (!(ev.buttons & 1)) {
            onUp(ev);
            return;
          }
          const dx = ev.clientX - startX;
          const dy = ev.clientY - startY;
          if (Math.abs(dx) > 3 || Math.abs(dy) > 3) {
            panDragRef.current!.moved = true;
          }
          if (!panDragRef.current?.moved) return;
          const maxSL = Math.max(0, totalVirtualWRef.current - seqAreaWRef.current);
          const sl = Math.max(0, Math.min(maxSL, startSL - dx));
          el.scrollTop = startST - dy;
          setScrollLeft(sl);
        };
        const onUp = (ev: MouseEvent) => {
          document.removeEventListener("mousemove", onMove);
          document.removeEventListener("mouseup", onUp);
          document.body.style.userSelect = "";
          setIsPanning(false);
          const moved = panDragRef.current?.moved ?? false;
          panDragRef.current = null;
          if (moved) return;
          // A clean click (no drag) re-centres the view on that locus.
          const saw = seqAreaWRef.current;
          const tvw = totalVirtualWRef.current;
          const virtualX = startSL + ev.clientX - rect.left - labelWidth;
          const sl = Math.max(0, Math.min(Math.max(0, tvw - saw), virtualX - saw / 2));
          setScrollLeft(sl);
        };
        document.addEventListener("mousemove", onMove);
        document.addEventListener("mouseup", onUp);
        e.preventDefault();
        return;
      }

      /* ── right-drag = zoom to range (ported) ── */
      if (e.button !== 2) return;
      const startC = Math.max(
        0,
        Math.min(anchorLen - 1, Math.floor((scrollLeft + mouseXCanvas - labelWidth) / cellW)),
      );
      setZoomRect({ a0: startC, a1: startC });

      const colAt = (clientX: number) =>
        Math.max(
          0,
          Math.min(
            anchorLen - 1,
            Math.floor((scrollLeft + clientX - rect.left - labelWidth) / cellW),
          ),
        );

      const onUp = (ev: MouseEvent) => {
        document.removeEventListener("mousemove", onMove);
        document.removeEventListener("mouseup", onUp);
        setZoomRect(null);
        const endC = colAt(ev.clientX);
        const s = Math.min(startC, endC);
        const e2 = Math.max(startC, endC);
        if (e2 - s >= 5) {
          const newVF = Math.max(minVF, (e2 - s + 1) / anchorLen);
          const newTotalW = seqAreaW / newVF;
          const newSL = (s / anchorLen) * newTotalW;
          setViewFraction(newVF);
          setScrollLeft(newSL);
        }
      };
      const onMove = (ev: MouseEvent) => {
        if (!(ev.buttons & 2)) {
          onUp(ev);
          return;
        }
        setZoomRect({ a0: startC, a1: colAt(ev.clientX) });
      };
      document.addEventListener("mousemove", onMove);
      document.addEventListener("mouseup", onUp);
      e.preventDefault();
    },
    [prepared, anchorLen, scrollLeft, cellW, seqAreaW, labelWidth, minVF],
  );

  /* ── position bar drag (ported) ── */
  const scrollTrackRef = useRef<HTMLDivElement>(null);
  const handleBarMouseDown = (e: React.MouseEvent) => {
    if (!scrollTrackRef.current || !scrollRef.current) return;
    setIsBarDragging(true);
    const performMove = (clientX: number) => {
      if (!scrollTrackRef.current || !scrollRef.current) return;
      const rect = scrollTrackRef.current.getBoundingClientRect();
      const x = clientX - rect.left;
      const frac = Math.max(0, Math.min(1, x / rect.width));
      let targetStart = frac - viewFractionRef.current / 2;
      targetStart = Math.max(0, Math.min(1 - viewFractionRef.current, targetStart));
      const sl = targetStart * totalVirtualWRef.current;
      setScrollLeft(sl);
    };
    const onMouseMove = (ev: MouseEvent) => performMove(ev.clientX);
    const onMouseUp = () => {
      setIsBarDragging(false);
      window.removeEventListener("mousemove", onMouseMove);
      window.removeEventListener("mouseup", onMouseUp);
      document.body.style.userSelect = "";
    };
    document.body.style.userSelect = "none";
    window.addEventListener("mousemove", onMouseMove);
    window.addEventListener("mouseup", onMouseUp);
    performMove(e.clientX);
  };

  /* ── label column resize (ported) ── */
  const handleResizeMouseDown = (e: React.MouseEvent) => {
    setIsResizingLabel(true);
    const startX = e.clientX;
    const startW = labelWidth;
    const onMouseMove = (ev: MouseEvent) => {
      setLabelWidth(Math.max(80, Math.min(600, startW + ev.clientX - startX)));
    };
    const onMouseUp = () => {
      setIsResizingLabel(false);
      window.removeEventListener("mousemove", onMouseMove);
      window.removeEventListener("mouseup", onMouseUp);
      document.body.style.userSelect = "";
    };
    document.body.style.userSelect = "none";
    window.addEventListener("mousemove", onMouseMove);
    window.addEventListener("mouseup", onMouseUp);
  };

  if (error) return <p className="text-red-700 py-4 dark:text-red-400">{error}</p>;
  if (!prepared)
    return (
      <div className="flex items-center gap-3 text-zinc-500 py-16 justify-center dark:text-zinc-400">
        <Spinner /> Preparing the alignment view
        {data ? "" : " (variant events are computed on first use)"}...
      </div>
    );

  return (
    <div className="space-y-3">
      <div className="overflow-hidden rounded-xl border border-zinc-200 bg-white dark:border-zinc-800 dark:bg-zinc-900">
        {/* header */}
        <div className="flex items-center justify-between flex-wrap gap-2 border-b border-zinc-200 px-4 py-2 dark:border-zinc-800">
          <h2 className="text-sm font-semibold text-zinc-900 dark:text-zinc-100">
            Alignment{" "}
            <span className="font-mono text-[13px] font-normal text-zinc-500 dark:text-zinc-400">
              ({prepared.rows.length + 1} rows, {anchorLen.toLocaleString("en-US")} bp)
            </span>
          </h2>
          <div className="flex items-center gap-3">
            <div className="flex overflow-hidden rounded-md border border-zinc-300 dark:border-zinc-700">
              <button
                onClick={() => {
                  setViewMode("bars");
                  setViewFraction(1);
                  setScrollLeft(0);
                }}
                className={`px-3 py-1 text-[13px] font-medium transition-colors ${
                  viewMode === "bars"
                    ? "bg-zinc-100 text-zinc-900 dark:bg-zinc-800 dark:text-zinc-100"
                    : "bg-white text-zinc-600 hover:bg-zinc-50 dark:bg-zinc-900 dark:text-zinc-400 dark:hover:bg-zinc-700"
                }`}
              >
                Overview
              </button>
              <button
                onClick={() => {
                  const currentTotalW = seqAreaW / viewFraction;
                  const centerFrac = (scrollLeft + seqAreaW / 2) / currentTotalW;
                  setViewMode("letters");
                  const newVF = Math.min(1, 100 / anchorLen);
                  const newTotalW = seqAreaW / newVF;
                  const newSL = Math.max(
                    0,
                    Math.min(newTotalW - seqAreaW, centerFrac * newTotalW - seqAreaW / 2),
                  );
                  setViewFraction(newVF);
                  setScrollLeft(newSL);
                }}
                className={`border-l border-zinc-300 px-3 py-1 text-[13px] font-medium transition-colors dark:border-zinc-700 ${
                  viewMode === "letters"
                    ? "bg-zinc-100 text-zinc-900 dark:bg-zinc-800 dark:text-zinc-100"
                    : "bg-white text-zinc-600 hover:bg-zinc-50 dark:bg-zinc-900 dark:text-zinc-400 dark:hover:bg-zinc-700"
                }`}
              >
                Sequence
              </button>
            </div>
            <div className="flex items-center gap-1.5">
              <button
                onClick={zoomOut}
                className="flex h-6 w-6 items-center justify-center rounded border border-zinc-300 text-sm font-bold text-zinc-500 transition-colors hover:bg-zinc-100 hover:text-zinc-700 dark:border-zinc-700 dark:hover:bg-zinc-800"
                title="Zoom out"
              >
                &minus;
              </button>
              <input
                type="range"
                min={minVF}
                max={1}
                step={(1 - minVF) / 200}
                value={viewFraction}
                onChange={(e) => {
                  const newVF = parseFloat(e.target.value);
                  applyZoom(newVF / viewFraction, 0);
                }}
                className="h-1.5 w-24"
                style={{ direction: "rtl" }}
                aria-label="Zoom"
              />
              <button
                onClick={zoomIn}
                className="flex h-6 w-6 items-center justify-center rounded border border-zinc-300 text-sm font-bold text-zinc-500 transition-colors hover:bg-zinc-100 hover:text-zinc-700 dark:border-zinc-700 dark:hover:bg-zinc-800"
                title="Zoom in"
              >
                +
              </button>
            </div>
            <span className="w-20 text-right font-mono text-[13px] text-zinc-400">
              {Math.round(visibleBases).toLocaleString("en-US")} bp
            </span>
          </div>
        </div>

        {/* legend */}
        <div className="flex items-center gap-4 border-b border-zinc-100 bg-white px-4 py-1.5 text-[13px] text-zinc-500 dark:border-zinc-800 dark:bg-zinc-900 dark:text-zinc-400">
          <span className="flex items-center gap-1">
            <span
              className="inline-block h-3 w-3 rounded-sm"
              style={{ background: isDark ? "#3f3f46" : "#e4e4e7" }}
            />
            Sequence / Match
          </span>
          <span className="flex items-center gap-1">
            <span className="inline-block h-3 w-3 rounded-sm" style={{ background: SNP_COLOR }} />
            Mismatch
          </span>
          <span className="flex items-center gap-1">
            <span className="inline-block h-3 w-1" style={{ background: INS_COLOR }} />
            Insertion
          </span>
          <span className="flex items-center gap-1">
            <span className="inline-block h-3 w-3 rounded-sm" style={{ background: DEL_COLOR }} />
            Deletion
          </span>
          <span className="ml-auto italic text-zinc-400">
            Drag to pan, click to center. Ctrl/&#8984; + scroll to zoom.
            Right-drag to zoom to a range.
          </span>
        </div>

        {/* scroll arrow bar */}
        <div
          className="flex items-center border-b border-zinc-100 bg-white py-1 dark:border-zinc-800 dark:bg-zinc-900"
          style={{ paddingRight: `${RIGHT_PADDING}px` }}
        >
          <div
            style={{
              width: `${labelWidth}px`,
              flexShrink: 0,
              display: "flex",
              alignItems: "center",
              justifyContent: "flex-end",
              gap: "4px",
              paddingRight: "4px",
            }}
          >
            <button
              onMouseDown={() => startPan(-1)}
              onMouseUp={stopPan}
              onMouseLeave={stopPan}
              className="flex h-7 w-7 select-none items-center justify-center rounded border border-zinc-300 font-bold text-zinc-500 hover:bg-zinc-100 dark:border-zinc-700 dark:hover:bg-zinc-800"
              title="Scroll left"
            >
              &larr;
            </button>
            <button
              onMouseDown={() => startPan(1)}
              onMouseUp={stopPan}
              onMouseLeave={stopPan}
              className="flex h-7 w-7 select-none items-center justify-center rounded border border-zinc-300 font-bold text-zinc-500 hover:bg-zinc-100 dark:border-zinc-700 dark:hover:bg-zinc-800"
              title="Scroll right"
            >
              &rarr;
            </button>
          </div>
          <div
            ref={scrollTrackRef}
            onMouseDown={handleBarMouseDown}
            className="relative h-1.5 flex-1 cursor-pointer overflow-hidden rounded-full bg-zinc-100 dark:bg-zinc-800"
          >
            <div
              className={`h-full rounded-full bg-zinc-400 dark:bg-zinc-500 ${isBarDragging ? "duration-0" : "duration-75"}`}
              style={{ marginLeft: `${startFrac * 100}%`, width: `${viewFraction * 100}%` }}
            />
          </div>
        </div>

        {/* scrollable canvas area */}
        <div
          ref={setScrollEl}
          /* horizontal movement is virtual (state), so no wide spacer and no
             native horizontal scrollbar; vertical stays native for the rows */
          className="thin-scroll relative overscroll-contain overflow-x-hidden overflow-y-auto"
          style={{ height: `${Math.min(totalH, MAX_VIEWER_HEIGHT)}px` }}
        >
          <div
            style={{
              width: "100%",
              height: `${totalH}px`,
              position: "absolute",
              pointerEvents: "none",
            }}
          />
          <div style={{ position: "sticky", top: 0, left: 0, zIndex: 10, width: "fit-content" }}>
            <canvas
              ref={canvasRef}
              style={{
                display: "block",
                cursor: isPanning
                  ? "grabbing"
                  : viewMode === "letters"
                    ? "pointer"
                    : "crosshair",
              }}
              onMouseDown={handleCanvasMouseDown}
              onMouseMove={handleCanvasMouseMove}
              onMouseLeave={handleCanvasMouseLeave}
              onContextMenu={(e) => e.preventDefault()}
            />
            <canvas
              ref={hoverOverlayRef}
              style={{
                display: "block",
                position: "absolute",
                top: 0,
                left: 0,
                width: "100%",
                height: "100%",
                pointerEvents: "none",
                zIndex: 20,
              }}
            />
            <div
              onMouseDown={handleResizeMouseDown}
              className={`absolute top-0 bottom-0 z-30 w-1.5 cursor-col-resize transition-colors hover:bg-zinc-400/50 ${isResizingLabel ? "bg-zinc-400/60" : "bg-transparent"}`}
              style={{ left: `${labelWidth - 3}px` }}
            />
          </div>
        </div>
      </div>

      {/* label tooltip */}
      {labelHover && (
        <div
          className="pointer-events-none fixed z-50 whitespace-nowrap rounded border border-zinc-600 bg-zinc-800 px-2 py-1 text-[13px] text-white shadow-lg"
          style={{ left: `${labelHover.x + 12}px`, top: `${labelHover.y + 12}px` }}
        >
          {labelHover.text}
        </div>
      )}
      {/* position / variant tooltip */}
      {hoverInfo && (
        <div
          className="pointer-events-none fixed z-50 whitespace-nowrap rounded border border-zinc-600 bg-zinc-800 px-2 py-1 text-[13px] text-white shadow-lg"
          style={
            hoverInfo.flipX
              ? {
                  left: `${hoverInfo.x - 12}px`,
                  top: `${hoverInfo.y + 12}px`,
                  transform: "translateX(-100%)",
                }
              : { left: `${hoverInfo.x + 12}px`, top: `${hoverInfo.y + 12}px` }
          }
        >
          <div className="font-medium">{hoverInfo.title}</div>
          {hoverInfo.lines.map((l, i) => (
            <div key={i} className="font-mono text-[12px] text-zinc-300">
              {l}
            </div>
          ))}
        </div>
      )}
    </div>
  );
}
