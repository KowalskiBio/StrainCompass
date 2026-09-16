import type { Range } from "./genomeRange";

/**
 * Wrapping the strain map onto stacked rows.
 *
 * The visible window stays one contiguous range; the rows only slice it, the
 * way a paragraph wraps. Row r covers [start + r*rowSpan, start + (r+1)*rowSpan)
 * and continues where row r-1 ended, so N rows draw the same window at N times
 * the horizontal resolution.
 */

/** Per-row heights in px. Slim, so three or four rows fit an ordinary window. */
export type RowMetrics = { rulerH: number; laneH: number; geneH: number };
export const ROW_METRICS: RowMetrics = { rulerH: 22, laneH: 60, geneH: 32 };

/** Height of one row: its own little ruler plus the lane the beads sit in. */
export const ROW_PITCH = ROW_METRICS.rulerH + ROW_METRICS.laneH;

export const MIN_ROWS = 1;
export const MAX_ROWS = 4;

/** Page left for the legend and the caption under the map, in px. */
const BOTTOM_RESERVE = 110;
/** Slack needed before the count grows, so a resize cannot flip it per pixel. */
const GROW_SLACK = 0.25;

/**
 * How many rows fit under a map whose top edge is at `top` in the viewport.
 * `current` is the count in use, kept when the difference is within the slack.
 */
export function rowsForHeight(
  top: number,
  innerHeight: number,
  current: number,
  pitch: number = ROW_PITCH,
): number {
  const raw = (innerHeight - top - BOTTOM_RESERVE) / pitch;
  if (!Number.isFinite(raw)) return current;
  const next =
    raw >= current + 1 + GROW_SLACK || raw < current ? Math.floor(raw) : current;
  return Math.min(MAX_ROWS, Math.max(MIN_ROWS, next));
}

export type WrapLayout = {
  rows: number;
  width: number;
  height: number;
  rowPitch: number;
  /** Bases one row covers. The scale everything on screen is really drawn at. */
  rowSpan: number;
  start: number;
  end: number;
  rowStart(row: number): number;
  rowEnd(row: number): number;
  rowY(row: number): number;
  baselineY(row: number): number;
  /** The row a base falls on, clamped to the stack. */
  rowOfBp(bp: number): number;
  /** x of a base within its own row. */
  xOfBp(bp: number): number;
  /** x of a base measured in `row`, which need not be the row it falls on. */
  xInRow(bp: number, row: number): number;
  rowAtY(y: number): number;
  /** The base under a point in the SVG. The inverse of xOfBp/baselineY. */
  bpAtPoint(x: number, y: number): number;
};

export function makeWrapLayout(
  range: Range,
  rows: number,
  width: number,
  m: RowMetrics = ROW_METRICS,
): WrapLayout {
  const n = Math.max(1, Math.round(rows));
  const w = Math.max(1, width);
  const rowSpan = Math.max(1, range.end - range.start) / n;
  const rowPitch = m.rulerH + m.laneH;
  const clampRow = (r: number) => Math.min(n - 1, Math.max(0, r));
  const xInRow = (bp: number, row: number) =>
    ((bp - range.start) / rowSpan - row) * w;
  const rowOfBp = (bp: number) => clampRow(Math.floor((bp - range.start) / rowSpan));

  return {
    rows: n,
    width: w,
    height: n * rowPitch,
    rowPitch,
    rowSpan,
    start: range.start,
    end: range.end,
    rowStart: (row) => range.start + row * rowSpan,
    rowEnd: (row) => range.start + (row + 1) * rowSpan,
    rowY: (row) => row * rowPitch,
    baselineY: (row) => row * rowPitch + m.rulerH + m.laneH / 2,
    rowOfBp,
    xOfBp: (bp) => xInRow(bp, rowOfBp(bp)),
    xInRow,
    rowAtY: (y) => clampRow(Math.floor(y / rowPitch)),
    bpAtPoint: (x, y) =>
      range.start + (clampRow(Math.floor(y / rowPitch)) + x / w) * rowSpan,
  };
}

/** One row's share of a feature that may straddle row boundaries. */
export type RowPiece = {
  row: number;
  x: number;
  w: number;
  /** The piece's own bounds, to tell a cut edge from the feature's real end. */
  from: number;
  to: number;
};

/**
 * Split [from, to) into one piece per row it crosses, clipped to the window.
 * A feature shorter than a base still yields one zero-width piece, which the
 * caller widens to its minimum mark.
 */
export function rowPieces(l: WrapLayout, from: number, to: number): RowPiece[] {
  const a = Math.max(from, l.start);
  const b = Math.min(to, l.end);
  if (b < a) return [];
  const first = l.rowOfBp(a);
  // A feature ending exactly on a boundary belongs to the row above it, not to
  // the next row it would otherwise open with a zero-width piece.
  const last = l.rowOfBp(Math.max(a, b - l.rowSpan * 1e-9));
  const out: RowPiece[] = [];
  for (let row = first; row <= last; row++) {
    const pf = Math.max(a, l.rowStart(row));
    const pt = Math.min(b, l.rowEnd(row));
    const x = l.xInRow(pf, row);
    out.push({ row, x, w: l.xInRow(pt, row) - x, from: pf, to: pt });
  }
  return out;
}
