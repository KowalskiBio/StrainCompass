/** A visible window on a sequence, in 1-based base pairs. */
export type Range = { start: number; end: number };

/** Closest zoom, in base pairs across the window. */
export const MIN_SPAN = 40;

/** Clamp a window of `span` bp so it stays inside the sequence. */
export function clampRange(start: number, span: number, seqLength: number): Range {
  const s = Math.min(Math.max(1, start), Math.max(1, seqLength - span));
  return { start: Math.round(s), end: Math.round(s + span) };
}

/**
 * Zoom by `factor` while keeping `anchorBp` under the same point on screen.
 * Anchoring to the cursor's *fraction* rather than re-centering on it is what
 * stops the view sliding sideways on every step.
 */
export function zoomRange(
  base: Range,
  factor: number,
  anchorBp: number,
  seqLength: number,
): Range {
  const curSpan = base.end - base.start;
  const span = Math.min(seqLength, Math.max(MIN_SPAN, curSpan * factor));
  const frac = curSpan > 0 ? (anchorBp - base.start) / curSpan : 0.5;
  return clampRange(anchorBp - frac * span, span, seqLength);
}

/** Slide the window by `dxBp` base pairs, keeping its span. */
export function panRange(base: Range, dxBp: number, seqLength: number): Range {
  const span = base.end - base.start;
  return clampRange(base.start + dxBp, span, seqLength);
}
