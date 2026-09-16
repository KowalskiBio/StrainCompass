import { useCallback, useEffect, useRef, useState } from "react";

/**
 * Wheel/pinch gestures with a listener React cannot give us.
 *
 * React registers `wheel` at the root as a *passive* listener and offers no way to
 * opt out, so `e.preventDefault()` inside an `onWheel` prop is always a no-op. On
 * macOS a trackpad pinch arrives as a wheel event with `ctrlKey` set, and leaving
 * its default action in place lets the browser zoom the whole page. The only fix is
 * to attach the listener ourselves with `{ passive: false }`.
 */
export type WheelGesture =
  | { kind: "zoom"; factor: number; clientX: number }
  | { kind: "pan"; dx: number };

/** What produced a wheel event. */
export type WheelSource = "pinch" | "wheel" | "scroll";

// deltaMode 1/2 report lines/pages; these bring them back to rough pixels.
const LINE_PX = 16;
const PAGE_PX = 400;

/**
 * Zoom gain per pixel of travel, and the most one event may zoom.
 *
 * Pinch and wheel need different gains: a trackpad pinch arrives as a stream of
 * small deltas (a few pixels each, ~60/s), while a mouse wheel sends a handful of
 * coarse notches (commonly 100-120 px each). One shared constant high enough to
 * make pinching feel quick would make every wheel notch jump the full clamp.
 */
// A half-second trackpad pinch is roughly 120 px of travel, so ln(20)/120 puts
// one pinch at about 20x: whole genome to ~145 kb, gene level in three.
const PINCH_PER_PX = 0.025;
// About 2.6x per wheel notch, in the same spirit as the 2x +/- buttons.
const WHEEL_PER_PX = 0.008;
const MAX_STEP = 6;

/** Extra gain on Safari's pinch scale, to match the wheel path's feel. */
const GESTURE_GAIN = 1.8;

/**
 * There is no browser API for telling a mouse wheel from a trackpad, so this is a
 * heuristic: mouse wheels emit coarse, whole-number, single-axis deltas (commonly
 * +/-100 or +/-120), trackpads emit fine and often fractional ones on both axes.
 * Kept in one place so the thresholds can be tuned without touching the views.
 */
export function classifyWheel(e: WheelEvent): WheelSource {
  if (e.ctrlKey || e.metaKey) return "pinch";
  if (e.deltaMode !== 0) return "wheel";
  if (e.deltaX === 0 && Math.abs(e.deltaY) >= 40 && Number.isInteger(e.deltaY))
    return "wheel";
  return "scroll";
}

/** Wheel deltas in pixels, whatever unit the event used. */
export function pixelDelta(e: WheelEvent): { dx: number; dy: number } {
  const unit = e.deltaMode === 1 ? LINE_PX : e.deltaMode === 2 ? PAGE_PX : 1;
  return { dx: e.deltaX * unit, dy: e.deltaY * unit };
}

/** Keep a zoom factor within one event's allowed range. */
function clampStep(factor: number): number {
  return Math.min(MAX_STEP, Math.max(1 / MAX_STEP, factor));
}

/**
 * Continuous zoom factor from a pixel delta. Positive delta (pinch in, wheel down)
 * gives a factor above 1, which widens the visible window: zoom out.
 */
export function zoomFactor(dy: number, source: WheelSource = "pinch"): number {
  const gain = source === "wheel" ? WHEEL_PER_PX : PINCH_PER_PX;
  return clampStep(Math.exp(dy * gain));
}

/** Safari reports trackpad pinch as GestureEvent, not as a ctrlKey wheel event. */
const hasGestureEvents =
  typeof window !== "undefined" && "ongesturestart" in window;

type SafariGestureEvent = Event & { scale: number; clientX: number };

/**
 * Returns a callback ref. Attach it to the element that should own the gestures.
 *
 * `capturePan: false` leaves unmodified scrolling to the browser, for elements that
 * already scroll natively.
 */
export function useWheelGestures<T extends HTMLElement>(
  onGesture: (g: WheelGesture) => void,
  opts: { capturePan?: boolean } = {},
): (node: T | null) => void {
  const capturePan = opts.capturePan ?? true;

  // A callback ref rather than a RefObject: it re-runs the effect when the element
  // actually mounts, which a RefObject silently fails to do when the component
  // renders a placeholder on its first pass.
  const [node, setNode] = useState<T | null>(null);

  const cb = useRef(onGesture);
  useEffect(() => {
    cb.current = onGesture;
  });

  useEffect(() => {
    if (!node) return;

    const onWheel = (e: WheelEvent) => {
      const source = classifyWheel(e);
      const { dx, dy } = pixelDelta(e);

      if (source === "scroll") {
        if (!capturePan) return;
        e.preventDefault();
        // A vertical two-finger swipe pans too: the map has only one axis to move on.
        cb.current({ kind: "pan", dx: dx || dy });
        return;
      }

      // Always prevent, even when Safari will deliver the gesture separately below:
      // this is what stops the browser zooming the page.
      e.preventDefault();
      if (source === "pinch" && hasGestureEvents) return;
      cb.current({ kind: "zoom", factor: zoomFactor(dy, source), clientX: e.clientX });
    };

    node.addEventListener("wheel", onWheel, { passive: false });

    if (!hasGestureEvents) {
      return () => node.removeEventListener("wheel", onWheel);
    }

    let lastScale = 1;
    const onGestureStart = (e: Event) => {
      e.preventDefault();
      lastScale = (e as SafariGestureEvent).scale || 1;
    };
    const onGestureChange = (e: Event) => {
      e.preventDefault();
      const ev = e as SafariGestureEvent;
      const scale = ev.scale || 1;
      if (scale <= 0 || lastScale <= 0) return;
      // Growing scale means pinch out, which narrows the window: factor below 1.
      const factor = clampStep((lastScale / scale) ** GESTURE_GAIN);
      lastScale = scale;
      cb.current({ kind: "zoom", factor, clientX: ev.clientX });
    };
    const onGestureEnd = (e: Event) => e.preventDefault();

    node.addEventListener("gesturestart", onGestureStart);
    node.addEventListener("gesturechange", onGestureChange);
    node.addEventListener("gestureend", onGestureEnd);

    return () => {
      node.removeEventListener("wheel", onWheel);
      node.removeEventListener("gesturestart", onGestureStart);
      node.removeEventListener("gesturechange", onGestureChange);
      node.removeEventListener("gestureend", onGestureEnd);
    };
  }, [node, capturePan]);

  return useCallback((n: T | null) => setNode(n), []);
}
