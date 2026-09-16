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

// Zoom per pixel of pinch travel, and the most one event may zoom.
const ZOOM_PER_PX = 0.0035;
const MAX_STEP = 2;

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

/**
 * Continuous zoom factor from a pixel delta. Positive delta (pinch in, wheel down)
 * gives a factor above 1, which widens the visible window: zoom out.
 */
export function zoomFactor(dy: number): number {
  return Math.min(MAX_STEP, Math.max(1 / MAX_STEP, Math.exp(dy * ZOOM_PER_PX)));
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
      cb.current({ kind: "zoom", factor: zoomFactor(dy), clientX: e.clientX });
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
      const factor = Math.min(MAX_STEP, Math.max(1 / MAX_STEP, lastScale / scale));
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
