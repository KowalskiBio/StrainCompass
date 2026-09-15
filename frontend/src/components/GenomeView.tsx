import { useEffect, useMemo, useRef, useState } from "react";
import { api } from "../api";
import type { Run, WgaBlock, WgaData, WgaGene } from "../types";
import { Spinner } from "./ui";

/**
 * Whole genome alignment view: reference ruler, gene track, and one
 * alignment track per query. Zoom with the wheel, pan by dragging, or
 * set the range by typing coordinates.
 */
export function GenomeView({
  run,
  initialRange,
  initialGene,
  onOpenGene,
  onRangeChange,
}: {
  run: Run;
  initialRange?: { seqid: string; start: number; end: number } | null;
  initialGene?: string | null;
  onOpenGene: (locus: string) => void;
  onRangeChange: (r: { seqid: string; start: number; end: number }) => void;
}) {
  const [data, setData] = useState<WgaData | null>(null);
  const [error, setError] = useState<string | null>(null);
  const [seqid, setSeqid] = useState(initialRange?.seqid ?? "");
  const [range, setRange] = useState<{ start: number; end: number } | null>(
    initialRange ? { start: initialRange.start, end: initialRange.end } : null,
  );
  const svgRef = useRef<SVGSVGElement>(null);
  const initialGeneRef = useRef(initialGene);
  const [popup, setPopup] = useState<{
    x: number;
    y: number;
    gene?: WgaGene;
    block?: WgaBlock & { queryName: string };
  } | null>(null);
  const dragRef = useRef<{ x: number; start: number } | null>(null);

  useEffect(() => {
    let cancelled = false;
    api
      .wga(run.id)
      .then((d) => {
        if (cancelled) return;
        setData(d);
        setSeqid((prev) => prev || (d.reference[0]?.[0] ?? ""));
      })
      .catch((e) => !cancelled && setError((e as Error).message));
    return () => {
      cancelled = true;
    };
  }, [run.id]);

  const seqLength = useMemo(
    () => data?.reference.find((r) => r[0] === seqid)?.[1] ?? 0,
    [data, seqid],
  );

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
  if (error) return <p className="text-red-700 py-4 dark:text-red-400">{error}</p>;
  if (!data || !range || !seqid)
    return (
      <div className="flex items-center gap-3 text-zinc-500 py-16 justify-center dark:text-zinc-400">
        <Spinner /> Loading the genome view...
      </div>
    );

  const genes = data.genes.filter((g) => g.seqid === seqid);
  const queries = data.queries.map((q) => ({
    name: q.query_name,
    blocks: q.blocks.filter((b) => b.ref_seqid === seqid),
  }));

  const width = 1100;
  const rulerH = 28;
  const geneTrackH = 34;
  const trackH = 34;
  const height = rulerH + geneTrackH + queries.length * trackH + 8;
  const bpToX = (bp: number) =>
    ((bp - range.start) / Math.max(1, range.end - range.start)) * width;

  function zoom(factor: number, anchorBp: number) {
    setRange((r) => {
      if (!r) return r;
      const half = ((r.end - r.start) / 2) * factor;
      let start = anchorBp - half;
      let end = anchorBp + half;
      if (start < 1) {
        start = 1;
        end = Math.min(seqLength, 1 + half * 2);
      }
      if (end > seqLength) {
        end = seqLength;
        start = Math.max(1, seqLength - half * 2);
      }
      if (end - start < 40) return r;
      return { start: Math.round(start), end: Math.round(end) };
    });
  }

  const bpAt = (clientX: number) => {
    const rect = svgRef.current!.getBoundingClientRect();
    const frac = (clientX - rect.left) / rect.width;
    return range.start + frac * (range.end - range.start);
  };

  // ticks for the ruler
  const span = range.end - range.start;
  const tickStep = niceStep(span);
  const firstTick = Math.ceil(range.start / tickStep) * tickStep;
  const ticks: number[] = [];
  for (let t = firstTick; t <= range.end; t += tickStep) ticks.push(t);

  const inRange = (s: number, e: number) => e >= range.start && s <= range.end;

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
      </div>

      <div
        className="border border-zinc-200 rounded-xl bg-white overflow-x-auto dark:border-zinc-800 dark:bg-zinc-900"
        onMouseLeave={() => setPopup(null)}
      >
        <svg
          ref={svgRef}
          width={width}
          height={height}
          viewBox={`0 0 ${width} ${height}`}
          className="block select-none"
          style={{ cursor: dragRef.current ? "grabbing" : "grab", minWidth: width }}
          onWheel={(e) => {
            e.preventDefault();
            const bp = bpAt(e.clientX);
            zoom(e.deltaY > 0 ? 1.2 : 0.83, bp);
          }}
          onMouseDown={(e) => {
            dragRef.current = { x: e.clientX, start: range.start };
          }}
          onMouseMove={(e) => {
            if (dragRef.current && svgRef.current) {
              const rect = svgRef.current.getBoundingClientRect();
              const shiftBp =
                ((dragRef.current.x - e.clientX) / rect.width) *
                (range.end - range.start);
              const start = Math.min(
                Math.max(1, dragRef.current.start + shiftBp),
                seqLength - (range.end - range.start),
              );
              setRange({ start, end: start + (range.end - range.start) });
            }
          }}
          onMouseUp={() => (dragRef.current = null)}
          onMouseLeave={() => (dragRef.current = null)}
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

          {/* gene track */}
          <g>
            {genes
              .filter((g) => inRange(g.start, g.end))
              .map((g) => {
                const x = bpToX(Math.max(g.start, range.start));
                const x2 = bpToX(Math.min(g.end, range.end));
                const w = Math.max(2, x2 - x);
                return (
                  <rect
                    key={g.locus_tag}
                    x={x}
                    y={rulerH + 6}
                    width={w}
                    height={geneTrackH - 14}
                    rx={1.5}
                    style={{ fill: geneColor(g) }}
                    className="cursor-pointer"
                    onClick={(e) => {
                      e.stopPropagation();
                      const rect = svgRef.current!.getBoundingClientRect();
                      setPopup({
                        x: e.clientX - rect.left,
                        y: rulerH + 20,
                        gene: g,
                      });
                    }}
                  />
                );
              })}
            <text x={4} y={rulerH + 2} fontSize="10" style={{ fill: "var(--gv-track-label)" }}>
              genes
            </text>
          </g>

          {/* query tracks */}
          {queries.map((q, i) => {
            const y = rulerH + geneTrackH + i * trackH + 4;
            return (
              <g key={q.name}>
                <rect
                  x={0}
                  y={y}
                  width={width}
                  height={trackH - 10}
                  style={{ fill: "var(--gv-track-bg)" }}
                />
                <text x={4} y={y + 9} fontSize="10" style={{ fill: "var(--gv-track-label)" }}>
                  {i === 0 ? "alignments: " : ""}
                </text>
                {q.blocks
                  .filter((b) => inRange(b.ref_start, b.ref_end))
                  .map((b, j) => {
                    const x = bpToX(Math.max(b.ref_start, range.start));
                    const x2 = bpToX(Math.min(b.ref_end, range.end));
                    const w = Math.max(2, x2 - x);
                    return (
                      <g key={j}>
                        <rect
                          x={x}
                          y={y + 12}
                          width={w}
                          height={trackH - 22}
                          fill={identityColor(b.identity)}
                          className="cursor-pointer"
                          onClick={(e) => {
                            e.stopPropagation();
                            setPopup({
                              x: e.clientX - svgRef.current!.getBoundingClientRect().left,
                              y: y + trackH,
                              block: { ...b, queryName: q.name },
                            });
                          }}
                        />
                        {b.qry_rev && w > 12 && (
                          <text
                            x={x + w / 2}
                            y={y + 10}
                            fontSize="8"
                            style={{ fill: "var(--gv-tick-label)" }}
                            textAnchor="middle"
                          >
                            rev
                          </text>
                        )}
                      </g>
                    );
                  })}
              </g>
            );
          })}
        </svg>
      </div>

      {/* legend */}
      <div className="flex flex-wrap items-center gap-4 text-xs text-zinc-500 dark:text-zinc-400">
        <span className="inline-flex items-center gap-1.5">
          <span className="w-3 h-3 rounded-sm" style={{ background: identityColor(100) }} />
          high identity
        </span>
        <span className="inline-flex items-center gap-1.5">
          <span className="w-3 h-3 rounded-sm" style={{ background: identityColor(85) }} />
          medium
        </span>
        <span className="inline-flex items-center gap-1.5">
          <span className="w-3 h-3 rounded-sm" style={{ background: identityColor(70) }} />
          low identity
        </span>
        <span className="inline-flex items-center gap-1.5">
          <span className="w-3 h-3 rounded-sm" style={{ background: "var(--gv-ruler-line)" }} />
          not aligned
        </span>
        <span>Scroll to zoom, drag to pan. Click a gene or a block for details.</span>
      </div>

      {popup && (
        <GenePopup
          popup={popup}
          onClose={() => setPopup(null)}
          onOpenGene={(locus) => {
            setPopup(null);
            onOpenGene(locus);
          }}
        />
      )}
    </div>
  );
}

function GenePopup({
  popup,
  onClose,
  onOpenGene,
}: {
  popup: { x: number; y: number; gene?: WgaGene; block?: WgaBlock & { queryName: string } };
  onClose: () => void;
  onOpenGene: (locus: string) => void;
}) {
  useEffect(() => {
    const t = setTimeout(() => {}, 0);
    return () => clearTimeout(t);
  }, []);
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
      style={{ left: Math.min(popup.x, 800), top: popup.y + 40 }}
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
          <button
            className="mt-3 w-full h-10 rounded-md bg-zinc-900 text-white text-sm hover:bg-zinc-700 dark:bg-zinc-100 dark:text-zinc-900 dark:hover:bg-zinc-300"
            onClick={() => onOpenGene(popup.gene!.locus_tag)}
          >
            Show alignment
          </button>
        </>
      )}
      {popup.block && (
        <>
          <p className="font-semibold text-[15px]">Alignment block</p>
          <p className="text-xs text-zinc-500 mt-1 dark:text-zinc-400">Query: {popup.block.queryName}</p>
          <p className="text-xs text-zinc-500 mt-1 font-mono dark:text-zinc-400">
            reference {popup.block.ref_start.toLocaleString("en-US")} -{" "}
            {popup.block.ref_end.toLocaleString("en-US")}
          </p>
          <p className="text-xs text-zinc-500 mt-1 font-mono dark:text-zinc-400">
            query {popup.block.qry_seqid}:{" "}
            {popup.block.qry_start.toLocaleString("en-US")} -{" "}
            {popup.block.qry_end.toLocaleString("en-US")}
            {popup.block.qry_rev ? " (reverse strand)" : ""}
          </p>
          <p className="text-xs text-zinc-500 mt-1 font-mono dark:text-zinc-400">
            identity {popup.block.identity.toFixed(2)}%
          </p>
        </>
      )}
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

/** Sequential single-hue scale for identity 0-100. */
function identityColor(identity: number): string {
  const t = Math.max(0, Math.min(100, identity)) / 100;
  // from light steel blue to deep blue
  const stops: [number, number, number][] = [
    [226, 232, 240],
    [145, 175, 212],
    [30, 95, 158],
  ];
  let c: [number, number, number];
  if (t < 0.5) {
    const f = t / 0.5;
    c = mix(stops[0], stops[1], f);
  } else {
    const f = (t - 0.5) / 0.5;
    c = mix(stops[1], stops[2], f);
  }
  return `rgb(${c.map((x) => Math.round(x)).join(",")})`;
}

function mix(
  a: [number, number, number],
  b: [number, number, number],
  t: number,
): [number, number, number] {
  return [a[0] + (b[0] - a[0]) * t, a[1] + (b[1] - a[1]) * t, a[2] + (b[2] - a[2]) * t];
}
