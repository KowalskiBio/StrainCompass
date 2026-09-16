import { useEffect, useMemo, useRef, useState } from "react";
import type { WgaGene } from "../types";
import { usePopoverDismiss } from "./ui";

/** Matches shown at once; the ranking puts the ones worth seeing on top. */
const MAX_HITS = 8;

/**
 * How well a gene answers `q`, lower is better; -1 means it does not.
 *
 * The ranking is what makes typing a full locus tag feel like an exact lookup
 * while a vague word still finds something: an exact id beats a prefix, a
 * prefix beats a substring, and the function annotation is searched last so a
 * gene named for the query always outranks one merely described by it.
 */
function score(g: WgaGene, q: string): number {
  const locus = g.locus_tag.toLowerCase();
  const symbol = g.symbol?.toLowerCase() ?? "";
  const product = g.product?.toLowerCase() ?? "";
  if (locus === q || symbol === q) return 0;
  if (symbol.startsWith(q)) return 1;
  if (locus.startsWith(q)) return 2;
  if (symbol.includes(q)) return 3;
  if (locus.includes(q)) return 4;
  if (product.includes(q)) return 5;
  return -1;
}

/**
 * Find a gene anywhere in the reference and jump the map to it.
 *
 * Searching by locus tag, symbol or function annotation all land in the same
 * box because a user looking for "the katA gene" rarely knows which of those
 * the annotation actually calls it.
 */
export function GeneSearch({
  genes,
  multiContig,
  onPick,
}: {
  genes: WgaGene[];
  /** Show each hit's contig, which only disambiguates on a multi-contig reference. */
  multiContig: boolean;
  onPick: (g: WgaGene) => void;
}) {
  const [text, setText] = useState("");
  const [open, setOpen] = useState(false);
  const [cursor, setCursor] = useState(0);
  const inputRef = useRef<HTMLInputElement>(null);
  const ref = usePopoverDismiss(open, () => setOpen(false));

  const hits = useMemo(() => {
    const q = text.trim().toLowerCase();
    if (!q) return [];
    const scored: { g: WgaGene; s: number }[] = [];
    for (const g of genes) {
      const s = score(g, q);
      if (s >= 0) scored.push({ g, s });
    }
    // Ties are broken by position so the list is stable and reads in genome
    // order rather than in whatever order the annotation happened to load.
    scored.sort((a, b) => a.s - b.s || a.g.start - b.g.start);
    return scored.slice(0, MAX_HITS).map((x) => x.g);
  }, [genes, text]);

  useEffect(() => setCursor(0), [text]);

  function pick(g: WgaGene | undefined) {
    if (!g) return;
    onPick(g);
    setOpen(false);
    // The locus tag, not the symbol: it is the one name every gene has and it
    // matches exactly, so searching the box's own contents finds this gene again.
    setText(g.locus_tag);
    inputRef.current?.blur();
  }

  return (
    <div className="relative" ref={ref}>
      <input
        ref={inputRef}
        value={text}
        onChange={(e) => {
          setText(e.target.value);
          setOpen(true);
        }}
        onFocus={() => setOpen(true)}
        onKeyDown={(e) => {
          if (e.key === "ArrowDown") {
            e.preventDefault();
            setOpen(true);
            setCursor((c) => Math.min(hits.length - 1, c + 1));
          } else if (e.key === "ArrowUp") {
            e.preventDefault();
            setCursor((c) => Math.max(0, c - 1));
          } else if (e.key === "Enter") {
            e.preventDefault();
            pick(hits[cursor]);
          } else if (e.key === "Escape") {
            setOpen(false);
          }
        }}
        placeholder="Find a gene..."
        aria-label="Find a gene by locus tag, symbol or function"
        className="h-11 w-52 px-3 rounded-lg border border-zinc-300 text-[15px] dark:border-zinc-700 dark:bg-zinc-900"
      />
      {open && text.trim() !== "" && (
        <div className="absolute left-0 top-12 z-30 w-80 bg-white border border-zinc-200 rounded-lg shadow-lg p-1 dark:bg-zinc-900 dark:border-zinc-800">
          {hits.length === 0 ? (
            <p className="px-3 py-2 text-sm text-zinc-400 dark:text-zinc-500">
              No gene matches that.
            </p>
          ) : (
            hits.map((g, i) => (
              <button
                key={`${g.seqid}:${g.locus_tag}`}
                // mousedown, not click: the dismiss listener fires on mousedown
                // and would close the list before a click could land.
                onMouseDown={(e) => {
                  e.preventDefault();
                  pick(g);
                }}
                onMouseEnter={() => setCursor(i)}
                className={`block w-full text-left px-3 py-1.5 rounded-md ${
                  i === cursor ? "bg-zinc-100 dark:bg-zinc-800" : ""
                }`}
              >
                <span className="flex items-baseline gap-2">
                  <span className="font-mono text-sm truncate">{g.locus_tag}</span>
                  {g.symbol && g.symbol !== g.locus_tag && (
                    <span className="text-sm font-medium truncate">{g.symbol}</span>
                  )}
                  <span className="ml-auto shrink-0 text-xs text-zinc-400 font-mono dark:text-zinc-500">
                    {multiContig ? `${g.seqid} ` : ""}
                    {Math.round(g.start).toLocaleString("en-US")}
                  </span>
                </span>
                {g.product && (
                  <span className="block text-xs text-zinc-400 truncate dark:text-zinc-500">
                    {g.product}
                  </span>
                )}
              </button>
            ))
          )}
        </div>
      )}
    </div>
  );
}
