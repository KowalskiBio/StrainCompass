# Genome View Overhaul: StrainMap + Alignment

## Goal

Redesign the "Genome view" tab with two user-selectable sub-modes:

- **StrainMap**: a BioCyc-style "beads on a string" gene map filling the whole
  rectangle. Genes are rectangles on a line, arrow direction = strand, colored by the
  presence call of the selected query (present = green, partial = amber, absent = red,
  n/a = gray). The line is **wrapped onto several rows** (as many as the page height
  allows, up to 4): row r continues where row r-1 ended, the way a paragraph wraps, so
  the same window is drawn at as many times the horizontal resolution as there are rows.
  Each row carries its own ruler and baseline, and a gene cut by a row boundary is drawn
  once per row with a flat edge at the cut. The old blue identity-alignment rows are
  removed. Gene labels (symbol, fallback locus_tag) are drawn inside the boxes when
  zoomed in enough. Variant markers (SNP = red tick, deletion = purple tick,
  insertion = blue tick) are drawn on the gene once a *row* covers less than ~20 kb.
  Positions shown as a ruler only. Every zoom and detail threshold is per row, so the
  density on screen is the same whatever the row count. Zoom and pan are unchanged:
  the window stays one contiguous range, and dragging a row's width across moves the
  window by what one row covers.
  Inspiration: https://biocyc.org/genbro/genbro.shtml?orgid=10403S_RAST&replicon=CDF
- **Alignment**: an Oligool-style MSA viewer. One concatenated reference coordinate space
  (all contigs joined, separator lines + contig names), a reference row plus one row per
  query FASTA. Colors: mismatch = red `#dc2626`, insertion = blue `#3b82f6`,
  deletion = purple `#9333ea`, match = gray. Zoom from whole genome down to base level;
  letters mode below ~85 visible bases with hysteresis at 100. Reverse-strand query
  bases are shown revcomp'd with a "rev" badge.

## Data design

- Variant data (SNPs, deletions, insertions per query per reference contig) is extracted
  from the kept `work/cmp.delta` nucmer files.
- **Precomputed by `run_comparison`** and written as `{run_dir}/queries/{qid}/variants.json`
  at run end, so the viewers never walk the delta on first open. Runs made before this
  still backfill on demand (parallel per query) and then serve from the same cache.
- Whole-genome events are shipped to the client at once (hundreds of KB is fine) and are
  cached per run id in `api.ts`, so the map and the alignment viewer share one fetch.
- Only reference bases for the visible window are fetched on demand via a `refseq`
  endpoint (seqid + start + end, ~5 kb cap) for letters mode.

### Types (crates/types/src/tables.rs)

```rust
pub struct AlignmentEvents {
    pub snps: Vec<(u64, u8, u8)>,      // (ref_pos_1based, ref_base, query_base)
    pub dels: Vec<(u64, u64)>,          // (ref_pos, length)
    pub ins: Vec<(u64, String)>,         // (ref_pos, inserted bases)
}
pub struct AlignmentQuery {
    pub query_id: i64,
    pub query_name: String,
    pub blocks: Vec<WgaBlock>,           // reuse existing block type
    pub events: BTreeMap<String, AlignmentEvents>, // key = reference contig seqid
}
pub struct AlignmentData {
    pub reference: Vec<(String, u64)>,   // (seqid, length) in reference order
    pub queries: Vec<AlignmentQuery>,
}
```

### Endpoints

- `GET /runs/{id}/alignment` returns `AlignmentDataColumnar`: per query and seqid the
  events are parallel arrays (`snp_pos/snp_ref/snp_qry`, `del_pos/del_len`,
  `ins_pos/ins_seq`) instead of one object per event — a divergent query carries
  >100k SNPs and the object form made ~60 MB responses that froze the browser's
  JSON parser. `api.ts` hydrates the arrays into the per-event objects the
  consumers index, once per fetch inside the cached promise. New runs read the
  persisted variants.json directly. On a cache miss (old runs) the events are
  computed per query in parallel (cpu_slots-bounded, reference fasta parsed once,
  per-query locks to avoid duplicate work) and written atomically (temp file +
  rename). variants.json on disk keeps the object form; only the HTTP response is
  columnar.
- `GET /runs/{id}/refseq?seqid=..&start=..&end=..` returns uppercase reference bases for
  the window (clamped to ~5 kb).

### Engine

- New module `crates/engine/src/variants.rs`: walks a nucmer delta file (same walking
  pattern as `msa::reconstruct` in `pipeline.rs` / `coverage.rs`), emits
  `BTreeMap<seqid, AlignmentEvents>` for one query. Overlapping nucmer blocks are
  deduplicated (first block covering a position wins).
- Reverse-strand blocks: query bases are revcomp'd before comparison so events are in
  reference orientation.

## Frontend structure

- Mode selector in `ProjectPage.tsx` (URL param `gv=strain|align`), locus preserved
  across mode switches.
- `GenomeView.tsx` evolves into StrainMap:
  - remove identity track rows and `identityColor`
  - genes on a single line, strand = arrow direction
  - label inside rectangle = symbol, fallback locus_tag
  - ruler for positions
  - variant markers on zoom (< ~20 kb window)
- New `AlignmentView.tsx` porting Oligool's zoom/pan:
  - viewFraction 0.005..1, spacer div + sticky canvas, Ctrl/Cmd+wheel zoom,
    zoom rect, DPR handling, virtualization, letters threshold ~85 bases
    with hysteresis 100.

## Skipped in v1

GC track, minimap, clean-region autofind, primers, NCBI modal, copy-selection.

## Implementation order

1. `variants.rs` + types + engine tests
2. Endpoints + API test
3. Mode selector + StrainMap rework
4. Gene labels + variant markers
5. AlignmentView bars + zoom
6. Letters mode
7. Cross-mode sync / legend / dark mode
8. fmt / clippy / test / build, deploy
