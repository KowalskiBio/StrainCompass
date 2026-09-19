# straincompass

Bacterial genome comparison workbench: compare one or more query genomes
(draft or complete assemblies) against an annotated reference, inspect the
results in a large interactive table, and explore an interactive whole
genome alignment (WGA) view.

## Vision

The app wraps and extends the logic of
`plato/radka/compare_genome_vs_reference.R`: reference genes scored as
present / PARTIAL / ABSENT by alignment coverage, unaligned gap reports,
and an optional strict BLAST recheck of a gene panel. Viruses will be
supported later; bacteria come first.

## Deployment strategy

1. **Phase 1 (now): web app on a VM.** Served from a Proxmox VM. Users
   access it through a browser. Jobs (nucmer, BLAST) run on the VM.
2. **Phase 2 (later): standalone app.** Same core, packaged as a desktop
   application for offline use.

The architecture must keep the analysis core decoupled from the UI so the
same engine serves both deployment modes.

## Core concepts

- **Project**: a unit of work, created by the user.
- **Organism type per project**: bacteria first (viruses later).
- **Reference**: user supplies a bacterial reference genome (FASTA) and
  its annotation (GFF). Optionally fetched from NCBI by accession
  (e.g. GCF_000196035.1); a user provided NCBI API key (Settings page)
  speeds these downloads up.
- **Inputs (queries)**: one or multiple query FASTA files (user strains,
  draft or complete).
- **Engine**: sanitize FASTA headers, run nucmer (query vs reference) and
  dnadiff, compute per gene alignment coverage with merged intervals,
  build the unaligned gap table, and optionally re BLAST a gene panel
  with strict presence thresholds (>= 90% coverage AND >= 90% identity).

## Outputs

### 1. Table view

A window with a huge, virtualized table, like `View(table)` in R but in
the browser. It shows the equivalent of the `.tsv` files the R script
produces:

- `genes_coverage.tsv`: every reference gene with locus_tag, symbol,
  biotype, seqid, start, end, length, cov_bp, cov_pct, call
  (present / PARTIAL / ABSENT)
- `unaligned_gaps.tsv`: unaligned regions (>= min_gap bp) with genes
  inside them
- `panel_recheck.tsv`: strict BLAST calls for a gene panel (optional)
- `dnadiff.report`: overall alignment statistics

Features: column visibility chosen by the user, sorting, filtering,
searching, multi query comparison (one column set per query or a joined
presence/absence matrix across all queries in the project).

### Export / download to the user's computer

Every result can leave the app as a file, with a big always visible
"Export" button:

- On each table (genes coverage, unaligned gaps, panel recheck,
  presence/absence matrix): "Export" with a format choice, TSV or CSV.
  Exports respect the user's current view: active column selection,
  sorting and filters are what gets written, so the downloaded file
  matches what the user sees (with a "current view" vs "all columns"
  choice, default current view).
- On the gene alignment (MSA) viewer: export the shown alignment as
  FASTA or clustal.
- On the run "Files" panel: each artifact (including the dnadiff report)
  downloadable as is.
- Downloads stream straight from the browser to the user's computer
  (standard browser download, no server side zip required; a "download
  all tables" button that bundles them in one zip is a nice extra).
  Served via `GET /runs/{id}/export/{table}?format=tsv|csv`, marked
  `Content-Disposition: attachment` so the file lands in the user's
  Downloads folder with a clear name like
  `LM259_vs_reference_genes_coverage.tsv`.

### 2. Interactive genome view

A window with an interactive genome built from the input vs reference
WGA:

- linear and/or circular reference genome rendering
- aligned blocks drawn as tracks per query, color coded by identity
- unaligned gaps highlighted, genes inside gaps annotated
- gene tracks from the reference GFF (protein coding, tRNA, rRNA)
- zoom, pan, click a gene or block for details, tooltip with locus tag,
  symbol, coverage, coordinates

### 3. Gene alignment (MSA) viewer

Drilling from "this gene is PARTIAL" down to the actual bases:

- **Hover** on any gene in the table or genome view shows a quick
  preview: coverage %, best block identity, mismatch and indel counts
  for each query.
- **Click** (or a "Show alignment" button) opens an MSA viewer dialog:
  a pairwise alignment of the reference gene versus the aligned query
  region(s), one row per query, rendered as a scrollable sequence
  track. Mismatches, insertions and deletions are highlighted in
  place; a coordinate ruler shows gene position; premature stop codons
  in the query are flagged. When a gene is split across several
  alignment blocks, blocks are shown in order with the unaligned
  stretches marked as gaps.
- **Data source**: computed from the run's cached delta/alignment
  artifacts (MUMmer's delta encodes the full alignment including
  indels; show-snps provides the variant list). No re-alignment needed
  for reasonable gene sizes; the endpoint renders on demand.
- **Export**: the gene alignment as FASTA / clustal for reporting.

## Architecture

Strict backend / frontend split. The frontend never runs any analysis;
it only collects parameters and renders results. The backend is Rust
end to end: one workspace, zero Python. The engine lives in the backend
as a self contained crate so it can later be embedded in the desktop
build unchanged (Tauri is Rust too, so phase 2 links the very same
crate into the app).

```
browser SPA  <── JSON over HTTP ──>  API server  ──>  engine
  (React)          (REST)          (Rust, axum)  (pipeline)
                                       │
                                       └──> job runner ──> nucmer / show-coords /
                                             (tokio tasks)    dnadiff / blastn
```

Honest note on speed: the wall clock of a run is dominated by the
external tools (nucmer, BLAST), which stay external regardless of
language. What Rust buys here is everything around them: FASTA/GFF/
delta parsing and coverage math in milliseconds instead of seconds, a
snappy API with low memory, a single static binary per service, no
runtime to install on the VM, and one language shared with the desktop
packaging. A bacterial genome (3k genes, ~3 Mb) parses and scores in
well under a second.

### Making the tools faster

- **Parallelize across queries, not within a run.** A single nucmer
  invocation is effectively single threaded: extra cores do NOT make one
  comparison faster. But a project with N queries means N independent
  nucmer runs, which the job runner executes in parallel on all cores.
  This is where "more cores" actually pays off, and it is the main
  scaling path in the app.
- **RAM is a non issue for bacteria.** The suffix tree scales with
  genome size; a 3-10 Mb bacterial genome needs a few hundred MB.
  Allocate cores, not RAM.
- **nucmer knobs**: default settings are already fast at bacterial
  scale (seconds to a minute per comparison). `--mum` (default, unique
  matches) is faster than `--maxmatch`; a larger `-l` cuts matching
  work. Do not tune for speed, tune for correctness.
- **Optional alternative aligner**: minimap2 (`-t` for real multi
  threading, PAF output) can serve as a second alignment backend for
  the big table if nucmer ever becomes the bottleneck; the delta
  pipeline stays the reference/default. Evaluate during engine work.
- **NCBI API key (user provided).** Reference genomes and annotations
  fetched by accession go through NCBI, which throttles anonymous
  requests hard (about 3/s). A user supplied NCBI API key raises the
  limit (about 10/s), making reference download and multi accession
  fetches much faster. The key is entered in Settings, stored server
  side (never logged, never returned in full by the API), and attached
  to every NCBI request the engine makes. No key: everything still
  works, just slower.

```
PUT    /settings/ncbi_api_key       store key (masked in responses)
DELETE /settings/ncbi_api_key       remove key
GET    /settings                    show masked key + status
```

### Backend (Rust workspace)

- **Crates**: `straincompass-engine` (pipeline, no web code), `straincompass-api`
  (HTTP server), `straincompass-types` (shared DTOs + parameter model). The
  engine has no dependency on axum or the DB, so the Tauri build uses it
  directly.
- **API server**: axum + tokio. Owns projects, files, runs, users.
  OpenAPI schema via utoipa; the frontend client is generated from it.
- **Engine**: a Rust port of the R script logic. Steps: FASTA header
  sanitization, nucmer, delta header rewrite, show-coords parsing,
  per gene coverage with merged intervals, gap computation, optional
  BLAST panel recheck. Rust makes the parsing and interval math
  (nom + custom iterators, sorted vec merging) effectively free
  compared to tool runtimes. Each step is a pure function of (input
  files, parameter set) where possible, so parameters can be changed
  without redoing expensive alignment work.
- **Job runner**: analysis runs are async jobs on tokio tasks with a
  SQLite backed queue (escalate to a dedicated queue later if needed).
  Long tool invocations stream stdout line by line into run logs.
  Status polled via `GET /runs/{id}` or pushed over WebSocket; logs
  streamed to the UI.
- **Storage**: SQLite via sqlx/rusqlite for metadata (later Postgres if
  multi user load demands it); one directory per project on disk for
  FASTA/GFF/delta/TSV artifacts; `runs` table keeps the exact parameter
  JSON for every run so any result is reproducible.
- **Web frontend serving**: in phase 1 the same axum binary serves the
  compiled SPA from its static dir, so the VM stack needs just one app
  container plus nginx.

### Engine parameters (tweakable from the frontend)

Every threshold the R script hardcodes becomes a named parameter with a
default, a range, and a UI control. The frontend sends a parameter JSON
with each run request; the backend validates it (serde + garde
validation attributes, rejects out of range values with a field level
error response) and stores it with the run.

| Parameter        | Default | Range    | Controls                                    |
|------------------|---------|----------|---------------------------------------------|
| `min_gap`        | 200     | 0-100000 | minimum unaligned region reported (bp)      |
| `present_cov`    | 95      | 0-100    | % of gene aligned to call it present       |
| `partial_cov`    | 1       | 0-100    | >0% and <present_cov = PARTIAL (implicit)  |
| `blast_cov`      | 90      | 0-100    | % query coverage for panel gene present    |
| `blast_pid`      | 90      | 0-100    | % identity for panel gene present          |
| `blast_evalue`   | 1e-10   | 1e-50-10 | BLAST e value cutoff                       |
| `nucmer_minmatch`| (nucmer default) | int | nucmer `-l` minimal match length |
| `nucmer_breaklen`| (nucmer default) | int | nucmer `-b` breakpoint distance |
| `dnadiff`        | on/off  | bool     | run overall alignment report or skip        |

Cheap re runs: `min_gap`, `present_cov` and the BLAST thresholds only
affect post processing of a finished alignment, so changing them
recomputes tables in seconds without re-running nucmer. Only
`nucmer_*` parameters and input changes trigger a real re-alignment.
The run model records which layer a parameter belongs to
(postprocess vs align) and the engine reuses cached artifacts
accordingly.

### API sketch

```
POST   /projects                      create project (name, organism=bacteria)
GET    /projects                      list projects
GET    /projects/{id}                 project detail
POST   /projects/{id}/reference       upload FASTA + GFF (or NCBI accession)
POST   /projects/{id}/queries         upload one or many query FASTAs
POST   /projects/{id}/panel           upload gene panel FASTA (optional)
POST   /projects/{id}/runs            start run: {query_ids, params{...}}
GET    /runs/{id}                     status, progress, logs
GET    /runs/{id}/genes_coverage      table data (paginated, sortable)
GET    /runs/{id}/unaligned_gaps      table data
GET    /runs/{id}/panel_recheck      table data (if panel given)
GET    /runs/{id}/wga                 alignment blocks + gaps + genes for viewer
GET    /runs/{id}/gene/{locus}       gene preview stats (hover) + full
                                      alignment rows for the MSA viewer
GET    /runs/{id}/export/{table}      TSV/CSV download
GET    /runs/{id}/params              parameter set used (reproducibility)
```

### Frontend

- **Stack**: React + TypeScript + Vite, Tailwind (or similar).
- **Target users**: mostly NOT IT experienced people (lab scientists,
  students). The UI must be usable by someone who has never seen a
  command line and does not know what nucmer or a delta file is.
- **Pages**:
  - project list / create project
  - project detail with tabs: Inputs (reference, queries, panel),
    Runs, Table view, Genome view
- **Parameter panel**: a form (sliders, number inputs, toggles) bound to
  the engine parameter schema above. Presets ("strict", "loose",
  "default") plus a diff against defaults. Submitting starts a run;
  re-submitting with only postprocess thresholds changed shows a "fast
  recompute" hint.
- **Table view**: virtualized grid (AG Grid or TanStack Table) fed by
  the paginated run endpoints; column visibility picker, global search,
  per column filters, sort. Handles 10k+ rows smoothly.
- **Genome view**: interactive WGA viewer. Reference coordinates on a
  horizontal axis (circular mode later), one alignment track per query,
  blocks colored by % identity, unaligned gaps shaded, gene track from
  the GFF. Zoom/pan, click for details (tooltip popover with locus tag,
  symbol, coverage, coordinates), deep links from table rows to the
  corresponding genome position and back.
- **Gene alignment (MSA) viewer**: hover preview plus a dialog MSA view
  for any PARTIAL gene (mismatches, indels, premature stops), fed by
  `GET /runs/{id}/gene/{locus}` (see Outputs, section 3).
- **Live updates**: run progress, tool logs and errors surface in a run
  drawer (WebSocket or polling).

### UI / UX principles (non IT audience)

The users are biologists, not developers. Every screen follows:

- **Clean and minimalistic**: one primary action per screen. No menus
  inside menus, no jargon in labels. Say "Compare genomes", not "Run
  nucmer pipeline". Technical terms live in tooltips / an info icon,
  never in button text.
- **Big and greatly visible**: large buttons, large fonts (16 px base,
  headings much bigger), generous spacing, high contrast colors.
  Comfortable on a laptop projector in a lab meeting. Touch friendly
  hit areas, nothing smaller than ~44 px.
- **Forgiving inputs**: drag and drop file upload with clear accepted
  formats ("FASTA file, e.g. LM259.fasta"), inline validation with
  plain language errors ("This file is not a FASTA file"), no way to
  reach a raw stack trace. Destructive actions ask for confirmation.
- **Guided flow**: a project is created through a small numbered wizard
  (1 reference, 2 queries, 3 optional panel, 4 run), shown as a nice
  dialog popup (see "Everything through the UI"). Defaults are
  always sensible: a user who never opens a parameter panel still gets
  a valid, correct result. Parameters are hidden behind an "Advanced
  settings" section, collapsed by default, with plain language
  descriptions next to each control and safe ranges enforced by the
  backend.
- **Obvious state**: running jobs show a big progress indicator and
  elapsed time, not a spinner alone. Failures say what to do next
  ("Reference file could not be read. Try uploading it again or
  contact support"), not an error code.
- **Forgiving navigation**: table and genome views are always reachable
  by big tabs; a selected row, filter or zoom state survives
  navigation and reloads; nothing is lost by clicking the wrong
  button.
- **Accessibility and low reading load**: icon + text on every button,
  color never the only signal (color blind safe palette for
  present / PARTIAL / ABSENT), keyboard support for the grid.
- **Test with real users**: before each release, watch one target user
  complete the full flow (create project, upload, run, find an absent
  gene) without help. If they hesitate, the UI is wrong, not them.

### Everything through the UI (no filesystem exposure)

The app is completely UI based. Users never see, type or browse a
server path, and never open a terminal. Concretely:

- **The only things a user ever provides**: their NCBI API key
  (optional, Settings), the reference genome (FASTA), its annotation
  (GFF), one or more query FASTAs, and optionally a gene panel FASTA.
  Exactly the inputs the R script demanded, nothing more. Tools
  (MUMmer, BLAST+) are pre installed in the app container: installing
  them is the operator's one time job, invisible to users.
- **Everything else is the app's business**: working directories,
  sanitized FASTAs, delta files, caches, intermediate files. They exist
  on disk but the user never needs to know where. No path inputs
  anywhere in the UI.
- **All outputs are reachable in the UI**: every table, the dnadiff
  report, run logs, and per run parameter sets are viewable in the app
  and downloadable as files (TSV/CSV exports, report as text). A run
  page shows its full artifact list; a "Files" panel lists what the run
  produced with friendly names ("Genes coverage table",
  "Alignment report"), each with View and Download buttons.
- **Projects are managed in the UI**: create, rename, delete (with
  confirmation), and revisit old runs with their exact parameters and
  results. Storage quotas, if any, are shown as a plain number
  ("This project uses 340 MB").
- **Guided input flow**: adding the reference, queries and panel is a
  friendly dialog popup (wizard): step 1 pick reference (upload files
  or fetch from NCBI by accession with a big search box), step 2 add
  query FASTAs (drag and drop, multiple at once), step 3 optional gene
  panel, step 4 review and "Compare genomes". Each step validates
  immediately and explains in plain language what is missing. The
  dialog pops up automatically when a project is created and can be
  reopened any time from a big "Add / change inputs" button.

Supporting endpoints (run artifacts surfaced through the API, not the
filesystem):

```
GET    /runs/{id}/files              list artifacts with friendly names
GET    /runs/{id}/files/{name}       view or download one artifact
GET    /projects/{id}/usage          storage used by the project
```

### Deployment (phase 1, VM)

- One docker compose stack on the Proxmox VM: `app` (single static Rust
  binary: API + job runner + SPA static serving), `db` (SQLite volume;
  Postgres only if multi user load demands it), nginx as reverse proxy
  with TLS. Because everything is one static binary, the container is
  distroless and starts in milliseconds.
- MUMmer and BLAST+ live in the app image (apt/bioconda), no
  auto-download at runtime.
- Phase 2 desktop build (Tauri): the engine crate links directly into
  the app process, no separate backend, storage moves to a local
  directory. The REST layer becomes an in process command interface;
  the frontend is the same codebase.

## Roadmap

1. [ ] Repo scaffolding (Rust workspace: engine, api, types + frontend),
       CI, tool container (MUMmer + BLAST+)
2. [ ] Engine crate: port the R logic (sanitize, nucmer, coverage, gaps,
       panel recheck) with a validated parameter model and layered
       caching (align vs postprocess artifacts)
3. [ ] API server (axum): projects, uploads, runs, parameter validation,
       job runner, result endpoints
4. [ ] Frontend shell: project pages, guided input wizard dialog (big
       minimal UI per the UX principles), run drawer with logs, run
       artifact Files panel (view + download, no filesystem exposure)
5. [ ] Table view: virtualized grid, column picker, filters, exports
       (TSV/CSV)
6. [ ] Multi query support: run several queries, presence/absence matrix
7. [ ] Interactive WGA view: reference genome rendering, alignment
       tracks, gap highlighting, gene tooltips, deep links from table
8. [ ] Gene alignment (MSA) viewer: hover preview, per gene pairwise
       alignment with mismatches/indels, export FASTA/clustal
9. [ ] Usability test: watch target (non IT) users complete the full
       flow, fix friction
10. [ ] VM deployment on Proxmox: per deployment.md (native binary +
        systemd + existing nginx, port 127.0.0.1:8010, backup before
        every push, never touch foreign prod services)
11. [ ] Hardening: result caching, quotas, users
12. [ ] Phase 2: standalone desktop packaging (Tauri, engine in
        process)
13. [ ] Virus support (organism type: virus)

## Known caveats carried over from the R script

- Collapsed multicopy repeats (rRNA operons) can show individual copies
  as unaligned: assembly artifacts, not true absences. Absent genes are
  reported split by gene_biotype so tRNA/rRNA artifacts can be spotted.
- Genes at contig boundaries can appear truncated or absent: confirm
  with the BLAST recheck panel.
- A fully aligned gene can still carry a premature stop codon; sequence
  presence does not prove a functional ORF.
- Gained regions are query stretches with *no alignment to the reference*,
  which is weaker than *absent from the reference*. nucmer anchors on
  matches unique to the reference (`--mumreference`), so a query copy of a
  multicopy reference family (rRNA operon, IS element, transposase) may
  have nothing unique to seed from and be reported as gained even though
  the reference carries it several times over. This is the query-side
  mirror of the first caveat above. The rows flag what they can -
  `at_contig_end`, `flanks_disagree`, and the complete-ORF count rather
  than the raw one - and the Gained tab's region card offers the
  confirming search on demand: the region's sequence searched back
  against the reference (`GET /runs/{id}/gained/verify`) three ways -
  nucleotide, translated, and the longest exact match - with nothing
  found reported as consistent with a true gain.
- A draft query assembly contributes an unaligned tip at every contig end,
  so `min_gained` (500 bp by default) and the `at_contig_end` flag are load
  bearing on fragmented input.

## Gained regions (added 2026-09-17)

Each query now also gets the mirror image of the unaligned-gaps table:
stretches of the *query* genome that no alignment covers, placed on the
reference by the alignment blocks flanking them, with the genes inside
them predicted by prodigal.

- `crates/engine/src/gained.rs` complements `DeltaFile::aligned_qry_intervals()`
  against the query contig lengths, then anchors each region as `Between`
  (both flanks agree), `Flank` (one usable flank, or a fallback from
  disagreeing ones) or `Unanchored` (the contig has no alignment at all -
  a plasmid).
- **prodigal is optional.** `ToolPaths::prodigal` is an `Option`, because
  `discover()` is all-or-nothing and no deployment made before this
  feature has the binary. Without it the regions are still computed and
  the run still succeeds; the ORF counts are `None`, which the tables and
  the map render as "not available" and never as zero.
- Surfaced as a per-query **Gained** tab (hidden on runs that predate the
  feature, via `RunDto.has_gained`) and a **gained** toggle in the strain
  map's legend, which draws a teal mark at each anchor. The parameters are
  `min_gained` and `gained_orfs`, both postprocess-layer, so changing them
  re-runs cheaply and does not invalidate a cached delta.
- Gained rows are per query in query coordinates, so they deliberately do
  not enter the presence/absence matrix, whose call arrays are positionally
  parallel to the reference gene list. "Which gains do strains A and B
  share?" is therefore not answerable yet; that needs cross-query
  clustering of the regions by sequence identity.
- **Reference back-check** (added 2026-09-18): the region card in the
  Gained tab can search the region's sequence against the whole
  reference genome on demand (`crates/engine/src/blast.rs::
  gained_verify`, `GET /runs/{id}/gained/verify`). Three searches are
  reported, because each covers a blind spot of the others:
  - A **nucleotide search** (blastn, `-task blastn` 11-mer seeding,
    `-dust no`) run once at a loose ceiling (E <= 10) and split into a
    strong tier (E <= 1e-5) the verdicts are built on, and a weak tier
    shown collapsed rather than silently dropped. It is deliberately
    *more* sensitive than the pipeline's aligner, because it exists to
    catch what the aligner's unique-reference anchors miss.
  - A **translated search** (tblastx, E <= 1e-5, `-seg no`) that runs
    only when the strong nucleotide tier is empty and the region is at
    most 50 kb: the second opinion on "found nothing", which finds
    divergent coding homologs below ~70% nucleotide identity. Over the
    cap it is skipped with a note instead of holding a cpu slot for
    minutes.
  - The **longest exact match** (`crates/engine/src/longest_match.rs`,
    a suffix automaton over the region scanning the reference): the
    longest run of bases the two share verbatim, cutoff-free, reported
    with example coordinates on both sides. Unrelated DNA of these
    sizes shares ~log4(n*m) bases by chance, so the number is read
    against that expectation, not against a threshold.
  Results are cached client-side per region; the server keeps no cache,
  so the check always runs against the reference as it stands now. The
  standing limit, stated in the UI copy: no sequence search proves
  absence - "not found, even translated" is the closest available
  evidence of a true gain. Every hit of every tier is labeled with the
  reference genes its interval overlaps, because a row of coordinates
  answers nothing until it says what it hits.
- **Gene identification** (added 2026-09-18): queries are unannotated
  draft assemblies, so their genes arrive from prodigal with coordinates
  and nothing else - the gained table's LM4B_RS... columns are the
  *reference* genes flanking the insertion point, not the gained genes,
  and nothing in the system could say what the gained genes are. The
  region card can now name them (`crates/engine/src/blast.rs::
  gained_identify`, `GET /runs/{id}/gained/identify`): each predicted
  ORF is searched, as translated DNA (blastx, E <= 1e-5), against a
  protein database built on the fly from the reference's own annotated
  CDSs (translated from ref.fa by the GFF, standard genetic code; genes
  whose translation is not a clean protein are skipped). This names the
  ORFs the reference does carry homologs of - interrupted copies,
  repeat-family members, the "probably not a true gain" cases. The ORFs
  it leaves unnamed are exactly the novel ones: for those the response
  carries their sequences, and the client links each to NCBI BLAST
  (blastx vs nr) and the whole region likewise, since no local database
  can name a gene the reference has never seen.
- **Names in the table itself** (added 2026-09-18, hours later): the
  naming pass now runs with the run, not on demand - after prodigal
  predicts the genes of all of a query's gained regions, one blastx run
  over all of them at once (`gained::gained_regions` calls
  `blast::name_gained_orfs`) fills each ORF's `best` and each row's
  `gene_names`. The gained table carries a "Genes inside (named)"
  column ("all novel" when predicted genes matched no reference
  protein), the search box matches names, and the TSV/CSV export has a
  `gene_names` column plus each ORF's name in the `orfs` column
  (`1101..1700(+) =LM4B_RS13330 creatininase...`). Like gene prediction,
  naming is decoration on an already-complete result: a failure leaves
  genes unnamed rather than failing the comparison. Runs computed before
  the pass have empty `gene_names` (serde default) and name nothing
  until re-run; the on-demand identification endpoint covers them.

## NCBI cross references (added 2026-09-15)

Reference genes now carry `protein_id`, `product` and `gene_id`, so the
results table has a "Protein (NCBI)" column and the alignment dialog links
the gene name straight to NCBI.

The catch the implementation has to handle: RefSeq GFFs put `protein_id`
and `product` on the **CDS** child, not on the `gene` feature. Both the R
script and the engine read gene-level features only, so neither could see
an accession. Both now join the CDS attributes back by the shared
`locus_tag` (first CDS wins for genes with several children).

Link priority, identical in `frontend/src/ncbi.ts` and the R script so the
TSV and the web UI point at the same record:

1. `Dbxref=GeneID:<n>` -> `/gene/<n>` (best, but older RefSeq GFFs lack it;
   NC_003212.1 has zero)
2. `protein_id=WP_...` -> `/protein/<acc>` (3089 of 3219 genes in EGD-e)
3. sequence accession + coordinates -> `/nuccore/<seqid>?from=&to=`
   (always resolves, and is what pseudogenes and RNAs fall back to)

Note the table column is written into each query's stored results at run
time, so runs made before this change show "-" until re-run. The alignment
dialog re-parses the GFF live and therefore works on old runs immediately.
