# bactiment

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

### 2. Interactive genome view

A window with an interactive genome built from the input vs reference
WGA:

- linear and/or circular reference genome rendering
- aligned blocks drawn as tracks per query, color coded by identity
- unaligned gaps highlighted, genes inside gaps annotated
- gene tracks from the reference GFF (protein coding, tRNA, rRNA)
- zoom, pan, click a gene or block for details, tooltip with locus tag,
  symbol, coverage, coordinates

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

- **Crates**: `bactiment-engine` (pipeline, no web code), `bactiment-api`
  (HTTP server), `bactiment-types` (shared DTOs + parameter model). The
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
8. [ ] Usability test: watch target (non IT) users complete the full
       flow, fix friction
9. [ ] VM deployment on Proxmox (docker compose, reverse proxy, TLS,
       auth)
10. [ ] Hardening: result caching, quotas, users
11. [ ] Phase 2: standalone desktop packaging (Tauri, engine in
        process)
12. [ ] Virus support (organism type: virus)

## Known caveats carried over from the R script

- Collapsed multicopy repeats (rRNA operons) can show individual copies
  as unaligned: assembly artifacts, not true absences. Absent genes are
  reported split by gene_biotype so tRNA/rRNA artifacts can be spotted.
- Genes at contig boundaries can appear truncated or absent: confirm
  with the BLAST recheck panel.
- A fully aligned gene can still carry a premature stop codon; sequence
  presence does not prove a functional ORF.
