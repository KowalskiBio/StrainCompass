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
  (e.g. GCF_000196035.1).
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

## Tech stack (proposal, to be decided in repo docs)

- Backend: Python (FastAPI) wrapping the pipeline; MUMmer (nucmer,
  show-coords, dnadiff) and NCBI BLAST+ as external tools, managed in a
  container on the VM
- Frontend: web UI with a virtualized data grid (e.g. AG Grid / TanStack
  Table) for the table view, and a genome browser library (e.g.
  PGDV.js / custom SVG canvas) for the WGA view
- Storage: project directory per project on the VM disk
- Phase 2: package backend + frontend core as desktop app (e.g. Tauri or
  Electron) with local tool bundling, reusing the same engine

## Roadmap

1. [ ] Repo scaffolding, CI, tool container (MUMmer + BLAST+)
2. [ ] Core engine: port the R logic (sanitize, nucmer, coverage, gaps,
       panel recheck) to a service API
3. [ ] Project + organism model: create project, attach reference
       (upload FASTA/GFF or NCBI accession), upload queries
4. [ ] Table view: virtualized grid, column picker, filters, exports
       (TSV/CSV)
5. [ ] Multi query support: run several queries, presence/absence matrix
6. [ ] Interactive WGA view: reference genome rendering, alignment
       tracks, gap highlighting, gene tooltips
7. [ ] VM deployment on Proxmox (reverse proxy, TLS, auth)
8. [ ] Hardening: job queue, result caching, quotas
9. [ ] Phase 2: standalone desktop packaging
10. [ ] Virus support (organism type: virus)

## Known caveats carried over from the R script

- Collapsed multicopy repeats (rRNA operons) can show individual copies
  as unaligned: assembly artifacts, not true absences. Absent genes are
  reported split by gene_biotype so tRNA/rRNA artifacts can be spotted.
- Genes at contig boundaries can appear truncated or absent: confirm
  with the BLAST recheck panel.
- A fully aligned gene can still carry a premature stop codon; sequence
  presence does not prove a functional ORF.
