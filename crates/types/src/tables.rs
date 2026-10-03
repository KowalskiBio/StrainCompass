//! Row types for the result tables and the viewer endpoints.

use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, Default)]
#[serde(rename_all = "UPPERCASE")]
pub enum Call {
    #[default]
    Present,
    Partial,
    Absent,
}

impl Call {
    pub fn as_str(&self) -> &'static str {
        match self {
            Call::Present => "PRESENT",
            Call::Partial => "PARTIAL",
            Call::Absent => "ABSENT",
        }
    }
}

/// One row of the genes coverage table (per query).
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct GeneCoverageRow {
    pub locus_tag: String,
    pub symbol: String,
    /// RefSeq protein accession, e.g. "WP_003759425.1"; empty when the gene
    /// has no coding child in the annotation.
    #[serde(default)]
    pub protein_id: String,
    pub biotype: String,
    pub seqid: String,
    pub start: u64,
    pub end: u64,
    pub length: u64,
    pub cov_bp: u64,
    pub cov_pct: f64,
    pub call: Call,
    /// Best block identity over this gene (0-100), 0 when nothing aligned.
    pub best_identity: f64,
    pub mismatches: u64,
    pub indels: u64,
    /// Where the gene lies in the query assembly: the span its aligned
    /// part maps to, as "contig:start-end(+|-)", several joined with "; "
    /// (most aligned bases first) when it is split or duplicated. Empty
    /// when nothing aligned.
    #[serde(default)]
    pub qry_loci: String,
}

/// One row of the unaligned gaps table.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct GapRow {
    pub seqid: String,
    pub start: u64,
    pub end: u64,
    pub length: u64,
    pub n_genes: usize,
    /// Locus tags of genes that overlap the gap.
    pub genes: Vec<String>,
}

/// Where a gained region sits relative to the reference. A gained region
/// has no reference coordinates of its own, so it is placed by the
/// alignment blocks that flank it on the query contig.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, Default)]
#[serde(rename_all = "snake_case")]
pub enum GainedAnchor {
    /// Both flanking blocks land close together on the same reference
    /// sequence: the region sits between `anchor_start` and `anchor_end`.
    Between,
    /// Only one flank is usable - the region is at a query contig end, or
    /// the two flanks pointed at different places and we fell back to the
    /// left one. The position is that flank's junction, half a placement.
    Flank,
    /// The query contig has no alignment to the reference at all, so the
    /// region cannot be placed. A whole extra replicon (a plasmid, a
    /// phage) lands here.
    #[default]
    Unanchored,
}

impl GainedAnchor {
    pub fn as_str(&self) -> &'static str {
        match self {
            GainedAnchor::Between => "BETWEEN",
            GainedAnchor::Flank => "FLANK",
            GainedAnchor::Unanchored => "UNANCHORED",
        }
    }
}

/// One gene predicted inside a gained region.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct GainedOrf {
    /// 1-based inclusive, in QUERY CONTIG coordinates (not relative to the
    /// region), so the ORF can be pulled straight out of the query fasta.
    pub start: u64,
    pub end: u64,
    pub strand: i8,
    /// The ORF runs off an edge of the region, so it is probably truncated
    /// rather than a whole gene.
    pub partial: bool,
    /// The gene finder's confidence in the call (0-100).
    pub confidence: f64,
    /// Best match among the reference's own proteins, when the run's
    /// naming pass found one. Absent in results computed before the
    /// pass existed, or when the gene is novel to the reference.
    #[serde(default)]
    pub best: Option<OrfMatch>,
    /// Best match in NCBI's nr, from the on-demand naming of novel
    /// genes (persisted per query, merged wherever rows are read).
    /// Absent when not yet searched or nothing significant was found.
    #[serde(default)]
    pub ncbi: Option<OrfMatch>,
}

/// One stretch of a query genome with no alignment to the reference, and
/// the genes predicted inside it.
///
/// "No alignment" is weaker than "not in the reference": nucmer anchors on
/// matches unique to the reference, so a query copy of a gene the reference
/// carries several times over may have nothing unique to seed from and land
/// here anyway. The wording everywhere downstream says "has no alignment to
/// the reference", never "is new".
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct GainedRow {
    /// Query contig the region sits on (sanitized id, as in query.fa).
    pub qry_seqid: String,
    /// 1-based inclusive coordinates on the query contig.
    pub start: u64,
    pub end: u64,
    pub length: u64,
    /// G+C percent over the unambiguous bases of the region. Far from the
    /// genome's own GC is the classic hint of horizontally acquired DNA.
    pub gc_pct: f64,
    /// True when the region touches either end of its query contig, where a
    /// draft assembly's breaks routinely look like gained sequence.
    pub at_contig_end: bool,
    pub anchor: GainedAnchor,
    /// Reference sequence the anchor is on; empty when unanchored.
    pub anchor_seqid: String,
    /// The reference interval the region is placed at (1-based inclusive).
    /// Both ends are the junction position for `Flank`; both are 0 when
    /// unanchored.
    pub anchor_start: u64,
    pub anchor_end: u64,
    /// Set when both flanks exist but tell different stories (different
    /// reference sequences, far apart, or opposite orientations): the anchor
    /// fell back to the left flank and should be read with care.
    pub flanks_disagree: bool,
    /// Locus tag of the reference gene at the left flank, and at the right.
    /// Empty when there is no gene there.
    pub left_gene: String,
    pub right_gene: String,
    /// Genes predicted inside the region. `None` means gene prediction did
    /// not run - which is not the same as having looked and found none.
    pub n_orfs: Option<u32>,
    /// Of those, the ones with both a start and a stop inside the region.
    /// This is the number worth quoting: an ORF running off an edge is
    /// usually a gene the insertion interrupted, or one a contig break cut.
    pub n_orfs_complete: Option<u32>,
    /// The predicted genes themselves.
    pub orfs: Vec<GainedOrf>,
    /// The names of the predicted genes the reference's own proteins
    /// could identify (symbol or product, as annotated), best-match
    /// first within each gene. Genes absent from this list are novel to
    /// the reference - not found by a translated search against every
    /// protein it encodes. Empty in results computed before the naming
    /// pass existed, and in runs whose gene prediction was switched off.
    #[serde(default)]
    pub gene_names: Vec<String>,
    /// Whether the naming pass ran for this row. An empty `gene_names`
    /// means "novel" only when this is true; false marks results from
    /// before the pass existed (or a naming failure), where empty means
    /// unknown, not unmatched. The same distinction `n_orfs: null`
    /// already makes for prediction.
    #[serde(default)]
    pub named: bool,
}

/// Whether the genes inside the gained regions could be predicted.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(tag = "state", content = "reason", rename_all = "snake_case")]
pub enum GainedOrfStatus {
    /// Predicted successfully.
    Predicted,
    /// Not attempted or not possible; carries a plain language reason.
    Unavailable(String),
    /// The run predates the feature, so nothing is known either way.
    #[default]
    Unknown,
}

/// One search hit of a gained region against the reference genome:
/// blastn for the nucleotide tiers, tblastx for the translated tier.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct GainedBlastHit {
    /// Reference sequence the hit is on.
    pub ref_seqid: String,
    /// 1-based inclusive reference interval of the hit, ordered
    /// low..high regardless of the hit's strand.
    pub ref_start: u64,
    pub ref_end: u64,
    /// Percent identity over the aligned length: nucleotide identity
    /// for blastn, amino-acid identity for tblastx.
    pub identity: f64,
    /// Alignment length in bases (nucleotide, even for tblastx,
    /// which reports the aligned span in nucleotide coordinates).
    pub length: u64,
    /// 1-based inclusive interval on the region itself (relative to the
    /// sequence that was searched, not the whole query contig).
    pub qry_start: u64,
    pub qry_end: u64,
    pub evalue: f64,
    pub bitscore: f64,
    /// The reference genes the hit's interval overlaps, formatted for
    /// display ("LM4B_RS13070 (rlmN)"). A hit the user must judge means
    /// little until it says what it hits.
    pub genes: Vec<String>,
}

/// The best reference-protein match of one predicted gene inside a
/// gained region: the closest thing to a name an unannotated query
/// genome can be given locally.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct OrfMatch {
    pub locus_tag: String,
    /// RefSeq protein accession when the reference annotation has one,
    /// for the NCBI protein link.
    pub protein_id: String,
    /// The gene symbol or product, whichever the annotation carries.
    pub label: String,
    /// Amino-acid identity over the aligned part, percent.
    pub identity: f64,
    /// Share of the predicted gene covered by the alignment, percent.
    pub coverage: f64,
    pub evalue: f64,
    /// Raw bitscore, for ranking matches of equal E-value.
    #[serde(default)]
    pub bitscore: f64,
}

/// One predicted gene of a gained region, identified or not, with the
/// nucleotide sequence it was called from so the client can link out
/// for the matches no local database can name.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct IdentifiedOrf {
    /// 1-based inclusive, in query contig coordinates, matching the
    /// ORFs of the gained row (positional, same order).
    pub start: u64,
    pub end: u64,
    pub strand: i8,
    /// Best reference-protein match; `None` means no similar protein
    /// in the reference, which is expected for the true gains.
    #[serde(rename = "match")]
    pub best: Option<OrfMatch>,
    /// The predicted gene's nucleotide sequence, plus strand.
    pub seq: String,
}

/// The gene-level answer to "what is gained here": each predicted gene
/// of one region with its best match among the reference's own
/// proteins, searched as translated DNA (blastx).
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct GainedIdentify {
    pub qry_seqid: String,
    /// 1-based inclusive coordinates on the query contig.
    pub start: u64,
    pub end: u64,
    /// Positionally parallel to the gained row's `orfs`.
    pub orfs: Vec<IdentifiedOrf>,
    /// The region's own sequence, for linking out to NCBI BLAST.
    pub region_seq: String,
}

/// The reference back-check of one gained region: the region's sequence
/// searched against the reference genome, which is the evidence "no
/// alignment to the reference" cannot supply on its own.
///
/// Three searches are reported. The nucleotide search comes in two tiers
/// split at E = 1e-5, so borderline matches can be seen instead of
/// silently dropped. The translated search (tblastx) only runs when the
/// strong nucleotide tier is empty, because its whole purpose is the
/// second opinion on "found nothing but I don't believe it" - a
/// divergent coding homolog that 11-mer nucleotide seeding misses. The
/// longest exact match is tool-free and cutoff-free: the longest stretch
/// of bases the region and the reference share anywhere, which bounds
/// every question about short exact remnants at once.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct GainedVerify {
    /// The region that was searched, echoed so a cached answer can be
    /// matched to its question.
    pub qry_seqid: String,
    /// 1-based inclusive coordinates on the query contig.
    pub start: u64,
    pub end: u64,
    /// Nucleotide hits with E <= 1e-5, ordered by bitscore, best first.
    /// Empty means the strong search found nothing similar.
    pub hits: Vec<GainedBlastHit>,
    /// Nucleotide hits with 1e-5 < E <= 10: below the detection
    /// threshold the verdicts are built on, kept because a scientist
    /// asked to believe an absence deserves to see what was almost
    /// there. Usually noise.
    pub weak_hits: Vec<GainedBlastHit>,
    /// Translated (tblastx) hits, E <= 1e-5. `None` means the search
    /// did not run - the nucleotide search already answered, or the
    /// region exceeds the length cap - which is not the same as having
    /// looked and found none.
    pub tx_hits: Option<Vec<GainedBlastHit>>,
    /// Plain-language reason the translated search did not run.
    pub tx_note: Option<String>,
    /// Longest run of bases the region shares, exactly, with any
    /// position of the reference genome. For unrelated DNA of these
    /// sizes chance alone gives ~log4(region * reference) bases, so a
    /// low number here is expected, not informative; a high one is.
    pub longest_exact_bp: u64,
    /// Where that longest exact match sits on the reference; empty
    /// seqid when the region shares nothing at all with the reference
    /// (not even one base).
    pub longest_exact_seqid: String,
    /// 1-based inclusive interval of the example match. Several may
    /// exist; one is reported.
    pub longest_exact_start: u64,
    pub longest_exact_end: u64,
    /// The same example match's interval on the region, so the match
    /// can be located without a search.
    pub longest_exact_qry_start: u64,
    pub longest_exact_qry_end: u64,
}

/// One row of the strict panel recheck (per query).
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct PanelRow {
    pub gene_id: String,
    pub qlen: u64,
    pub cov_pct: f64,
    pub identity: f64,
    pub best_evalue: String,
    pub call: Call,
    /// The best hit's place in the query, "contig:start-end(+|-)"; just
    /// the contig for runs made before positions were kept, empty when
    /// the gene was not found.
    #[serde(default)]
    pub qry_locus: String,
    /// The panel record behind the best hit: the gene's own id, or one of
    /// its other variants ("cadA__v2"); empty when nothing hit.
    #[serde(default)]
    pub variant: String,
    /// How many sequences (variants) the panel holds for this gene.
    #[serde(default)]
    pub n_variants: usize,
    /// Where the matching variant's sequence came from ("AMRFinderPlus
    /// (cadA_Lm)", "NCBI Nucleotide (L28104.1:158-2293)"), filled when
    /// the table is read; empty when nothing hit.
    #[serde(default)]
    pub variant_source: String,
    /// The protein-level search, run for a gene not found in full at DNA
    /// level: a different variant of the gene shows up here when its DNA
    /// has drifted too far for the nucleotide search.
    #[serde(default)]
    pub protein: Option<ProteinHit>,
    /// Not found in full at DNA level, yet the protein search found a
    /// close relative over most of its length: possibly a variant of the
    /// gene the panel does not hold.
    #[serde(default)]
    pub variant_warning: bool,
}

/// A panel gene's best protein-level (tblastn) match in one query.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
pub struct ProteinHit {
    /// The panel record whose translation matched.
    pub variant: String,
    /// Amino-acid identity over the aligned part, percent.
    pub identity: f64,
    /// Share of the panel protein covered, percent.
    pub coverage: f64,
    /// "contig:start-end(+|-)" in the query.
    pub locus: String,
    /// What the match is when it is a known gene rather than an unknown
    /// variant: "emrC, another panel gene", or a gene in DNA the strain
    /// shares with the reference. Empty when unexplained.
    #[serde(default)]
    pub explained_by: String,
}

/// Protein identity and coverage at which a gene missing at DNA level is
/// flagged as possibly present in another variant. Family relatives (SMR
/// pumps such as emrC and sugE) sit near 40 %, other variants of one gene
/// (the cadA of Tn5422 against cadA_Lm) near 70 %.
pub const VARIANT_MIN_PROTEIN_IDENTITY: f64 = 60.0;
/// A partial DNA match shorter than this share of the gene, whose protein
/// match is only a distant relative, is a conserved stretch of another
/// gene of the family (the ATPase core cadA shares with every Listeria
/// metal pump), and the gene is called absent.
pub const SHORT_MATCH_MAX_COVERAGE: f64 = 30.0;
pub const VARIANT_MIN_PROTEIN_COVERAGE: f64 = 80.0;

impl ProteinHit {
    pub fn suggests_variant(&self) -> bool {
        self.identity >= VARIANT_MIN_PROTEIN_IDENTITY
            && self.coverage >= VARIANT_MIN_PROTEIN_COVERAGE
    }
}

/// Presence/absence across all queries of a run.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct MatrixRow {
    pub locus_tag: String,
    pub symbol: String,
    pub biotype: String,
    pub seqid: String,
    pub start: u64,
    pub end: u64,
    /// One entry per query, same order as the run's query list.
    pub calls: Vec<Call>,
    pub cov_pcts: Vec<f64>,
}

/// One panel gene across all queries of a run.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct PanelMatrixRow {
    pub gene_id: String,
    pub qlen: u64,
    /// One entry per query, same order as the run's query list.
    pub calls: Vec<Call>,
    pub cov_pcts: Vec<f64>,
    pub identities: Vec<f64>,
    pub loci: Vec<String>,
    /// Per query: possibly present as a variant the panel lacks.
    #[serde(default)]
    pub variant_warnings: Vec<bool>,
    /// Per query, for a gene not present in full: the protein identity of
    /// the closest related gene, percent; None when there is none.
    #[serde(default)]
    pub related_identities: Vec<Option<f64>>,
}

/// One query contig and how much of it aligns to the reference.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct ContigStat {
    pub seqid: String,
    pub length: u64,
    pub aligned_bp: u64,
}

/// A predicted gene near a panel hit, taken from the gained region ORFs.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct ContextGene {
    /// Query contig coordinates, 1-based inclusive.
    pub start: u64,
    pub end: u64,
    pub strand: i8,
    /// Reference-given name, else the NCBI-given one, else empty (novel).
    /// The predicted gene that is the panel hit carries the panel gene's
    /// name instead.
    pub label: String,
    /// Where the gene comes from: "annotation" (a reference gene placed
    /// through the alignment), the panel hit itself ("panel"), or a gene
    /// predicted in a stretch the reference lacks, named from the
    /// reference proteome ("reference"), by NCBI ("ncbi"), or unnamed ("").
    pub source: String,
    /// The protein a predicted gene was named after. Set on the panel hit,
    /// whose own label is the panel gene: there it is only the closest
    /// relative found, worth showing with its identity.
    #[serde(default)]
    pub match_label: String,
    /// Amino-acid identity of the predicted gene to the protein it was
    /// named after, percent.
    #[serde(default)]
    pub match_identity: Option<f64>,
    /// Share of the predicted gene covered by that match, percent, when
    /// the search reported it.
    #[serde(default)]
    pub match_coverage: Option<f64>,
    /// The reference locus tag, for "annotation" genes.
    #[serde(default)]
    pub locus_tag: String,
    /// Protein accession of an "annotation" gene, to look it up in the
    /// reference library.
    #[serde(skip)]
    pub protein_id: String,
    /// Gene name of the library variant group the gene belongs to, when
    /// the annotation gives it none of its own ("gad").
    #[serde(default)]
    pub group_name: String,
    /// The same group's gene(s) in the species' RefSeq reference genome
    /// ("lmo2363"), telling apart genes of one function.
    #[serde(default)]
    pub reference_loci: Vec<String>,
    /// That genome's strain ("EGD-e").
    #[serde(default)]
    pub reference_strain: String,
    /// Bases between this gene and the panel hit; 0 when they overlap.
    pub distance: u64,
    /// The gene overlaps the panel hit itself.
    pub is_hit: bool,
    /// The name points at a mobile element (transposase, resolvase,
    /// integrase, plasmid replication...).
    pub mobile: bool,
    pub partial: bool,
}

/// Where a panel gene sits in one query genome.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct PanelContext {
    pub query_id: i64,
    pub query_name: String,
    pub gene_id: String,
    /// Where the panel's sequence for the gene came from: "AMRFinderPlus
    /// (emrC_Lis)", "VFDB (...)", "NCBI Nucleotide (...)", "reference
    /// genome (lmo0200)" or "your panel FASTA"; empty when unknown.
    #[serde(default)]
    pub panel_source: String,
    /// What the hit is when it is not a full match: a partial DNA match,
    /// or only a protein-level relative (a possible other variant).
    /// Empty for a gene found in full.
    #[serde(default)]
    pub match_note: String,
    /// For a gene not present in full: the name of the closest related
    /// gene where the strain's annotation or prediction gives one.
    #[serde(default)]
    pub related_gene: String,
    /// Its protein identity to the gene of interest, percent.
    #[serde(default)]
    pub related_identity: Option<f64>,
    pub call: Call,
    pub cov_pct: f64,
    pub identity: f64,
    /// The hit, empty contig when the gene was not found.
    pub contig: Option<ContigStat>,
    pub hit_start: u64,
    pub hit_end: u64,
    pub hit_strand: i8,
    pub n_contigs: usize,
    pub genome_bp: u64,
    /// The gained region (query stretch with no reference alignment)
    /// containing the hit, if any.
    pub region: Option<GainedRow>,
    /// Predicted genes within the context window, in contig order.
    pub genes: Vec<ContextGene>,
    /// Half-width of the window `genes` covers, bp.
    pub window: u64,
    /// "plasmid", "chromosome_insertion", "chromosome_shared",
    /// "unplaced" or "not_found".
    pub verdict: String,
    /// The plain language verdict with the numbers behind it.
    pub verdict_text: String,
    /// Names of the mobile-element genes near the hit.
    pub mobile_markers: Vec<String>,
    /// Why `genes` is empty although the hit is in a gained region
    /// (the gene finder did not run), empty otherwise.
    #[serde(default)]
    pub genes_note: String,
}

/// A panel gene aligned base by base to where it was found in one strain,
/// for comparing a Present strain with a Partial one side by side.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct PanelAlignmentView {
    pub query_id: i64,
    pub query_name: String,
    pub gene_id: String,
    pub call: Call,
    /// The panel record aligned, and where its sequence came from.
    pub variant: String,
    pub variant_source: String,
    pub panel_len: u64,
    pub panel_seq: String,
    /// Where in the strain: the DNA hit, or the protein-level relative.
    pub contig: String,
    pub strand: i8,
    /// 1-based inclusive range of the panel gene that aligned.
    pub panel_start: u64,
    pub panel_end: u64,
    /// The aligned stretch on the contig, 1-based inclusive.
    pub contig_start: u64,
    pub contig_end: u64,
    pub identity: f64,
    /// Share of the panel gene inside the aligned stretch, percent.
    pub panel_coverage: f64,
    pub mismatches: u64,
    pub gap_bases: u64,
    /// Aligned rows, '-' for gaps; panel gene on top.
    pub panel_row: String,
    pub strain_row: String,
    /// The dialog's note on a match that is not a full one.
    pub match_note: String,
}

/// How much of an element one query genome holds.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct ElementHit {
    pub query_id: i64,
    pub query_name: String,
    /// The query's panel call for the gene, when the run has a panel.
    pub call: Option<Call>,
    pub covered_bp: u64,
    pub covered_pct: f64,
    /// Separate stretches of the element found (after merging).
    pub pieces: usize,
    pub largest_piece: u64,
    /// Identity over the matched bases, length-weighted, percent.
    pub identity: f64,
    /// Query contigs the element matched, most matched bases first.
    pub contigs: Vec<String>,
    pub genome_bp: u64,
    /// For a gene not present in full in this query: the protein identity
    /// of the closest related gene, percent.
    #[serde(default)]
    pub related_identity: Option<f64>,
    /// The stretches of the element found, in element coordinates, in
    /// element order.
    #[serde(default)]
    pub covered: Vec<ElementPiece>,
}

/// One stretch of an element found in a genome.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, Default)]
pub struct ElementPiece {
    /// 1-based inclusive, on the element.
    pub start: u64,
    pub end: u64,
    /// Identity over the stretch, length-weighted, percent.
    pub identity: f64,
}

/// An annotated gene of an element record (a library or NCBI plasmid).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, Default)]
pub struct ElementGene {
    /// 1-based inclusive, on the element.
    pub start: u64,
    pub end: u64,
    /// 1 or -1.
    pub strand: i8,
    pub name: String,
    pub locus_tag: String,
    pub product: String,
}

/// A whole-element comparison across the queries of a run.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct ElementReport {
    pub gene_id: String,
    pub source_query_id: i64,
    /// "contig" (the contig carrying the gene), "region" (the inserted
    /// stretch carrying it, "contig:start-end", when the contig is
    /// otherwise shared with the reference) or "accession".
    pub element_kind: String,
    /// The contig name or the accession.
    pub element_name: String,
    /// NCBI title for an accession, empty for a contig.
    pub element_title: String,
    pub element_len: u64,
    pub hits: Vec<ElementHit>,
    /// The record's annotated genes; empty for a strain's own contig,
    /// which carries no annotation.
    #[serde(default)]
    pub genes: Vec<ElementGene>,
}

/// A reference-library record carrying a panel gene: a whole element to
/// compare the strains against.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, Default)]
pub struct ElementRecord {
    pub accession: String,
    /// "Listeria monocytogenes L1551 plasmid pLM33"
    pub title: String,
    /// "plasmid" or "chromosome".
    pub kind: String,
    pub length: u64,
    /// The library genomes it stands for (itself and its near-duplicates).
    pub genomes: u64,
}

/// The library records carrying a panel gene.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, Default)]
pub struct ElementRecords {
    /// "Listeria library 2026-10-02.2"; empty without a library.
    pub library: String,
    pub records: Vec<ElementRecord>,
    /// Largest record the comparison takes; a chromosome is bigger.
    pub max_bp: u64,
    /// Why the list could not be made, when it could not.
    pub note: String,
}

/// One NCBI record carrying a panel gene, from the gene origin lookup.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct OriginRecord {
    pub accession: String,
    pub title: String,
    pub length: u64,
    /// "plasmid", "chromosome" or "contig" (a draft-assembly piece that
    /// NCBI does not place on either).
    pub kind: String,
    pub identity: f64,
    /// Share of the panel gene the best match covers, percent.
    pub coverage: f64,
}

/// Where a panel gene usually occurs, from an NCBI BLAST search of its
/// sequence. The answer is evidence about other genomes, never about the
/// strains of the run.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct GeneOrigin {
    pub gene_id: String,
    /// The panel sequence searched, by where it came from ("NCBI
    /// Nucleotide (L28104.1:158-2293)"): the gene's variant that matched
    /// most often in the run, so the answer is about the strains' copy.
    #[serde(default)]
    pub searched_with: String,
    /// "running", "done" or "error".
    pub state: String,
    /// What was searched: a genus ("Listeria") or "all bacteria".
    pub scope: String,
    /// Why the scope is what it is (e.g. widened for too few matches).
    pub scope_note: String,
    /// Records matching the gene over (nearly) its whole length.
    pub n_matches: usize,
    pub n_plasmid: usize,
    pub n_chromosome: usize,
    pub n_contig: usize,
    /// Complete plasmids carrying the gene, smallest first: ready-made
    /// elements for the whole-element comparison.
    pub plasmids: Vec<OriginRecord>,
    /// A few chromosome records carrying it, for contrast.
    pub chromosomes: Vec<OriginRecord>,
    pub message: String,
    /// The genus held too few matches; searching all bacteria (slow) is
    /// offered.
    #[serde(default)]
    pub can_widen: bool,
}

/// One gene found by the resistance/virulence screen.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct ScreenHit {
    /// "AMRFinderPlus" or "VFDB".
    pub source: String,
    pub gene: String,
    pub product: String,
    /// AMRFinderPlus: "AMR", "STRESS" or "VIRULENCE"; VFDB: "VIRULENCE".
    pub kind: String,
    /// AMRFinderPlus subtype (ACID, BIOCIDE, METAL, POINT ...) or the
    /// VFDB category ("Exotoxin", "Invasion" ...).
    pub category: String,
    /// AMRFinderPlus class ("QUATERNARY AMMONIUM", "CADMIUM", ...) or
    /// the VFDB virulence factor ("Listeriolysin O (LLO)").
    pub class: String,
    pub contig: String,
    pub start: u64,
    pub end: u64,
    pub strand: i8,
    pub identity: f64,
    /// Share of the reference gene covered, percent.
    pub coverage: f64,
    /// The reference the hit was matched to (accession + name).
    pub reference: String,
    /// AMRFinderPlus method (EXACTX, BLASTX, PARTIALX, HMM ...) or "blastn".
    pub method: String,
}

/// Whether the screen ran for a query, and why not.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(tag = "state", content = "reason", rename_all = "snake_case")]
pub enum ScreenStatus {
    Done,
    /// One of the two resources is missing; carries which.
    Partial(String),
    /// Switched off, or the tools/database are missing; carries why.
    Unavailable(String),
    /// The run predates the screen.
    #[default]
    Unknown,
}

/// One screened gene across all queries of a run.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct ScreenMatrixRow {
    pub source: String,
    pub gene: String,
    pub product: String,
    pub kind: String,
    pub category: String,
    pub class: String,
    /// Identity of the best hit per query, same order as the run's query
    /// list; None where the gene was not found.
    pub identities: Vec<Option<f64>>,
}

/// Alignment block on the reference, for the genome viewer.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WgaBlock {
    pub ref_seqid: String,
    pub ref_start: u64,
    pub ref_end: u64,
    pub qry_seqid: String,
    pub qry_start: u64,
    pub qry_end: u64,
    pub qry_rev: bool,
    /// 0-100
    pub identity: f64,
}

/// Reference gene for the genome viewer track.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct WgaGene {
    pub locus_tag: String,
    pub symbol: String,
    pub biotype: String,
    pub seqid: String,
    pub start: u64,
    pub end: u64,
    pub strand: i8,
    /// Function annotation from the GFF `product` attribute.
    #[serde(default)]
    pub product: String,
    /// RefSeq/GenBank protein accession of the coding gene ("WP_...").
    #[serde(default)]
    pub protein_id: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WgaQuery {
    pub query_id: i64,
    pub query_name: String,
    pub blocks: Vec<WgaBlock>,
    /// Presence call per gene, aligned with `WgaData::genes`.
    #[serde(default)]
    pub calls: Vec<Call>,
    /// Coverage percent per gene, aligned with `WgaData::genes`.
    #[serde(default)]
    pub cov_pcts: Vec<f64>,
    /// Best block identity per gene (0-100), aligned with `WgaData::genes`.
    #[serde(default)]
    pub identities: Vec<f64>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WgaData {
    /// Reference sequence lengths.
    pub reference: Vec<(String, u64)>,
    pub genes: Vec<WgaGene>,
    pub queries: Vec<WgaQuery>,
}

/// A single nucleotide polymorphism at a reference position: the query
/// carries a different base. Bases are ASCII bytes, and the query base
/// is already in reference orientation (reverse complemented when the
/// aligning block is on the reverse strand).
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct SnpEvent {
    /// 1-based reference position.
    pub pos: u64,
    pub r: u8,
    pub q: u8,
}

/// A stretch of reference bases missing from the query (a gap in the
/// query row inside an aligned block).
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct DelEvent {
    /// 1-based reference position of the first deleted base.
    pub pos: u64,
    pub len: u64,
}

/// Query bases inserted between two adjacent aligned reference
/// positions. `pos` is the reference position after which the bases
/// sit (0 = before the first reference base of the block). The
/// sequence is in reference orientation.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct InsEvent {
    pub pos: u64,
    pub seq: String,
}

/// Variant events of one query against one reference sequence.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct AlignmentEvents {
    pub snps: Vec<SnpEvent>,
    pub dels: Vec<DelEvent>,
    pub ins: Vec<InsEvent>,
}

/// Wire form of [AlignmentEvents]: parallel arrays instead of one
/// object per event. A divergent query carries >100k SNPs, and the
/// object form made the alignment endpoint answer ~4 MB per query
/// (~60 MB per run) that the browser's JSON parser ground through on
/// the main thread for seconds. Arrays of plain numbers serialize
/// smaller and parse several times faster. Only the HTTP response
/// uses this; variants.json on disk keeps the object form.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
pub struct AlignmentEventsColumnar {
    pub snp_pos: Vec<u64>,
    pub snp_ref: Vec<u8>,
    pub snp_qry: Vec<u8>,
    pub del_pos: Vec<u64>,
    pub del_len: Vec<u64>,
    pub ins_pos: Vec<u64>,
    pub ins_seq: Vec<String>,
}

impl From<AlignmentEvents> for AlignmentEventsColumnar {
    fn from(ev: AlignmentEvents) -> Self {
        let mut out = AlignmentEventsColumnar {
            snp_pos: Vec::with_capacity(ev.snps.len()),
            snp_ref: Vec::with_capacity(ev.snps.len()),
            snp_qry: Vec::with_capacity(ev.snps.len()),
            del_pos: Vec::with_capacity(ev.dels.len()),
            del_len: Vec::with_capacity(ev.dels.len()),
            ins_pos: Vec::with_capacity(ev.ins.len()),
            ins_seq: Vec::with_capacity(ev.ins.len()),
        };
        for s in ev.snps {
            out.snp_pos.push(s.pos);
            out.snp_ref.push(s.r);
            out.snp_qry.push(s.q);
        }
        for d in ev.dels {
            out.del_pos.push(d.pos);
            out.del_len.push(d.len);
        }
        for i in ev.ins {
            out.ins_pos.push(i.pos);
            out.ins_seq.push(i.seq);
        }
        out
    }
}

/// One query with its alignment blocks and per-base variant events
/// for the whole-genome alignment viewer.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct AlignmentQuery {
    pub query_id: i64,
    pub query_name: String,
    pub blocks: Vec<WgaBlock>,
    /// Events keyed by reference seqid.
    #[serde(default)]
    pub events: std::collections::BTreeMap<String, AlignmentEvents>,
}

/// Wire form of [AlignmentQuery] (see [AlignmentEventsColumnar]):
/// same fields, columnar events.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct AlignmentQueryColumnar {
    pub query_id: i64,
    pub query_name: String,
    pub blocks: Vec<WgaBlock>,
    #[serde(default)]
    pub events: std::collections::BTreeMap<String, AlignmentEventsColumnar>,
}

/// Wire form of [AlignmentData] (see [AlignmentEventsColumnar]).
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct AlignmentDataColumnar {
    /// Reference sequence lengths (seqid, length), sorted.
    pub reference: Vec<(String, u64)>,
    pub queries: Vec<AlignmentQueryColumnar>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct AlignmentData {
    /// Reference sequence lengths (seqid, length), sorted.
    pub reference: Vec<(String, u64)>,
    pub queries: Vec<AlignmentQuery>,
}

/// Hover preview + alignment rows for one gene (MSA viewer).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GeneBlock {
    /// Reference coordinates of the aligned slice (1-based inclusive).
    pub ref_start: u64,
    pub ref_end: u64,
    /// Query coordinates of the aligned slice (on the forward query
    /// sequence; original orientation preserved through qry_rev).
    pub qry_start: u64,
    pub qry_end: u64,
    pub qry_rev: bool,
    pub identity: f64,
    pub ref_seq: String,
    pub qry_seq: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GeneQueryAlignment {
    pub query_id: i64,
    pub query_name: String,
    pub call: Call,
    pub cov_pct: f64,
    pub best_identity: f64,
    pub mismatches: u64,
    pub indels: u64,
    /// Alignment blocks overlapping the gene, ordered by reference position.
    pub blocks: Vec<GeneBlock>,
    /// Unaligned reference stretches between blocks, (start, end) inclusive.
    pub unaligned: Vec<(u64, u64)>,
    /// Premature stop codons found in the query sequence of this gene.
    pub premature_stops: Vec<PrematureStop>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PrematureStop {
    /// Codon index (0-based from the gene start).
    pub codon_index: u64,
    /// Amino acid position (1-based) where the stop appears.
    pub aa_position: u64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GeneDetail {
    pub locus_tag: String,
    pub symbol: String,
    #[serde(default)]
    pub protein_id: String,
    pub biotype: String,
    pub seqid: String,
    pub start: u64,
    pub end: u64,
    pub strand: i8,
    pub length: u64,
    pub reference_seq: String,
    pub queries: Vec<GeneQueryAlignment>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum RunStatus {
    Queued,
    Running,
    Succeeded,
    Failed,
}

impl RunStatus {
    pub fn as_str(&self) -> &'static str {
        match self {
            RunStatus::Queued => "queued",
            RunStatus::Running => "running",
            RunStatus::Succeeded => "succeeded",
            RunStatus::Failed => "failed",
        }
    }
}

/// Table sorting / filtering / pagination query, shared by the table
/// endpoints and the export endpoints so exports match the visible view.
#[derive(Debug, Clone, Deserialize)]
pub struct TableQuery {
    #[serde(default)]
    pub query_id: Option<i64>,
    #[serde(default = "default_page")]
    pub page: u64,
    #[serde(default = "default_page_size")]
    pub page_size: u64,
    #[serde(default)]
    pub sort_by: Option<String>,
    #[serde(default)]
    pub sort_dir: Option<String>,
    #[serde(default)]
    pub search: Option<String>,
    /// Filter by call: present | partial | absent.
    ///
    /// Tables without a call column reuse this field for their own one
    /// filter rather than growing the struct the export endpoints share:
    /// the matrix takes `not_present`, and gained regions take
    /// `anchored` | `unanchored`.
    #[serde(default)]
    pub call: Option<String>,
    /// For exports: comma separated column keys; empty = all.
    #[serde(default)]
    pub cols: Option<String>,
}

fn default_page() -> u64 {
    0
}
fn default_page_size() -> u64 {
    200
}

/// A page of table rows plus total count, so the UI can paginate.
#[derive(Debug, Clone, Serialize)]
pub struct Page<T> {
    pub rows: Vec<T>,
    pub total: u64,
}
