//! GFF parsing tests (pure, no external tools).

use straincompass_engine::gff::parse_gff_str;

#[test]
fn gene_features_get_product_from_matching_cds() {
    let gff = "##gff-version 3
chr1\tsrc\tgene\t101\t500\t.\t+\t.\tlocus_tag=TAG1;gene=symA;gene_biotype=protein_coding
chr1\tsrc\tCDS\t101\t500\t.\t+\t0\tlocus_tag=TAG1;product=beta-lactamase family hydrolase;protein_id=WP_000001.1
chr1\tsrc\tgene\t601\t900\t.\t-\t.\tlocus_tag=TAG2;gene_biotype=protein_coding;product=toxin ABC transporter
";
    let genes = parse_gff_str(gff).unwrap();
    assert_eq!(genes.len(), 2);
    assert_eq!(genes[0].locus_tag, "TAG1");
    assert_eq!(genes[0].product, "beta-lactamase family hydrolase");
    assert_eq!(genes[0].protein_id, "WP_000001.1");
    // product directly on the gene feature is kept
    assert_eq!(genes[1].product, "toxin ABC transporter");
    // no CDS child -> no protein accession
    assert_eq!(genes[1].protein_id, "");
    assert_eq!(genes[1].strand, -1);
}

#[test]
fn cds_fallback_merges_by_locus_tag_and_keeps_product() {
    let gff = "##gff-version 3
chr1\tsrc\tCDS\t101\t300\t.\t+\t0\tlocus_tag=TAG9;product=hypothetical protein;protein_id=WP_000009.1
chr1\tsrc\tCDS\t401\t600\t.\t+\t0\tlocus_tag=TAG9;product=hypothetical protein
chr1\tsrc\ttRNA\t1001\t1080\t.\t+\t.\tlocus_tag=TRNA1;product=tRNA-Ala
";
    let genes = parse_gff_str(gff).unwrap();
    assert_eq!(genes.len(), 2);
    let tag9 = genes.iter().find(|g| g.locus_tag == "TAG9").unwrap();
    assert_eq!(tag9.start, 101);
    assert_eq!(tag9.end, 600);
    assert_eq!(tag9.product, "hypothetical protein");
    assert_eq!(tag9.protein_id, "WP_000009.1");
    let trna = genes.iter().find(|g| g.locus_tag == "TRNA1").unwrap();
    assert_eq!(trna.biotype, "tRNA");
    assert_eq!(trna.product, "tRNA-Ala");
    assert_eq!(trna.protein_id, "");
}

#[test]
fn missing_product_is_empty() {
    let gff = "##gff-version 3
chr1\tsrc\tgene\t101\t500\t.\t+\t.\tlocus_tag=TAG1;gene_biotype=protein_coding
";
    let genes = parse_gff_str(gff).unwrap();
    assert_eq!(genes[0].product, "");
}
