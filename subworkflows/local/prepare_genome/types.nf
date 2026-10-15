// Documentation only. Casting record(...) to these types corrupts remote Path fields.
record PreparedGenome {
    alignment_index: Path?
    bbsplit_index: Path?
    dict: Path?
    fai: Path?
    bcftools_annotations_tbi: Path?
    dbsnp_tbi: Path?
    germline_resource_tbi: Path?
    known_indels_tbi: Path?
    known_snps_tbi: Path?
    pon_tbi: Path?
    msisensor2_models: Path?
    msisensorpro_scan: Path?
    chr_dir: Path?
}
