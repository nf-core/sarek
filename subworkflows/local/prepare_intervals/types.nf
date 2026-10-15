// Documentation only. Casting record(...) to these types corrupts remote Path fields.
record PreparedIntervals {
    split_bed: List<Path>?
    split_bed_gz: List<Path>?
    combined_bed: Path?
    combined_bed_gz: Path?
}
