// Documentation only. Casting record(...) to these types corrupts remote Path fields.
record SarekMultiqc {
    id: String
    meta: Map
    report: Path
    data: Path
    plots: Path?
}
