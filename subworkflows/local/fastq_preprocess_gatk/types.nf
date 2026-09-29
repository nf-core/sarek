// Documentation only. Casting record(...) to these types corrupts remote Path fields.
record FastpResult {
    id: String
    meta: Map
    reads: List<Path>
    html: Path
    json: Path
    log: Path
}

record BbsplitResult {
    id: String
    meta: Map
    reads: List<Path>
    stats: Path
}

record PreprocessingAlignment {
    id: String
    meta: Map
    alignment: Path
    index: Path
}

record MarkduplicatesAlignment {
    id: String
    meta: Map
    directory: String
    alignment: Path
    index: Path
}

record RecalibrationTable {
    id: String
    meta: Map
    table: Path
}
