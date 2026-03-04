#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
BIN="${PROJECT_DIR}/target/release/ironqc"
INPUT_DIR="${SCRIPT_DIR}/input/raredisease"
OUTPUT_DIR="${SCRIPT_DIR}/output/raredisease"
REF_DIR="${SCRIPT_DIR}/reference/raredisease"

REF="${INPUT_DIR}/reference.fasta"
FAI="${INPUT_DIR}/reference.fasta.fai"

if [ ! -f "$BIN" ]; then
    echo "Release binary not found. Building..."
    (cd "$PROJECT_DIR" && cargo build --release)
fi

if [ ! -d "$INPUT_DIR" ]; then
    echo "Test data not found at $INPUT_DIR"
    echo "Download raredisease test data first:"
    echo "  mkdir -p $INPUT_DIR"
    echo "  curl -sL -o $INPUT_DIR/hugelymodelbat_sorted_md.bam https://raw.githubusercontent.com/nf-core/test-datasets/raredisease/testdata/hugelymodelbat_sorted_md.bam"
    echo "  curl -sL -o $INPUT_DIR/hugelymodelbat_sorted_md.bam.bai https://raw.githubusercontent.com/nf-core/test-datasets/raredisease/testdata/hugelymodelbat_sorted_md.bam.bai"
    echo "  curl -sL -o $INPUT_DIR/reference.fasta https://raw.githubusercontent.com/nf-core/test-datasets/raredisease/reference/reference.fasta"
    echo "  curl -sL -o $INPUT_DIR/reference.fasta.fai https://raw.githubusercontent.com/nf-core/test-datasets/raredisease/reference/reference.fasta.fai"
    exit 1
fi

mkdir -p "$OUTPUT_DIR/stats" "$OUTPUT_DIR/mosdepth"

echo "=== ironqc benchmark ==="
echo "Binary: $BIN"
echo "Input:  $INPUT_DIR"
echo ""

for bam in "$INPUT_DIR"/*_sorted_md.bam; do
    [ -f "$bam" ] || continue
    name=$(basename "$bam" .bam)
    echo "--- $name ---"

    echo -n "  stats:    "
    /usr/bin/time -p "$BIN" stats \
        --reference "$REF" \
        --prefix "$OUTPUT_DIR/stats/$name" \
        "$bam" 2>&1 | grep real | awk '{print $2 "s"}'

    echo -n "  mosdepth: "
    /usr/bin/time -p "$BIN" mosdepth \
        --fasta "$REF" \
        "$bam" \
        "$OUTPUT_DIR/mosdepth/$name" 2>&1 | grep real | awk '{print $2 "s"}'

    mkdir -p "$OUTPUT_DIR/bundle/${name}_indexcov"
    echo -n "  bundle:   "
    /usr/bin/time -p "$BIN" bundle \
        --reference "$REF" \
        --fai "$FAI" \
        --prefix "$OUTPUT_DIR/bundle/$name" \
        --indexcov-dir "$OUTPUT_DIR/bundle/${name}_indexcov" \
        "$bam" 2>&1 | grep real | awk '{print $2 "s"}'

    echo ""
done

if [ -d "$REF_DIR" ]; then
    echo "=== Validation against upstream ==="
    for stats_file in "$OUTPUT_DIR"/stats/*.stats; do
        [ -f "$stats_file" ] || continue
        name=$(basename "$stats_file" .stats)
        ref_file="$REF_DIR/${name%.sorted_md}.samtools.stats"
        if [ -f "$ref_file" ]; then
            echo "Comparing stats: $name"
            diff <(grep '^SN' "$stats_file" | cut -f2,3) \
                 <(grep '^SN' "$ref_file" | cut -f2,3) || echo "  DIFFERENCES FOUND"
        fi
    done
fi

echo "=== Done ==="
