#!/usr/bin/env bash
#
# Build a small DeepVariant joint-genotyping review dataset from 1000 Genomes Phase 3.
#
# Downloads 3 unrelated samples, subsets each alignment to a single chromosome, indexes
# them, and writes a sarek samplesheet. No data is committed to the repo - run this to
# (re)create the review dataset locally.
#
# 1000 Genomes Phase 3 alignments are GRCh37 / human_g1k_v37 (numeric chromosome names,
# e.g. "20"), so run sarek with --genome GATK.GRCh37.
#
# Requirements: samtools (with remote/htslib support), curl.
#
# Usage:
#   tests/scripts/download_joint_test_data.sh [OUTDIR] [CHROMOSOME]
#   tests/scripts/download_joint_test_data.sh test_data_1000g 20
#
set -euo pipefail

OUTDIR="${1:-test_data_1000g}"
REGION="${2:-20}"                          # single chromosome (GRCh37 naming, no "chr")
BASE="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data"

# 3 unrelated GBR samples, all low-coverage Phase 3 alignments
SAMPLES=(HG00096 HG00097 HG00099)
SEXES=(XY XX XX)

command -v samtools >/dev/null 2>&1 || { echo "ERROR: samtools is required" >&2; exit 1; }

mkdir -p "$OUTDIR"
ABS_OUTDIR="$(cd "$OUTDIR" && pwd)"
CSV="$OUTDIR/joint_germline_deepvariant_1000g.csv"
echo "patient,sex,status,sample,bam,bai" > "$CSV"

for i in "${!SAMPLES[@]}"; do
    S="${SAMPLES[$i]}"
    SEX="${SEXES[$i]}"
    BAM_URL="$BASE/$S/alignment/$S.mapped.ILLUMINA.bwa.GBR.low_coverage.20120522.bam"
    OUT="$ABS_OUTDIR/${S}.chr${REGION}.bam"

    echo ">> ${S}: extracting chromosome ${REGION} from ${BAM_URL}"
    # htslib streams the remote BAM + its .bai and pulls only the requested region,
    # so we never download the whole genome-wide alignment.
    samtools view -b "$BAM_URL" "$REGION" > "$OUT"
    samtools index "$OUT"

    echo "${S},${SEX},0,${S},${OUT},${OUT}.bai" >> "$CSV"
done

echo
echo "Samplesheet written to: ${CSV}"
echo
echo "Run the DeepVariant joint-genotyping pipeline (GRCh37) with, e.g.:"
echo "  nextflow run . -profile docker \\"
echo "    --input ${CSV} \\"
echo "    --step variant_calling \\"
echo "    --tools deepvariant --joint_genotype \\"
echo "    --genome GATK.GRCh37 \\"
echo "    --outdir results_1000g_joint"
echo
echo "Outputs (two distinct locations):"
echo "  - Raw joint (GLnexus) cohort VCF:"
echo "      results_1000g_joint/variant_calling/deepvariant/joint_variant_calling/"
echo "  - Final VEP-annotated cohort VCF (the deliverable):"
echo "      results_1000g_joint/annotation/deepvariant/joint_variant_calling/   (*_VEP.ann.vcf.gz)"
echo
echo "Verify it is multi-sample (3 samples):"
echo "  bcftools query -l results_1000g_joint/variant_calling/deepvariant/joint_variant_calling/*.vcf.gz"
echo "Verify VEP annotation is present on the final VCF:"
echo "  bcftools view -h results_1000g_joint/annotation/deepvariant/joint_variant_calling/*_VEP.ann.vcf.gz | grep 'ID=CSQ'"
