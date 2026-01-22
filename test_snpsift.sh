#!/bin/bash

# Simple test script for SnpSift annotation
# This will annotate test.vcf.gz with test2.vcf.gz using SnpSift annotateMem

cd /Users/friederike.hanssen/Projects/sarek-snpsift

nextflow run main.nf \
    -profile test,docker \
    --input tests/csv/3.0/vcf_single.csv \
    --step annotate \
    --tools snpsift \
    --snpsift_databases tests/config/snpsift_test_databases.json \
    --snpsift_create_dbs true \
    --outdir results_snpsift_test \
    --validationSchemaIgnoreParams input,snpsift_databases
