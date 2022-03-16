#!/usr/bin/env bash
set -euo pipefail

# Build and push all containers

build_push() {
    GENOME=$1
    SNPEFF_CACHE_VERSION=$2
    SNPEFF_TAG=$3

    docker build \
        -t nfcore/snpeff:${SNPEFF_TAG}.${GENOME} \
        software/snpeff/. \
        --build-arg GENOME=${GENOME} \
        --build-arg SNPEFF_CACHE_VERSION=${SNPEFF_CACHE_VERSION}

    docker push nfcore/snpeff:${SNPEFF_TAG}.${GENOME}
}

build_push "GRCh37"    "75" "5.0"
build_push "GRCh38"    "99" "5.0"
build_push "GRCm38"    "99" "5.0"
build_push "CanFam3.1" "99" "5.0"
build_push "WBcel235"  "99" "5.0"
