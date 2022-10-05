#!/usr/bin/env bash
set -euo pipefail

# Build and push all containers

build_push() {
    GENOME=$1
    SNPEFF_CACHE_VERSION=$2
    SNPEFF_VERSION=$3

    docker build \
        . \
        -t nfcore/snpeff:${SNPEFF_VERSION}.${GENOME} \
        --build-arg GENOME=${GENOME} \
        --build-arg SNPEFF_CACHE_VERSION=${SNPEFF_CACHE_VERSION} \
        --build-arg SNPEFF_VERSION=${SNPEFF_VERSION}

    docker push nfcore/snpeff:${SNPEFF_VERSION}.${GENOME}
}

build_push "GRCh37"    "87"  "5.1"
build_push "GRCh38"    "105" "5.1"
build_push "GRCm38"    "99"  "5.1"
build_push "GRCm39"    "105" "5.1"
build_push "CanFam3.1" "99"  "5.1"
build_push "WBcel235"  "105" "5.1"
