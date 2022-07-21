#!/usr/bin/env bash
set -euo pipefail

# Build and push all containers

build_push() {
    GENOME=$1
    SPECIES=$2
    VEP_CACHE_VERSION=$3
    VEP_VERSION=$4

    docker build \
        . \
        -t nfcore/vep:${VEP_VERSION}.${GENOME} \
        --build-arg GENOME=${GENOME} \
        --build-arg SPECIES=${SPECIES} \
        --build-arg VEP_CACHE_VERSION=${VEP_CACHE_VERSION} \
        --build-arg VEP_VERSION=${VEP_VERSION}

    docker push nfcore/vep:${VEP_VERSION}.${GENOME}
}

build_push "GRCh37"    "homo_sapiens"           "106" "106.1"
build_push "GRCh38"    "homo_sapiens"           "106" "106.1"
build_push "GRCm38"    "mus_musculus"           "102" "106.1"
build_push "GRCm39"    "mus_musculus"           "106" "106.1"
build_push "CanFam3.1" "canis_lupus_familiaris" "104" "106.1"
build_push "WBcel235"  "caenorhabditis_elegans" "106" "106.1"
