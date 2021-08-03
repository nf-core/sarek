#!/usr/bin/env bash
set -euo pipefail

# Build and push all containers

build_push() {
    GENOME=$1
    SPECIES=$2
    VEP_VERSION=$3
    VEP_TAG=$4

    docker build \
        -t nfcore/vep:${VEP_TAG}.${GENOME} \
        software/vep/. \
        --build-arg GENOME=${GENOME} \
        --build-arg SPECIES=${SPECIES} \
        --build-arg VEP_VERSION=${VEP_VERSION}

    docker push nfcore/vep:${VEP_TAG}.${GENOME}
}

build_push "GRCh37"    "homo_sapiens"           "104" "104.3"
build_push "GRCh38"    "homo_sapiens"           "104" "104.3"
build_push "GRCm38"    "mus_musculus"           "102" "104.3"
build_push "GRCm39"    "mus_musculus"           "104" "104.3"
build_push "CanFam3.1" "canis_lupus_familiaris" "104" "104.3"
build_push "WBcel235"  "caenorhabditis_elegans" "104" "104.3"
