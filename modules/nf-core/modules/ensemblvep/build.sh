#!/usr/bin/env bash
set -euo pipefail

# Build and push all containers

build_push() {
    GENOME=$1
    SPECIES=$2
    VEP_VERSION=$3
    VEP_TAG=$4

    docker build \
        . \
        -t nfcore/vep:${VEP_TAG}.${GENOME} \
        --build-arg GENOME=${GENOME} \
        --build-arg SPECIES=${SPECIES} \
        --build-arg VEP_VERSION=${VEP_VERSION} \
        --build-arg VEP_TAG=${VEP_TAG}

    docker push nfcore/vep:${VEP_TAG}.${GENOME}
}

build_push "GRCh37"    "homo_sapiens"           "105" "105.0"
build_push "GRCh38"    "homo_sapiens"           "105" "105.0"
build_push "GRCm38"    "mus_musculus"           "102" "105.0"
build_push "GRCm39"    "mus_musculus"           "105" "105.0"
build_push "CanFam3.1" "canis_lupus_familiaris" "104" "105.0"
build_push "WBcel235"  "caenorhabditis_elegans" "105" "105.0"
