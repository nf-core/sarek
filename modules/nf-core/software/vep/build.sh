export GENOME=GRCh37
export SPECIES=homo_sapiens
export VEP_VERSION=104
export VEP_TAG=104.3

docker build -t nfcore/vep:${VEP_TAG}.${GENOME} modules/nf-core/software/vep/. --build-arg GENOME=${GENOME} --build-arg SPECIES=${SPECIES} --build-arg VEP_VERSION=${VEP_VERSION}
docker push nfcore/vep:${VEP_TAG}.${GENOME}

export GENOME=GRCh38
export SPECIES=homo_sapiens
export VEP_VERSION=104
export VEP_TAG=104.3

docker build -t nfcore/vep:${VEP_TAG}.${GENOME} modules/nf-core/software/vep/. --build-arg GENOME=${GENOME} --build-arg SPECIES=${SPECIES} --build-arg VEP_VERSION=${VEP_VERSION}
docker push nfcore/vep:${VEP_TAG}.${GENOME}

export GENOME=GRCm39
export SPECIES=mus_musculus
export VEP_VERSION=104
export VEP_TAG=104.3

docker build -t nfcore/vep:${VEP_TAG}.${GENOME} modules/nf-core/software/vep/. --build-arg GENOME=${GENOME} --build-arg SPECIES=${SPECIES} --build-arg VEP_VERSION=${VEP_VERSION}
docker push nfcore/vep:${VEP_TAG}.${GENOME}

export GENOME=GRCm38
export SPECIES=mus_musculus
export VEP_VERSION=102
export VEP_TAG=104.3

docker build -t nfcore/vep:${VEP_TAG}.${GENOME} modules/nf-core/software/vep/. --build-arg GENOME=${GENOME} --build-arg SPECIES=${SPECIES} --build-arg VEP_VERSION=${VEP_VERSION}
docker push nfcore/vep:${VEP_TAG}.${GENOME}

export GENOME=CanFam3.1
export SPECIES=canis_lupus_familiaris
export VEP_VERSION=104
export VEP_TAG=104.3

docker build -t nfcore/vep:${VEP_TAG}.${GENOME} modules/nf-core/software/vep/. --build-arg GENOME=${GENOME} --build-arg SPECIES=${SPECIES} --build-arg VEP_VERSION=${VEP_VERSION}
docker push nfcore/vep:${VEP_TAG}.${GENOME}

export GENOME=WBcel235
export SPECIES=caenorhabditis_elegans
export VEP_VERSION=104
export VEP_TAG=104.3

docker build -t nfcore/vep:${VEP_TAG}.${GENOME} modules/nf-core/software/vep/. --build-arg GENOME=${GENOME} --build-arg SPECIES=${SPECIES} --build-arg VEP_VERSION=${VEP_VERSION}
docker push nfcore/vep:${VEP_TAG}.${GENOME}
