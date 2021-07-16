export GENOME=GRCh37
export SNPEFF_CACHE_VERSION=75
export SNPEFF_TAG=5.0

docker build -t nfcore/snpeff:${SNPEFF_TAG}.${GENOME} modules/nf-core/software/snpeff/. --build-arg GENOME=${GENOME} --build-arg SNPEFF_CACHE_VERSION=${SNPEFF_CACHE_VERSION}
docker push nfcore/snpeff:${SNPEFF_TAG}.${GENOME}

export GENOME=GRCh38
export SNPEFF_CACHE_VERSION=99
export SNPEFF_TAG=5.0

docker build -t nfcore/snpeff:${SNPEFF_TAG}.${GENOME} modules/nf-core/software/snpeff/. --build-arg GENOME=${GENOME} --build-arg SNPEFF_CACHE_VERSION=${SNPEFF_CACHE_VERSION}
docker push nfcore/snpeff:${SNPEFF_TAG}.${GENOME}

export GENOME=GRCm38
export SNPEFF_VERSION=99
export SNPEFF_TAG=5.0

docker build -t nfcore/snpeff:${SNPEFF_TAG}.${GENOME} modules/nf-core/software/snpeff/. --build-arg GENOME=${GENOME} --build-arg SNPEFF_CACHE_VERSION=${SNPEFF_CACHE_VERSION}
docker push nfcore/snpeff:${SNPEFF_TAG}.${GENOME}

export GENOME=CanFam3.1
export SNPEFF_VERSION=99
export SNPEFF_TAG=5.0

docker build -t nfcore/snpeff:${SNPEFF_TAG}.${GENOME} modules/nf-core/software/snpeff/. --build-arg GENOME=${GENOME} --build-arg SNPEFF_CACHE_VERSION=${SNPEFF_CACHE_VERSION}
docker push nfcore/snpeff:${SNPEFF_TAG}.${GENOME}

export GENOME=WBcel235
export SNPEFF_VERSION=99
export SNPEFF_TAG=5.0

docker build -t nfcore/snpeff:${SNPEFF_TAG}.${GENOME} modules/nf-core/software/snpeff/. --build-arg GENOME=${GENOME} --build-arg SNPEFF_CACHE_VERSION=${SNPEFF_CACHE_VERSION}
docker push nfcore/snpeff:${SNPEFF_TAG}.${GENOME}
