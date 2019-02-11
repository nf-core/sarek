#!/bin/bash
set -xeuo pipefail

CODENAME=''
RELEASE=''

while [[ $# -gt 0 ]]
do
  key=$1
  case $key in
    -c|--codename)
    CODENAME=$2
    shift # past argument
    shift # past value
    ;;
    -r|--release)
    RELEASE=$2
    shift # past argument
    shift # past value
    ;;
  esac
done

if [[ $RELEASE == "" ]]
then
  echo "No release specified"
  exit
fi

if [[ $CODENAME == "" ]]
then
  echo "Preparing release $RELEASE"
  sed -i "s/\[Unreleased\]/[$RELEASE] - $(date +'%Y-%m-%d')/g" CHANGELOG.md
else
  echo "Preparing release $RELEASE - $CODENAME"
  sed -i "s/\[Unreleased\]/[$RELEASE] - $CODENAME - $(date +'%Y-%m-%d')/g" CHANGELOG.md
fi

sed -i "s/sarek-[0-9\.]\+/sarek-$RELEASE/g" Dockerfile
sed -i "s/sarek-[0-9\.]\+/sarek-$RELEASE/g" environment.yml
sed -i "s/sarek-snpeff-[0-9\.]\+/sarek-snpeff-$RELEASE/g" containers/snpeffgrch37/Dockerfile
sed -i "s/sarek-snpeff-[0-9\.]\+/sarek-snpeff-$RELEASE/g" containers/snpeffgrch37/environment.yml
sed -i "s/sarek-snpeff-[0-9\.]\+/sarek-snpeff-$RELEASE/g" containers/snpeffgrch38/Dockerfile
sed -i "s/sarek-snpeff-[0-9\.]\+/sarek-snpeff-$RELEASE/g" containers/snpeffgrch38/environment.yml
sed -i "s/sarek-vep-[0-9\.]\+/sarek-vep-$RELEASE/g" containers/vepgrch37/Dockerfile
sed -i "s/sarek-vep-[0-9\.]\+/sarek-vep-$RELEASE/g" containers/vepgrch37/environment.yml
sed -i "s/sarek-vep-[0-9\.]\+/sarek-vep-$RELEASE/g" containers/vepgrch38/Dockerfile
sed -i "s/sarek-vep-[0-9\.]\+/sarek-vep-$RELEASE/g" containers/vepgrch38/environment.yml
sed -i "s/sarek-[0-9\.]\+/sarek-$RELEASE/g" Singularity
sed -i "s/VERSION [0-9\.]\+/VERSION $RELEASE/g" Singularity
sed -i "s/version = '[0-9\.]\+'/version = '$RELEASE'/g" nextflow.config

git commit CHANGELOG.md Dockerfile environment.yml Singularity nextflow.config -m "preparing release $RELEASE [skip ci]"
