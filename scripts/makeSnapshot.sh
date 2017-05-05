#!/bin/bash -x
TAGSTR=`git describe --tags`
echo $TAGSTR
PWD=`pwd`
CURRENTDIR=${PWD##*/}
tar czv \
--exclude=.git* \
--exclude=.nextflow* \
--exclude=*swp \
--exclude=doc/SprintReview* \
--exclude=Preprocessing* \
--exclude=References* \
--exclude=Reports* \
--exclude=timeline* \
--exclude=trace* \
--exclude=VariantCalling* \
--exclude=work* \
-f ../CAW.${TAGSTR}.tgz ../${CURRENTDIR}
