#!/bin/bash -x
TAGSTR=`git describe --tags`
echo $TAGSTR
PWD=`pwd`
CURRENTDIR=${PWD##*/}
tar czv --exclude=.git* --exclude=.next* --exclude=work* --exclude=Preproc* --exclude=Variant* --exclude=Reports* --exclude=timeline* --exclude=trace* --exclude=*swp -f ../CAW.${TAGSTR}.tgz ../${CURRENTDIR}
