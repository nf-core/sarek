#! /bin/bash

#Input data
vcf=$1
for chr in `seq 1 22` X Y M
do
echo "chr "$chr
echo grep "$chr" $vcf
done
