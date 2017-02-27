#!/bin/bash -x
# Command line options should be:
# $1 string for sample name
#
# the script working directory have to be one above the VariantCalling directory
#
awk '$7 !~/REJECT/{print}' ./VariantCalling/MuTect1/MuTect1_*vcf > ./VariantCalling/MuTect1/MuTect1_NR.vcf
#VCFS="./VariantCalling/MuTect1/MuTect1_NR.vcf ./VariantCalling/MuTect2/MuTect2_*.vcf ./VariantCalling/Strelka/strelka/results/all.somatic.*.vcf ./VariantCalling/HaplotypeCaller/HaplotypeCaller_*_normal*.vcf ./VariantCalling/HaplotypeCaller/HaplotypeCaller_*_tumor*.vcf ./VariantCalling/Manta/*.candidateSmallIndels.vcf "

# this set is for testing
VCFS="./VariantCalling/MuTect1/MuTect1_NR.vcf ./VariantCalling/MuTect2/MuTect2_*.vcf ./VariantCalling/Strelka/strelka/results/all.somatic.*.vcf ./VariantCalling/Manta/*.candidateSmallIndels.vcf "

for f in $VCFS; do {
    echo -n "compressing $f ... "
    bgzip -i -@4 -f -c $f > ${f}.gz; 
    echo -n " indexing ... "
    bcftools index -f ${f}.gz 
    echo "done"
} done
VCFSTOMERGE=`for f in $VCFS; do echo -n ${f}.gz" "; done`
OUTVCF=union${1}.vcf

vcf-merge $VCFSTOMERGE > VariantCalling/$OUTVCF 2>/dev/null

# get an intersect, where we are reporting only those that are reported by both MuTects and Strelka
# but we are ignoring HaplotypeCaller and the less sensitive Manta indel caller
awk '/^#/ {print} !/^#/{ if($10 !~/^\.$/ && $12 !~/^\.$/ && ($15 !~/^\.$/ || $17 !~/^\.$/) ) print }' VariantCalling/${OUTVCF} > VariantCalling/intersect${1}.vcf

# get a set that dos not contain REJECT entries from MuTect1
awk '/^#/ {print} $7 !~/REJECT/ && !/^#/{print}' VariantCalling/${OUTVCF} > VariantCalling/nonreject${1}.vcf
echo "Union and intersection at VariantCalling/$OUTVCF and VariantCalling/intersect${1}.vcf, also union without REJECTS at VariantCalling/nonreject${1}.vcf"
echo "Calculating statistics"
for s in union intersect nonreject; do
    PREFIX=${s}${1}
    echo -n "${PREFIX} ... "
    bcftools stats VariantCalling/${PREFIX}.vcf > VariantCalling/${PREFIX}.stats
    plot-vcfstats -P -s -p VariantCalling/${PREFIX} VariantCalling/${PREFIX}.stats
    echo "done"
done

