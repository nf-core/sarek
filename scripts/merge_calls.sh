#!/bin/bash -x

TUMOR_SAMPLE=$1
NORMAL_SAMPLE=$2
SOURCE_DIR=$3

module load picard
module load htslib

# first make a list of files to merge 
mkdir MERGING
for f in `find ${SOURCE_DIR}/ -name "*_vs_*.vcf" -or -name "*_vs_*.vcf.gz" | egrep -v "candidate|diploid|freebayes"` ; do cp -v $f MERGING ; done
find ${SOURCE_DIR} -name "free*filtered.vcf" | xargs -I {} cp -v {} MERGING/
find MERGING -type f > filestomerge.list
# this is only to check which are those where we have to fix the NORMAL/TUMOR string
#for f in `cat filestomerge.list`; do echo $f;zgrep -m1 CHROM $f; done

# fix the string for MuTect2, Strelka SNPs and indels
MUTECT2=`grep mutect2 filestomerge.list`
STRELKA_SNVS=`grep snvs filestomerge.list`
STRELKA_INDELS=`grep indels filestomerge.list`
for f in $MUTECT2 $STRELKA_SNVS $STRELKA_INDELS; do 
	echo "Fixing $f"
	zcat $f | perl -p -e 's/NORMAL/'$NORMAL_SAMPLE'/' | perl -p -e 's/TUMOR/'$TUMOR_SAMPLE'/' > ${f%.gz}
	bgzip -f ${f%.gz}
	tabix $f
done

# double-check
echo "Now the lines should not contain NORMAL/TUMOR strings, but the corresponding sample names"
for f in `cat filestomerge.list`; do echo $f;zgrep -m1 CHROM $f; done


# next is remove the AD field and fix genotypes for FreeBayes
# this step actually can be done when filtering FB variants, but now we are over that step
FB_FILE=`grep filtered filestomerge.list`
awk -f removeADfromFB.awk $FB_FILE | bgzip -c > ${FB_FILE}.gz
tabix ${FB_FILE}.gz
# add .gz to the fixed gzipped file
perl -p -i -e 's/filtered.vcf/filtered.vcf.gz/' filestomerge.list

# This is the actual merging step
java -Xmx7g -jar $PICARD_HOME/picard.jar MergeVcfs I=filestomerge.list O=${TUMOR_SAMPLE}_${NORMAL_SAMPLE}_all_somatic_calls.vcf 
awk '/^#/||$7~/PASS/{print}' ${TUMOR_SAMPLE}_${NORMAL_SAMPLE}_all_somatic_calls.vcf > ${TUMOR_SAMPLE}_${NORMAL_SAMPLE}_filtered_somatic_calls.vcf

echo "Results written to ${TUMOR_SAMPLE}_${NORMAL_SAMPLE}_all_somatic_calls.vcf (all the calls) ${TUMOR_SAMPLE}_${NORMAL_SAMPLE}_filtered_somatic_calls.vcf (only the PASS-ed ones) files"
