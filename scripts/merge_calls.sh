#!/bin/bash

TUMOR_SAMPLE=$1
NORMAL_SAMPLE=$2
SOURCE_DIR=$3

module load GATK
module load tabix

# first make a list of files to merge 
find ${SOURCE_DIR}/ -name "*_vs_*.vcf" -or -name "*_vs_*.vcf.gz" |\
		egrep -v "candidate|diploid" > filestomerge.list

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
java -Xmx7g -jar $GATK_HOME/GenomeAnalysisTK.jar MergeVcfs -I filestomerge.list -O all_somatic_calls.vcf 
awk '/^#/||$7~/PASS/{print}' all_somatic_calls.vcf > filtered_somatic_calls.vcf

echo "Results written to all_somatic_calls.vcf (all the calls) filtered_somatic_calls.vcf (only the PASS-ed ones) files"
