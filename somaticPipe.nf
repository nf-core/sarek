#!/usr/bin/env nextflow
 
/*
 * Defines pipeline parameters in order to specify the refence genomes
 * and read pairs by using the command line options
 */

params.sample = "tcga.cl" // override with --sample <SAMPLE>
params.tpair1 = "data/${params.sample}.tumor_R1.fastq.gz"
params.npair1 = "data/${params.sample}.normal_R1.fastq.gz"
params.tpair2 = "data/${params.sample}.tumor_R2.fastq.gz"
params.npair2 = "data/${params.sample}.normal_R2.fastq.gz"
params.genome = "/proj/b2011196/nobackup/data/reference/GATK/bundle_2_8/b37/human_g1k_v37.fasta"
params.genomeidx = "${params.genome}.fai"
params.out = "$PWD"
//params.genome = "/sw/data/uppnex/ToolBox/ReferenceAssemblies/hg38make/bundle/2.8/b37/human_g1k_v37_decoy.fasta"

  

/*
 * validate input
 */


genome_file = file(params.genome)
genome_index = file(params.genomeidx)
tp1 = file(params.tpair1)
tp2 = file(params.tpair2)
np1 = file(params.npair1)
np2 = file(params.npair2)


if( !genome_file.exists() ) exit 1, "Missing reference: ${genome_file}"
if( !genome_index.exists() ) exit 1, "Missing index: ${genome_index}"
if( !tp1.exists() ) exit 2, "Missing read ${tp1}"
if( !tp2.exists() ) exit 2, "Missing read ${tp2}"
if( !np1.exists() ) exit 2, "Missing read ${np1}"
if( !np2.exists() ) exit 2, "Missing read ${np2}"


/*
 * processes, would be nice to refactor to have same process for tumor normal
 *
 * maybe create a channel that spits out pairs of tumor/normal fastq pairs
 * but Im not figuring out how to right now, so stick with full 4 files (tumor/normal pair)
 */	

process mapping_tumor_bwa {

	module 'bioinfo-tools'
	module 'bwa'
	module 'samtools/1.3'


	cpus 1

	input:
	file genome_file
	
	file tp1
	file tp2

	output:
	file '*.bam' into tumor_bam

	"""
	bwa mem -R "@RG\\tID:${params.sample}.tumor\\tSM:${params.sample}\\tLB:${params.sample}.tumor\\tPL:illumina" -B 3 -t ${task.cpus} -M ${params.genome} ${tp1} ${tp2} | samtools view -bS -t ${genome_index} - | samtools sort - > ${params.sample}.tumor.bam
	"""	

}

process mapping_normal_bwa {

	module 'bioinfo-tools'
	module 'bwa'
	module 'samtools/1.3'

	cpus 1


	input:
	file genome_file
	file np1
	file np2

	output:
	file '*.bam' into normal_bam // can we now use bam in the following?

	"""
	bwa mem -R "@RG\\tID:${params.sample}.normal\\tSM:${params.sample}\\tLB:${params.sample}.normal\\tPL:illumina" -B 3 -t ${task.cpus} -M ${params.genome} ${np1} ${np2} | samtools view -bS -t ${genome_index} - | samtools sort - > ${params.sample}.normal.bam
	"""	

}

/*
 * mark duplicates
 */



process mark_duplicates_tumor {

	module 'bioinfo-tools'
	module 'picard'


	input:
	file tumor_bam
	
	output:
	file '*md.bam' into tumor_md_bam


	"""
	java -Xmx7g -jar /sw/apps/bioinfo/picard/1.118/milou/MarkDuplicates.jar INPUT=${tumor_bam} METRICS_FILE=${tumor_bam}.metrics TMP_DIR=. ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=TRUE OUTPUT=${params.sample}.tumor.md.bam	
	"""



}

process mark_duplicates_normal {

	module 'bioinfo-tools'
	module 'picard'


	input:
	file normal_bam
	
	output:
	file '*md.bam' into normal_md_bam


	"""
	java -Xmx7g -jar /sw/apps/bioinfo/picard/1.118/milou/MarkDuplicates.jar INPUT=${normal_bam} METRICS_FILE=${normal_bam}.metrics TMP_DIR=. ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=TRUE OUTPUT=${params.sample}.normal.md.bam	
	"""

}

gatk_bundle  = "/sw/data/uppnex/ToolBox/ReferenceAssemblies/hg38make/bundle/2.8/b37"
1000g_indels = file("${gatk_bundle}/1000G_phase1.indels.b37.vcf")
dbsnp        = file("${gatk_bundle}/dbsnp_138.b37.vcf")
mills_indels = file("${gatk_bundle}/Mills_and_1000G_gold_standard.indels.b37.vcf")



/*
 * here we may need to merge tumor/normal???
 */

process make_intervals_tumor {


	input:
	file tumor_md_bam
	file normal_md_bam
	file genome_file
	file 1000g_indels
	file mills_indels

	output:
	file '*.intervals' into intervals_tumor



	"""
	java -Xmx7g -jar /sw/apps/bioinfo/GATK/3.3.0/GenomeAnalysisTK.jar -T RealignerTargetCreator -I ${tumor_md_bam} -I ${normal_md_bam} -R ${genome_file} -known ${1000g_indels} -known ${mills_indels} -o ${params.sample}.intervals
	"""

}



// realigntarget

/*
java -jar $GATK_HOME/GenomeAnalysisTK.jar \
-T RealignerTargetCreator \
-R $REFERENCE \
-I $NBAM -I $TBAM \
-known $KGINDELS \
-known $MILLS \
-o $SPAIR".intervals" \
-nt 16
if test $? = 0;then
    echo "RealignerTargetCreator went well!"
fi

*/
// realign
// fixmate
// recal


/*

GATK_HOME:=/sw/apps/bioinfo/GATK/3.3.0
GATK_BUNDLE:=/sw/data/uppnex/ToolBox/ReferenceAssemblies/hg38make/bundle/2.8/b37
DBSNP:=$(GATK_BUNDLE)/dbsnp_138.b37.vcf
1000_G_P1_INDELS:=$(GATK_BUNDLE)/1000G_phase1.indels.b37.vcf
MILLS_1000G_GS_INDELS:=$(GATK_BUNDLE)/Mills_and_1000G_gold_standard.indels.b37.vcf

.PRECIOUS: %.sorted.md.real.fm.recal.table

# realign, first create intervals
%.intervals : %.bam %.bam.bai
        java -Xmx$(MEM)g -jar $(GATK_HOME)/GenomeAnalysisTK.jar -T RealignerTargetCreator -I $< -R $(REF_FASTA) -known $(1000_G_P1_INDELS) -known $(MILLS_1000G_GS_INDELS) -o $@.tmp && mv $@.tmp $@

%.real.bam : %.bam %.intervals %.bam.bai
        java -Xmx$(MEM)g -Djava.io.tmpdir=$$SNIC_TMP -jar $(GATK_HOME)/GenomeAnalysisTK.jar -T IndelRealigner -I $< -R $(REF_FASTA) -known $(1000_G_P1_INDELS) -known $(MILLS_1000G_GS_INDELS) -targetIntervals $*.intervals -o $@.tmp && mv $@.tmp $@


# recalibrate bases, gatk style:
# first create recal table:
%.recal.table : %.bam %.bam.bai
        java -Xmx$(MEM)g -Djava.io.tmpdir=$$SNIC_TMP -jar $(GATK_HOME)/GenomeAnalysisTK.jar -l INFO -R $(REF_FASTA) -I $< -T BaseRecalibrator -knownSites $(DBSNP) -knownSites $(1000_G_P1_INDELS) -knownSites $(MILLS_1000G_GS_INDELS) -o $@.tmp && mv $@.tmp $@

# then recal bam file:
%.recal.bam : %.bam %.recal.table %.bam.bai
        java -Xmx$(MEM)g -Djava.io.tmpdir=$$SNIC_TMP -jar $(GATK_HOME)/GenomeAnalysisTK.jar -R $(REF_FASTA) -I $< -T PrintReads --BQSR $*.recal.table -o $@.tmp && mv $@.tmp $@ && mv $@.tmp.bai $@.bai


# make gvcf:
%.g.vcf : %.bam
        java -Xmx$(MEM)g -Djava.io.tmpdir=$$SNIC_TMP -jar $(GATK_HOME)/GenomeAnalysisTK.jar -R $(REF_FASTA) -T HaplotypeCaller -I $< --emitRefConfidence GVCF --variant_index_type LINEAR --dbsnp $(DBSNP) --variant_index_parameter 128000 -nct $(THREADS) -o $@.tmp && mv $@.tmp $@ && mv $@.tmp.idx $@.idx

%.g.vcf.gz : %.bam
        java -Xmx$(MEM)g -Djava.io.tmpdir=$$SNIC_TMP -jar $(GATK_HOME)/GenomeAnalysisTK.jar -R $(REF_FASTA) -T HaplotypeCaller -I $< --emitRefConfidence GVCF --variant_index_type LINEAR --dbsnp $(DBSNP) --variant_index_parameter 128000 -nct $(THREADS) -o $@.tmp.gz && mv $@.tmp.gz $@ && mv $@.tmp.tbi $@.tbi


*/
