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
//params.genome = "/proj/b2011196/nobackup/data/reference/GATK/bundle_2_8/b37/human_g1k_v37.fasta"
params.genome = "/sw/data/uppnex/ToolBox/ReferenceAssemblies/hg38make/bundle/2.8/b37/human_g1k_v37_decoy.fasta"
params.genomeidx = "${params.genome}.fai"
params.out = "$PWD"
params.kgindels = "/sw/data/uppnex/ToolBox/ReferenceAssemblies/hg38make/bundle/2.8/b37/1000G_phase1.indels.b37.vcf"
params.dbsnp = "/sw/data/uppnex/ToolBox/ReferenceAssemblies/hg38make/bundle/2.8/b37/dbsnp_138.b37.vcf"
params.millsindels = "/sw/data/uppnex/ToolBox/ReferenceAssemblies/hg38make/bundle/2.8/b37/Mills_and_1000G_gold_standard.indels.b37.vcf"



  

/*
 * validate input
 */


genome_file = file(params.genome)
genome_index = file(params.genomeidx)
kgindels = file(params.kgindels)
dbsnp = file(params.dbsnp)
millsindels = file(params.millsindels)

tp1 = file(params.tpair1)
tp2 = file(params.tpair2)
np1 = file(params.npair1)
np2 = file(params.npair2)


if( !genome_file.exists() ) exit 1, "Missing reference: ${genome_file}"
if( !genome_index.exists() ) exit 1, "Missing index: ${genome_index}"
if( !kgindels.exists() ) exit 1, "Missing index: ${kgindels}"
if( !dbsnp.exists() ) exit 1, "Missing index: ${dbsnp}"
if( !millsindels.exists() ) exit 1, "Missing index: ${millsindels}"
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
	file '*.tumor.bam' into tumor_bam

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
	file '*.normal.bam' into normal_bam // can we now use bam in the following?

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
	file '*.tumor.md.bam' into tumor_md_bam


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
	file '*.normal.md.bam' into normal_md_bam


	"""
	java -Xmx7g -jar /sw/apps/bioinfo/picard/1.118/milou/MarkDuplicates.jar INPUT=${normal_bam} METRICS_FILE=${normal_bam}.metrics TMP_DIR=. ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=TRUE OUTPUT=${params.sample}.normal.md.bam	
	"""

}



/*
 * // realigntarget
 * here we may need to merge tumor/normal???
 */

process create_intervals {

	// cpus = 16

	input:
	file tumor_md_bam
	file normal_md_bam
	file genome_file
	file genome_idx
	file kgindels
	file millsindels

	output:
	file '*.intervals' into intervals



	"""
	java -Xmx7g -jar /sw/apps/bioinfo/GATK/3.3.0/GenomeAnalysisTK.jar -T RealignerTargetCreator -I ${tumor_md_bam} -I ${normal_md_bam} -R ${genome_file} -known ${kgindels} -known ${millsindels} -o ${params.sample}.intervals
	"""

}

/*
 * realign
 *
 * here we need to split into tumor/normal again
 */

process realign {


	input:
	file tumor_md_bam
	file normal_md_bam
	file genome_file
	file kgindels
	file millsindels
	file genome_idx
	file intervals 

	output:
	file '*.tumor.real.bam' into real_tumor_bam
	file '*.normal.real.bam' into real_normal_bam


	"""
	java -Xmx7g -jar /sw/apps/bioinfo/GATK/3.3.0/GenomeAnalysisTK.jar -T IndelRealigner -I ${tumor_md_bam} -I ${normal_md_bam} -R ${genome_file} -targetIntervals ${intervals} -known ${kgindels} -known ${millsindels} -nWayOut '.real.bam'
	"""

}


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
