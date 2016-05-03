#!/usr/bin/env nextflow

/*
 * Sample run data for Variant Calling
 * use like (on milou):
 * ./nextflow run variantCalling.nf --tumor_bam ~/dev/chr17_testdata/HCC1143.tumor.bam --normal_bam ~/dev/chr17_testdata/HCC1143.normal.bam
 */

// First checking the existence of the BAM file and its index for the tumor
params.tumor_bam = "tumor.bam"
tumor_bam = file(params.tumor_bam)
if(!tumor_bam.exists()) exit 1, "Missing tumor file ${tumor_bam}; please specify --tumor_bam <SAMPLE> --normal_bam <SAMPLE> " 
// the index
params.tumor_bai = params.tumor_bam.replaceFirst(/.bam/,".bam.bai")
tumor_bai = file(params.tumor_bai)
if(!tumor_bai.exists()) exit 1, "Missing tumor file index ${tumor_bai}" 

// Ditto for the normal
params.normal_bam = "normal.bam"
normal_bam = file(params.normal_bam)
if(!normal_bam.exists()) exit 1, "Missing normal file ${normal_bam}; please specify --tumor_bam <SAMPLE> --normal_bam <SAMPLE>"  
// the normal index
params.normal_bai = params.normal_bam.replaceFirst(/.bam/,".bam.bai")
normal_bai = file(params.normal_bai)
if(!normal_bai.exists()) exit 1, "Missing normal file index ${normal_bai}"
 
params.genome		= "/sw/data/uppnex/ToolBox/ReferenceAssemblies/hg38make/bundle/2.8/b37/human_g1k_v37_decoy.fasta"
params.genomeidx	= "${params.genome}.fai"
params.cosmic		= "/sw/data/uppnex/ToolBox/ReferenceAssemblies/hg38make/bundle/2.8/b37/b37_cosmic_v54_120711.vcf"
params.cosmicidx	= "/sw/data/uppnex/ToolBox/ReferenceAssemblies/hg38make/bundle/2.8/b37/b37_cosmic_v54_120711.vcf.idx"
params.dbsnp		= "/sw/data/uppnex/ToolBox/ReferenceAssemblies/hg38make/bundle/2.8/b37/dbsnp_138.b37.vcf"
params.dbsnpidx		= "/sw/data/uppnex/ToolBox/ReferenceAssemblies/hg38make/bundle/2.8/b37/dbsnp_138.b37.vcf.idx"
params.mutect1_home	= "/sw/apps/bioinfo/mutect/1.1.5/milou"

genome_file		= file(params.genome)
genome_index	= file(params.genomeidx)
cosmic			= file(params.cosmic)
cosmicidx		= file(params.cosmicidx)
dbsnp			= file(params.dbsnp)
dbsnpidx		= file(params.dbsnpidx)

process MuTect1 {
	"""
	nextflow run $baseDir/MuTect1.nf \
	-w ~/dev/CAW/work \
	--tumor_bam ${params.tumor_bam} \
	--normal_bam ${params.normal_bam} \
	--genome ${params.genome} \
	--genomeidx ${params.genomeidx} \
	--cosmic ${params.cosmic} \
	--cosmicidx ${params.cosmicidx} \
	--dbsnp ${params.dbsnp} \
	--dbsnpidx ${params.dbsnpidx} \
	--mutect1_home ${params.mutect1_home}
	"""
}

// process FreeBayes {
// 	"""
// 	nextflow run $baseDir/FreeBayes.nf \
// 	-w ~/workspace/nextflow/work \
// 	--tumor_bam params.tumor_bam
// 	--normal_bam params.normal_bam
// 	"""
// }

// process VarDictJava {
// 	"""
// 	nextflow run $baseDir/VarDictJava.nf \
// 	-w ~/workspace/nextflow/work \
// 	--tumor_bam params.tumor_bam
// 	--normal_bam params.normal_bam
// 	"""
// }
