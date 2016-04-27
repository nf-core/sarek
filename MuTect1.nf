#!/usr/bin/env nextflow
/*
 *      Sample run data for VarDictJava
 * use like (on milou):
 * ./nextflow run VarDictJava.nf --tumor_bam ~/dev/chr17_testdata/HCC1143.tumor.bam --normal_bam ~/dev/chr17_testdata/HCC1143.normal.bam
 */

//
// First checking the existence of the BAM file and its index for the tumor
params.tumor_bam = "tumor.bam" // override with --tumor_bam <SAMPLE>
tumor_bam = file(params.tumor_bam)
if(!tumor_bam.exists()) exit 1, "Missing tumor file ${tumor_bam}" 
// the index
params.tumor_bai = params.tumor_bam.replaceFirst(/.bam/,".bam.bai")
tumor_bai = file(params.tumor_bai)
if(!tumor_bai.exists()) exit 1, "Missing tumor file index ${tumor_bai}" 

// Ditto for the normal
params.normal_bam = "normal.bam" // override with --normal_bam <SAMPLE>
normal_bam = file(params.normal_bam)
if(!normal_bam.exists()) exit 1, "Missing normal file ${normal_bam}" 
// the normal index
params.normal_bai = params.normal_bam.replaceFirst(/.bam/,".bam.bai")
normal_bai = file(params.normal_bai)
if(!normal_bai.exists()) exit 1, "Missing normal file index ${normal_bai}"
 
params.genome = "/sw/data/uppnex/ToolBox/ReferenceAssemblies/hg38make/bundle/2.8/b37/human_g1k_v37_decoy.fasta"
params.genomeidx = "${params.genome}.fai"
params.cosmic = "/sw/data/uppnex/ToolBox/ReferenceAssemblies/hg38make/bundle/2.8/b37/b37_cosmic_v54_120711.vcf"
params.cosmicidx = "/sw/data/uppnex/ToolBox/ReferenceAssemblies/hg38make/bundle/2.8/b37/b37_cosmic_v54_120711.vcf.idx"
params.dbsnp = "/sw/data/uppnex/ToolBox/ReferenceAssemblies/hg38make/bundle/2.8/b37/dbsnp_138.b37.vcf"
params.dbsnpidx = "/sw/data/uppnex/ToolBox/ReferenceAssemblies/hg38make/bundle/2.8/b37/dbsnp_138.b37.vcf.idx"

genome_file = file(params.genome)
genome_index = file(params.genomeidx)
cosmic = file(params.cosmic)
cosmicidx = file(params.cosmicidx)
dbsnp = file(params.dbsnp)
dbsnpidx = file(params.dbsnpidx)

process MuTect1 {
    module 'bioinfo-tools'
    module 'mutect'

    println $MUTECT_HOME


    cpus 2

    input:
    file genome_file
    file genome_index
    
    file tumor_bam
    file tumor_bai
    file normal_bam
    file normal_bai

    output:
    file '*.VarDict.vcf' into vardict_vc_call

    """
    VarDict -G ${params.genome} -f 0.01 -N TSN -b "${normal_bam}|${tumor_bam}" -z 1 -F 0x500  -c 1 -S 2 -E 3 -g 4 -R chr17:1000000-1100000 | testsomatic.R | var2vcf_somatic.pl > test.VarDict.vcf
    """
}
