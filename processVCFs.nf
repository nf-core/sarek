#!/usr/bin/env nextflow

/*
vim: syntax=groovy
-*- mode: groovy;-*-
================================================================================
=           combining and rough filtering VCF results of CAW                   =
================================================================================
*/

if( !params.vcflist) {
    help_message()
    exit 1, 'Missing list of VCFs to merge'
}
if( !params.callName) {
    help_message()
    exit 1, 'Missing call name'
}

// make a channel from command-line entries
vcfMappingChannel = makeVCFList(params.vcflist)

// In this process we are adding a tag to each VCFs, so after the merge process we can see which caller called what.
// Also applying a rough filter, leaving only records containing "PASS" in the 7th field in the VCF
process addCallerInfo {
    input:
    set callerTag, file(rawVCF) from vcfMappingChannel
    output:
    file "tagged*vcf" into taggedVCF

    script:
    """
    awk -f ${workflow.projectDir}/scripts/addCallerINFO.awk -v caller=${callerTag} ${rawVCF} > tagged${rawVCF}
    """
}

// we will need tabix, bcftools and vcftools to do the job
process compressAndIndexFilteredVCFs {

    module 'htslib/1.3'
    module 'bcftools/1.2'

    input:
    set file(vcfToCompress) from taggedVCF

    output:
    file "tagged*vcf.gz" into compressedVCF
    file "tagged*vcf.gz.csi" into csiVCF

    script:
    """
    bgzip ${vcfToCompress}
    bcftools index ${vcfToCompress}.gz
    """
    
}

// we will need the names of the compressed VCFs in two channels: one to make links, and the other to the command line
allVCFs = compressedVCF.reduce {a,b -> return a + " " + b}
process mergeVCFs {
    publishDir "mergedVariants",mode: 'copy'

    module 'htslib/1.3'
    module 'vcftools/0.1.14'

    input:
    set cVCFs from allVCFs

    output:
    file "${params.callName}_merged.vcf" into mergedVCF

    script:
    """
    vcf-merge ${cVCFs} > ${params.callName}_merged.vcf
    """
}

/*
    We are running the script like

    nextflow run processVCFs.nf --vcflist Manta=here/manta.vcf,MuTect1=there/mutect1.vcf    

    so, have an input list like "Manta=here/manta.vcf,MuTect1=there/mutect1.vcf" and we want to have a channel like
    [Manta, here/manta.vcf]
    [MuTect1, there/mutect1.vcf]
*/
def makeVCFList(vcflist) {
    // split the line at commas, and feed the list to a channel
    entries = vcflist.split(",")
    callChannel = Channel.from(entries)
    // each entry in the channel looks like "tag=file.vcf" . Now split at the = sign, save the tag and the filename
    vcfMappings = callChannel.map { line -> 
                        list    = line.split("=") 
                        tag     = list[0]
                        vcfFile = file(list[1])
                        [tag,vcfFile]
                    }
    
    return vcfMappings
}

def help_message() { // Display help message
    log.info "CANCER ANALYSIS WORKFLOW ~ processing VCFs "
    log.info "  Usage:"
    log.info "          nextflow run processVCF.nf -profile localhost --vcflist CALLER1=caller1.vcf,CALLER2=caller2.vcf[,...] --callName P12345 "
    log.info "  Example:"
    log.info "          nextflow run processVCF.nf -profile localhost --vcflist STRELKA=VariantCalling/Strelka/strelka/results/passed.somatic.snvs.vcf,MUTECT1=VariantCalling/MuTect1/results.vcf --callName P1234"
}
