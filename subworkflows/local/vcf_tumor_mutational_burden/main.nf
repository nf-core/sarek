include { BCFTOOLS_NORM         } from '../../../modules/nf-core/bcftools/norm/main'
include { TMB                   } from '../../../modules/local/tmb/main'

workflow TUMOR_MUTATIONAL_BURDEN {

    take:
    vcf         // channel: [ val(meta), vcf.gz, vcf.gz.tbi ]
    fasta
    target_bed

    main:

    versions = Channel.empty()

    vcf_in = vcf.map{ meta, vcf, tbi -> [
        meta + [ dbconfig: getTMBconfig(meta.annotation),
                varconfig: getTMBconfig(meta.variantcaller) ],
        vcf, tbi ] }

    BCFTOOLS_NORM(vcf_in, fasta)


    tmb_in = BCFTOOLS_NORM.out.vcf.map{ meta, vcf -> [ meta, vcf, meta.dbconfig, meta.varconfig ] }

    TMB(tmb_in, target_bed)

    versions = versions.mix(BCFTOOLS_NORM.out.versions)
    versions = versions.mix(TMB.out.versions)

    emit:
    versions
}

//
// Function to retrieve config file for tumor mutational burden
//
def getTMBconfig(tool) {
    def configFile  = "$projectDir/assets/tmb/${tool}.yml"
    def file = new File(configFile)
    return file.exists() ? configFile : null
}
