include { BWA_INDEX as BWAMEM1_INDEX                } from '../../../modules/nf-core/bwa/index'
include { BWAMEM2_INDEX                             } from '../../../modules/nf-core/bwamem2/index'
include { DRAGMAP_HASHTABLE                         } from '../../../modules/nf-core/dragmap/hashtable'
include { GATK4_CREATESEQUENCEDICTIONARY            } from '../../../modules/nf-core/gatk4/createsequencedictionary'
include { MSISENSORPRO_SCAN                         } from '../../../modules/nf-core/msisensorpro/scan'
include { MSISENSOR2_SCAN                           } from '../../../modules/nf-core/msisensor2/scan'
include { SAMTOOLS_FAIDX                            } from '../../../modules/nf-core/samtools/faidx'
include { TABIX_TABIX as TABIX_BCFTOOLS_ANNOTATIONS } from '../../../modules/nf-core/tabix/tabix'
include { TABIX_TABIX as TABIX_DBSNP                } from '../../../modules/nf-core/tabix/tabix'
include { TABIX_TABIX as TABIX_GERMLINE_RESOURCE    } from '../../../modules/nf-core/tabix/tabix'
include { TABIX_TABIX as TABIX_KNOWN_INDELS         } from '../../../modules/nf-core/tabix/tabix'
include { TABIX_TABIX as TABIX_KNOWN_SNPS           } from '../../../modules/nf-core/tabix/tabix'
include { TABIX_TABIX as TABIX_PON                  } from '../../../modules/nf-core/tabix/tabix'
include { UNTAR as UNTAR_CHR_DIR                    } from '../../../modules/nf-core/untar'
include { UNTAR as UNTAR_MSISENSOR2_MODELS          } from '../../../modules/nf-core/untar'
include { UNZIP as UNZIP_ALLELES                    } from '../../../modules/nf-core/unzip'
include { UNZIP as UNZIP_GC                         } from '../../../modules/nf-core/unzip'
include { UNZIP as UNZIP_LOCI                       } from '../../../modules/nf-core/unzip'
include { UNZIP as UNZIP_RT                         } from '../../../modules/nf-core/unzip'

workflow PREPARE_GENOME {
    take:
    ascat_alleles_in            // params.ascat_alleles
    ascat_loci_in               // params.ascat_loci
    ascat_loci_gc_in            // params.ascat_loci_gc
    ascat_loci_rt_in            // params.ascat_loci_rt
    bcftools_annotations_in     // params.bcftools_annotations
    bcftools_annotations_tbi_in // params.bcftools_annotations
    bwa_in                      // params.bwa
    bwamem2_in                  // params.bwamem2
    chr_dir_in                  // params.chr_dir
    dbsnp_in                    // params.dbsnp
    dbsnp_tbi_in                // params.dbsnp_tbi
    dict_in                     // params.dict
    dragmap_in                  // params.dragmap
    fasta_in                    // params.fasta
    fasta_fai_in                // params.fasta_fai
    germline_resource_in        // params.germline_resource
    germline_resource_tbi_in    // params.germline_resource_tbi
    known_indels_in             // params.known_indels
    known_indels_tbi_in         // params.known_indels_tbi
    known_snps_in               // params.known_snps
    known_snps_tbi_in           // params.known_snps_tbi
    msisensor2_models_in        // channel: [optional]  msisensor2_models
    msisensor2_scan_in          // channel: [optional]  msisensor2_scan
    msisensorpro_scan_in        // channel: [optional]  msisensorpro_scan
    pon_in                      // params.pon
    pon_tbi_in                  // params.pon_tbi
    aligner                     // params.aligner
    step                        // params.step
    tools                       // params.tools
    vep_include_fasta           // params.vep_include_fasta

    main:
    versions = Channel.empty()

    // TODO: EXTRACT FASTA FILE?
    fasta = fasta_in ? Channel.fromPath(fasta_in).map { fasta -> [[id: fasta.baseName], fasta] }.collect() : Channel.empty()
    vep_fasta = vep_include_fasta ? fasta : [[id: 'null'], []]

    if (step == 'mapping') {
        if (!bwa_in && (aligner == "bwa-mem" || aligner == "sentieon-bwamem" || aligner == "parabricks")) {
            BWAMEM1_INDEX(fasta)
            index_alignment = BWAMEM1_INDEX.out.index.collect()
            versions = versions.mix(BWAMEM1_INDEX.out.versions)
        }
        else if (aligner == "bwa-mem" || aligner == "sentieon-bwamem" || aligner == "parabricks") {
            index_alignment = Channel.fromPath(bwa_in).map { index -> [[id: 'bwa'], index] }.collect()
        }
        else if (!bwamem2_in && aligner == 'bwa-mem2') {
            BWAMEM2_INDEX(fasta)
            index_alignment = BWAMEM2_INDEX.out.index.collect()
            versions = versions.mix(BWAMEM2_INDEX.out.versions)
        }
        else if (aligner == 'bwa-mem2') {
            index_alignment = Channel.fromPath(bwamem2_in).map { index -> [[id: 'bwamem2'], index] }.collect()
        }
        else if (!dragmap_in && aligner == 'dragmap') {
            DRAGMAP_HASHTABLE(fasta)
            index_alignment = DRAGMAP_HASHTABLE.out.hashmap.collect()
            versions = versions.mix(DRAGMAP_HASHTABLE.out.versions)
        }
        else if (aligner == 'dragmap') {
            index_alignment = Channel.fromPath(dragmap_in).map { index -> [[id: 'dragmap'], index] }.collect()
        }
    }
    else {
        index_alignment = Channel.empty()
    }

    if (!dict_in && step != "annotate") {
        GATK4_CREATESEQUENCEDICTIONARY(fasta)
        dict = GATK4_CREATESEQUENCEDICTIONARY.out.dict.collect()
        versions = versions.mix(GATK4_CREATESEQUENCEDICTIONARY.out.versions)
    }
    else if (dict_in) {
        dict = Channel.fromPath(dict_in).map { it -> [[id: 'dict'], it] }.collect()
    }
    else {
        dict = Channel.empty()
    }

    if (!fasta_fai_in && step != "annotate") {
        SAMTOOLS_FAIDX(fasta, [[id: 'no_fai'], []], false)
        fasta_fai = SAMTOOLS_FAIDX.out.fai.collect()
        versions = versions.mix(SAMTOOLS_FAIDX.out.versions)
    }
    else if (fasta_fai_in) {
        fasta_fai = Channel.fromPath(fasta_fai_in).map { it -> [[id: 'fai'], it] }.collect()
    }
    else {
        fasta_fai = Channel.empty()
    }

    bcftools_annotations = bcftools_annotations_in ? Channel.fromPath(bcftools_annotations_in).collect() : Channel.value([])
    bcftools_annotations_tbi = bcftools_annotations_tbi_in ? Channel.fromPath(bcftools_annotations_tbi_in).collect() : Channel.value([])

    if (!bcftools_annotations_tbi_in && bcftools_annotations_in) {
        TABIX_BCFTOOLS_ANNOTATIONS(bcftools_annotations.flatten().map { vcf -> [[id: vcf.baseName], vcf] })
        bcftools_annotations_tbi = TABIX_BCFTOOLS_ANNOTATIONS.out.tbi.map { _meta, tbi -> [tbi] }.collect()
        versions = versions.mix(TABIX_BCFTOOLS_ANNOTATIONS.out.versions)
    }

    dbsnp = dbsnp_in ? Channel.fromPath(dbsnp_in).collect() : Channel.value([])
    dbsnp_tbi = dbsnp_tbi_in ? Channel.fromPath(dbsnp_tbi_in).collect() : Channel.value([])

    if (!dbsnp_tbi_in && dbsnp_in && ((step == "mapping" || step == "markduplicates" || step == "prepare_recalibration") || (tools.split(',').contains('controlfreec') || tools.split(',').contains('haplotypecaller') || tools.split(',').contains('sentieon_haplotyper') || tools.split(',').contains('sentieon_dnascope') || tools.split(',').contains('muse') || tools.split(',').contains('mutect2')))) {
        TABIX_DBSNP(dbsnp.flatten().map { vcf -> [[id: vcf.baseName], vcf] })
        dbsnp_tbi = TABIX_DBSNP.out.tbi.map { _meta, tbi -> [tbi] }.collect()
        versions = versions.mix(TABIX_DBSNP.out.versions)
    }

    germline_resource = germline_resource_in ? Channel.fromPath(germline_resource_in).collect() : Channel.value([])
    germline_resource_tbi = germline_resource_tbi_in ? Channel.fromPath(germline_resource_tbi_in).collect() : Channel.value([])

    if (!germline_resource_tbi_in && germline_resource_in && (tools.split(',').contains('mutect2') || tools.split(',').contains('sentieon_tnscope'))) {
        TABIX_GERMLINE_RESOURCE(germline_resource.flatten().map { vcf -> [[id: vcf.baseName], vcf] })
        germline_resource_tbi = TABIX_GERMLINE_RESOURCE.out.tbi.map { _meta, tbi -> [tbi] }.collect()
        versions = versions.mix(TABIX_GERMLINE_RESOURCE.out.versions)
    }

    known_indels = known_indels_in ? Channel.fromPath(known_indels_in).collect() : Channel.value([])
    known_indels_tbi = known_indels_tbi_in ? Channel.fromPath(known_indels_tbi_in).collect() : Channel.value([])

    if (!known_indels_tbi_in && known_indels_in && (step == 'mapping' || step == "markduplicates" || step == 'prepare_recalibration' || (tools.split(',').contains('haplotypecaller') || tools.split(',').contains('sentieon_haplotyper') || tools.split(',').contains('sentieon_dnascope')))) {
        TABIX_KNOWN_INDELS(known_indels.flatten().map { vcf -> [[id: vcf.baseName], vcf] })
        known_indels_tbi = TABIX_KNOWN_INDELS.out.tbi.map { _meta, tbi -> [tbi] }.collect()
        versions = versions.mix(TABIX_KNOWN_INDELS.out.versions)
    }

    known_snps = known_snps_in ? Channel.fromPath(known_snps_in).collect() : Channel.value([])
    known_snps_tbi = known_snps_tbi_in ? Channel.fromPath(known_snps_tbi_in).collect() : Channel.value([])

    if (!known_snps_tbi_in && known_snps_in && (step == 'mapping' || step == "markduplicates" || step == 'prepare_recalibration' || (tools.split(',').contains('haplotypecaller') || tools.split(',').contains('sentieon_haplotyper')))) {
        TABIX_KNOWN_SNPS(known_snps.flatten().map { vcf -> [[id: vcf.baseName], vcf] })
        known_snps_tbi = TABIX_KNOWN_SNPS.out.tbi.map { _meta, tbi -> [tbi] }.collect()
        versions = versions.mix(TABIX_KNOWN_SNPS.out.versions)
    }

    pon = pon_in ? Channel.fromPath(pon_in).collect() : Channel.value([])
    pon_tbi = pon_tbi_in ? Channel.fromPath(pon_tbi_in).collect() : Channel.value([])

    if (!pon_tbi_in && pon_in && tools.split(',').contains('mutect2')) {
        TABIX_PON(pon.flatten().map { vcf -> [[id: vcf.baseName], vcf] })
        pon_tbi = TABIX_PON.out.tbi.map { _meta, tbi -> [tbi] }.collect()
        versions = versions.mix(TABIX_PON.out.versions)
    }

    // known_sites is made by grouping both the dbsnp and the known snps/indels resources
    // Which can either or both be optional
    known_sites_indels = dbsnp.concat(known_indels).collect()
    known_sites_indels_tbi = dbsnp_tbi.concat(known_indels_tbi).collect()
    known_sites_snps = dbsnp.concat(known_snps).collect()
    known_sites_snps_tbi = dbsnp_tbi.concat(known_snps_tbi).collect()

    // MSI
    if (msisensor2_models_in && msisensor2_models_in.endsWith(".tar.gz") && tools.split(',').contains('msisensor2')) {
        UNTAR_MSISENSOR2_MODELS(Channel.fromPath(file(msisensor2_models_in)).map { archive -> [[id: archive.baseName], archive] })
        msisensor2_models = UNTAR_MSISENSOR2_MODELS.out.untar.map { _meta, extracted_archive -> extracted_archive }.collect()
        versions = versions.mix(UNTAR_MSISENSOR2_MODELS.out.versions)
    }
    else if (msisensor2_models_in && tools.split(',').contains('msisensor2')) {
        msisensor2_models = Channel.fromPath(msisensor2_models_in).collect()
    }
    else {
        msisensor2_models = Channel.value([])
    }

    if (msisensor2_scan_in) {
        msisensor2_scan = Channel.fromPath(msisensor2_scan_in)
    }
    else if (tools.split(',').contains('msisensor2')) {
        MSISENSOR2_SCAN(fasta)
        msisensor2_scan = MSISENSOR2_SCAN.out.scan.map { _meta, list -> [list] }.collect()
        versions = versions.mix(MSISENSOR2_SCAN.out.versions)
    }
    else {
        msisensor2_scan = Channel.value([])
    }

    if (msisensorpro_scan_in) {
        msisensorpro_scan = Channel.fromPath(msisensorpro_scan_in)
    }
    else if (tools.split(',').contains('msisensorpro')) {
        MSISENSORPRO_SCAN(fasta)
        msisensorpro_scan = MSISENSORPRO_SCAN.out.list.map { _meta, list -> [list] }.collect()
        versions = versions.mix(MSISENSORPRO_SCAN.out.versions)
    }
    else {
        msisensorpro_scan = Channel.value([])
    }

    // prepare ascat and controlfreec reference files
    if (!ascat_alleles_in) {
        ascat_alleles = Channel.empty()
    }
    else if (ascat_alleles_in.endsWith(".zip") && tools.split(',').contains('ascat')) {
        UNZIP_ALLELES(Channel.fromPath(file(ascat_alleles_in)).map { archive -> [[id: archive.baseName], archive] })
        ascat_alleles = UNZIP_ALLELES.out.unzipped_archive.map { _meta, extracted_archive -> extracted_archive }.collect()
        versions = versions.mix(UNZIP_ALLELES.out.versions)
    }
    else {
        ascat_alleles = Channel.fromPath(ascat_alleles_in).collect()
    }

    if (!ascat_loci_in) {
        ascat_loci = Channel.empty()
    }
    else if (ascat_loci_in.endsWith(".zip") && tools.split(',').contains('ascat')) {
        UNZIP_LOCI(Channel.fromPath(file(ascat_loci_in)).map { archive -> [[id: archive.baseName], archive] })
        ascat_loci = UNZIP_LOCI.out.unzipped_archive.map { _meta, extracted_archive -> extracted_archive }.collect()
        versions = versions.mix(UNZIP_LOCI.out.versions)
    }
    else {
        ascat_loci = Channel.fromPath(ascat_loci_in).collect()
    }

    if (!ascat_loci_gc_in) {
        ascat_loci_gc = Channel.value([])
    }
    else if (ascat_loci_gc_in.endsWith(".zip") && tools.split(',').contains('ascat')) {
        UNZIP_GC(Channel.fromPath(file(ascat_loci_gc_in)).map { archive -> [[id: archive.baseName], archive] })
        ascat_loci_gc = UNZIP_GC.out.unzipped_archive.map { _meta, extracted_archive -> extracted_archive }.collect()
        versions = versions.mix(UNZIP_GC.out.versions)
    }
    else {
        ascat_loci_gc = Channel.fromPath(ascat_loci_gc_in).collect()
    }

    if (!ascat_loci_rt_in) {
        ascat_loci_rt = Channel.value([])
    }
    else if (ascat_loci_rt_in.endsWith(".zip") && tools.split(',').contains('ascat')) {
        UNZIP_RT(Channel.fromPath(file(ascat_loci_rt_in)).map { archive -> [[id: archive.baseName], archive] })
        ascat_loci_rt = UNZIP_RT.out.unzipped_archive.map { _meta, extracted_archive -> extracted_archive }.collect()
        versions = versions.mix(UNZIP_RT.out.versions)
    }
    else {
        ascat_loci_rt = Channel.fromPath(ascat_loci_rt_in).collect()
    }

    if (!chr_dir_in) {
        chr_dir = Channel.value([])
    }
    else if (chr_dir_in.endsWith(".tar.gz") && tools.split(',').contains('controlfreec')) {
        UNTAR_CHR_DIR(Channel.fromPath(file(chr_dir_in)).map { archive -> [[id: archive.baseName], archive] })
        chr_dir = UNTAR_CHR_DIR.out.untar.map { _meta, extracted_archive -> extracted_archive }.collect()
        versions = versions.mix(UNTAR_CHR_DIR.out.versions)
    }
    else {
        chr_dir = Channel.fromPath(chr_dir_in).collect()
    }

    emit:
    ascat_alleles            // Channel: [ascat_alleles]
    ascat_loci               // Channel: [ascat_loci]
    ascat_loci_gc            // Channel: [ascat_loci_gc]
    ascat_loci_rt            // Channel: [ascat_loci_rt]
    bcftools_annotations     // Channel: [bcftools_annotations]
    bcftools_annotations_tbi // Channel: [bcftools_annotations_tbi]
    chr_dir                  // Channel: [chr_dir/]
    dbsnp                    // Channel: [dbsnp]
    dbsnp_tbi                // Channel: [dbsnp_tbi]
    dict                     // Channel: [meta, dict]
    fasta                    // Channel: [meta, fasta]
    fasta_fai                // Channel: [meta, fasta_fai]
    germline_resource        // Channel: [germline_resource]
    germline_resource_tbi    // Channel: [germline_resource_tbi]
    index_alignment          // Channel: [meta, index_alignment/] either bwa/, bwamem2/ or dragmap/
    known_indels             // Channel: [known_indels]
    known_indels_tbi         // Channel: [known_indels_tbi]
    known_sites_indels       // Channel: [known_sites_indels]
    known_sites_indels_tbi   // Channel: [known_sites_indels_tbi]
    known_sites_snps         // Channel: [known_sites_snps]
    known_sites_snps_tbi     // Channel: [known_sites_snps_tbi]
    known_snps               // Channel: [known_snps]
    known_snps_tbi           // Channel: [known_snps_tbi]
    msisensor2_models        // Channel: [models/]
    msisensor2_scan          // Channel: [genome_msi.list]
    msisensorpro_scan        // Channel: [genome_msi.list]
    pon                      // Channel: [pon]
    pon_tbi                  // Channel: [pon_tbi]
    vep_fasta                // Channel: [meta, vep_fasta]
    versions                 // channel: [versions.yml]
}
