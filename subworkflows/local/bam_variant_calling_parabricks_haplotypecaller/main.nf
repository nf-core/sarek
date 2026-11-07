//
// PARABRICKS HAPLOTYPECALLER germline variant calling (GPU-accelerated)
//

include { PARABRICKS_DBSNP                            } from '../../../modules/nf-core/parabricks/dbsnp'
include { PARABRICKS_HAPLOTYPECALLER                  } from '../../../modules/nf-core/parabricks/haplotypecaller'
include { HTSLIB_BGZIPTABIX as TABIX_BGZIPTABIX_VCF   } from '../../../modules/nf-core/htslib/bgziptabix'
include { HTSLIB_BGZIPTABIX as TABIX_BGZIPTABIX_GVCF  } from '../../../modules/nf-core/htslib/bgziptabix'

workflow BAM_VARIANT_CALLING_PARABRICKS_HAPLOTYPECALLER {
    take:
    cram                    // channel: [mandatory] [ meta, cram, crai ]
    fasta                   // channel: [mandatory] [ meta, fasta ]
    dbsnp                   // channel: [optional]  [] or [ dbsnp.vcf.gz ]
    dbsnp_tbi               // channel: [optional]  [] or [ dbsnp.vcf.gz.tbi ]
    intervals               // channel: [mandatory] [ interval.bed, num_intervals ], used to scatter joint genotyping
    intervals_bed_combined  // channel: [optional]  [] or [ intervals.bed ]

    main:
    // Combine each sample with the (optional) intervals list
    // intervals_bed_combined emits [] (no intervals) or [file] (one combined BED)
    // When no_intervals, the empty list contributes 0 elements to the combined tuple,
    // so check the tuple size to safely extract the optional 4th element.
    cram_intervals = cram
        .combine(intervals_bed_combined)
        .map { cram_combined ->
            def (meta, cram_, crai) = cram_combined
            def intervals_ = cram_combined.size() > 3 ? cram_combined[3] : []
            [ meta + [ variantcaller:'parabricks_haplotypecaller' ], cram_, crai, intervals_ ]
        }

    PARABRICKS_HAPLOTYPECALLER(cram_intervals, fasta)

    // pbrun haplotypecaller does not accept --dbsnp, so rsIDs are added in a separate step.
    // dbsnp/dbsnp_tbi are wrapped so an empty list stays a single element and can be branched on.
    vcf_dbsnp = PARABRICKS_HAPLOTYPECALLER.out.vcf
        .combine(dbsnp.map { dbsnp_ -> [ dbsnp_ ] })
        .combine(dbsnp_tbi.map { dbsnp_tbi_ -> [ dbsnp_tbi_ ] })
        .branch { _meta, _vcf, dbsnp_, _dbsnp_tbi ->
            annotate: dbsnp_
            skip: true
        }

    PARABRICKS_DBSNP(vcf_dbsnp.annotate)

    vcf_raw = PARABRICKS_DBSNP.out.vcf
        .mix(vcf_dbsnp.skip.map { meta, vcf_, _dbsnp, _dbsnp_tbi -> [ meta, vcf_ ] })

    // Depending on ext.args (e.g. '--gvcf'), PARABRICKS_HAPLOTYPECALLER emits either an
    // uncompressed VCF or an already-bgzipped gVCF, never both. HTSLIB_BGZIPTABIX detects the
    // input compression itself, so it both compresses the plain VCF and indexes
    // the already-compressed gVCF.
    TABIX_BGZIPTABIX_VCF(vcf_raw.map { meta, vcf_ -> [ meta, vcf_, [], [] ] }, 'compress', true, 'vcf')
    TABIX_BGZIPTABIX_GVCF(PARABRICKS_HAPLOTYPECALLER.out.gvcf.map { meta, gvcf_ -> [ meta, gvcf_, [], [] ] }, 'compress', true, 'g.vcf')

    vcf_tbi_joined  = TABIX_BGZIPTABIX_VCF.out.output.join(TABIX_BGZIPTABIX_VCF.out.index)
    gvcf_tbi_joined = TABIX_BGZIPTABIX_GVCF.out.output.join(TABIX_BGZIPTABIX_GVCF.out.index)

    vcf      = vcf_tbi_joined.map  { meta, vcf_,  _tbi -> [ meta, vcf_  ] }
    tbi      = vcf_tbi_joined.map  { meta, _vcf,  tbi_ -> [ meta, tbi_  ] }
    gvcf     = gvcf_tbi_joined.map { meta, gvcf_, _tbi -> [ meta, gvcf_ ] }
    gvcf_tbi = gvcf_tbi_joined.map { meta, _gvcf, gvcf_tbi_ -> [ meta, gvcf_tbi_ ] }

    // The gVCF covers all intervals in one file; GenomicsDBImport restricts each import to its own interval
    gvcf_tbi_intervals = gvcf_tbi_joined
        .combine(intervals)
        .map { meta, gvcf_, gvcf_tbi_, intervals_, num_intervals -> [ meta + [ num_intervals:num_intervals ], gvcf_, gvcf_tbi_, intervals_ ] }

    emit:
    vcf                 // channel: [ val(meta), vcf.gz ]
    tbi                 // channel: [ val(meta), vcf.gz.tbi ]
    gvcf                // channel: [ val(meta), g.vcf.gz ]
    gvcf_tbi            // channel: [ val(meta), g.vcf.gz.tbi ]
    gvcf_tbi_intervals  // channel: [ val(meta), g.vcf.gz, g.vcf.gz.tbi, interval.bed ], for joint genotyping
}
