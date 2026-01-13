//
// Run SnpSift to annotate VCF files with multiple databases
//

include { SNPSIFT_ANNOTATE_MULTI } from '../../../modules/local/snpsift/annotate_multi/main'

workflow VCF_ANNOTATE_SNPSIFT {
    take:
    vcf_tbi          // channel: [ val(meta), path(vcf), path(tbi) ]
    snpsift_db_list  // value: list of database configurations

    main:
    ch_versions = Channel.empty()

    // Stage all database files
    db_files = Channel.fromList(snpsift_db_list.collect { file(it.file, checkIfExists: true) }).collect()
    db_tbis = Channel.fromList(snpsift_db_list.collect { file(it.file_tbi, checkIfExists: true) }).collect()

    // Prepare database list with staged filenames
    db_configs = snpsift_db_list.collect { db ->
        [
            name: db.name,
            file: file(db.file).name,  // Use just the filename for the staged file
            fields: db.fields ?: '',
            id: db.id ?: false,
            prefix: db.prefix ?: ''
        ]
    }

    // Annotate with all databases
    SNPSIFT_ANNOTATE_MULTI(
        vcf_tbi,
        db_configs,
        db_files,
        db_tbis
    )

    ch_versions = ch_versions.mix(SNPSIFT_ANNOTATE_MULTI.out.versions)

    emit:
    vcf_tbi  = SNPSIFT_ANNOTATE_MULTI.out.vcf  // channel: [ val(meta), path(vcf), path(tbi) ]
    versions = ch_versions                      // channel: [ path(versions.yml) ]
}
