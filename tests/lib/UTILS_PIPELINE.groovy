class UTILS_PIPELINE {
    // These are the files to exclude when we want to snapshot
    static List<String> exclusionRegexesForUnstableFileContents = [
        // To exclude the pipeline_software_mqc_versions.yml file that contains the Nextflow version
        /nf_core_.*_software_mqc_versions\.yml/,

        // To exclude csv files created by the pipeline that contains absolute paths
        /(markduplicates|markduplicates_no_table|recalibrated|variantcalled)\.csv/,

        // To exclude from multiqc
        /multiqc\.log/,
        /multiqc_data\.json/,
        /multiqc_sources\.txt/,
        /multiqc_report\.html/,

        // To exclude cram and index files
        /.*\.(md|recal)\.cram/,

        // To exclude markduplicates metrics
        /.*\.md\.cram\.metrics/,

        // To exclude vcf and index files
        /.*\.vcf\.gz/,
        /.*\.vcf\.gz\.tbi/,

        // To exclude fastqc
        /.*_fastq\.(html|zip)/,

        // To exclude trimgalore
        /.*\.fastq\.gz_trimming_report\.txt/
    ]
}
