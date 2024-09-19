class UTILS_PIPELINE {
    // These are the files to exclude when we want to snapshot
    static List<String> exclusionRegexesForUnstableFileContents = [
        // To exclude the pipeline_software_mqc_versions.yml file that contains the Nextflow version
        /nf_core_.*_software_mqc_versions\.yml/,

        // To exclude trimgalore
        /.*\.fastq\.gz_trimming_report\.txt/
    ]
}
