// Helper functions for pipeline tests

class UTILS {

    public static def get_assertion = { outdir, stub ->
        // stable_name: All files + folders in ${outdir}/ with a stable name
        def stable_name = getAllFilesFromDir(outdir, relative: true, includeDir: true, ignore: ['pipeline_info/*.{html,json,txt}'])
        // stable_content: All files in ${outdir}/ with stable content
        def stable_content = getAllFilesFromDir(outdir, ignoreFile: 'tests/.nftignore')
        // bam_files: All bam files
        def bam_files = getAllFilesFromDir(outdir, include: ['**/*.bam'])
        // cram_files: All cram files
        def cram_files = getAllFilesFromDir(outdir, include: ['**/*.cram'])
        // Fasta file for cram verification with nft-bam
        def fasta_base = 'https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/'
        def fasta = fasta_base + 'genomics/homo_sapiens/genome/genome.fasta'
        // vcf_files: All vcf files
        def vcf_files = getAllFilesFromDir(outdir, include: ['**/*.vcf{,.gz}'])

        if (stub) {
            return [
                removeFromYamlMap("${outdir}/pipeline_info/nf_core_sarek_software_mqc_versions.yml", "Workflow"),
                stable_name
            ]
        } else {
            return [
                removeFromYamlMap("${outdir}/pipeline_info/nf_core_sarek_software_mqc_versions.yml", "Workflow"),
                stable_name,
                stable_content.isEmpty() ? 'No stable content' : stable_content,
                bam_files.isEmpty() ? 'No BAM files' : bam_files.collect { file -> file.getName() + ":md5," + bam(file.toString()).readsMD5 },
                cram_files.isEmpty() ? 'No CRAM files' : cram_files.collect { file -> file.getName() + ":md5," + cram(file.toString(), fasta).readsMD5 },
                vcf_files.isEmpty() ? 'No VCF files' : vcf_files.collect { file -> file.getName() + ":md5," + path(file.toString()).vcf.variantsMD5 }
            ]
        }
    }

    public static def get_test = { scenario ->
        return {
            if (scenario.stub) {
                options "-stub"
                tag "stub"
            } else {
                tag "no_stub"
            }

            if (scenario.tag) {
                tag scenario.tag
            }

            when {
                params {
                    outdir = "${outputDir}"
                    // Apply scenario-specific parameters
                    scenario.params.each { key, value ->
                        delegate."$key" = value
                    }
                }
            }

            then {
                assert workflow.success
                assertAll(
                    { assert snapshot(
                        workflow.trace.succeeded().size(),
                        *UTILS.get_assertion(params.outdir, scenario.stub)
                    ).match() }
                )
            }
        }
    }
}
