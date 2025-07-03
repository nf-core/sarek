// Helper functions for pipeline tests

class UTILS {

    public static def get_assertion = { Map args ->
        def outdir = args.outdir
        def stub = args.stub
        def include_txt = args.include_txt
        def vcf_gzip_lines = args.vcf_gzip_lines

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
        // txt_files: MuSE txt files
        def txt_files = getAllFilesFromDir(outdir, include: ['**/*.MuSE.txt'])
        // vcf_files: All vcf files
        def vcf_files = getAllFilesFromDir(outdir, include: ['**/*.vcf{,.gz}'], ignore: ['**/test{N,T}.germline.vcf{,.gz}'])

        def assertion = [
            removeFromYamlMap("${outdir}/pipeline_info/nf_core_sarek_software_mqc_versions.yml", "Workflow"),
            stable_name
        ]

        if (!stub) {
            assertion.add(stable_content.isEmpty() ? 'No stable content' : stable_content)
            assertion.add(bam_files.isEmpty() ? 'No BAM files' : bam_files.collect { file -> file.getName() + ":md5," + bam(file.toString()).readsMD5 })
            assertion.add(cram_files.isEmpty() ? 'No CRAM files' : cram_files.collect { file -> file.getName() + ":md5," + cram(file.toString(), fasta).readsMD5 })
            if (include_txt) {
                assertion.add(txt_files.isEmpty() ? 'No TXT files' : txt_files.collect{ file -> file.getName() + ":md5," + file.readLines()[2..-1].join('\n').md5() })
            }
            if (vcf_gzip_lines) {
                assertion.add(vcf_files.isEmpty() ? 'No VCF files' : vcf_files.collect { file -> [file.getName(), path(file.toString()).linesGzip[vcf_gzip_lines], path(file.toString()).vcf.summary] })
            } else {
                assertion.add(vcf_files.isEmpty() ? 'No VCF files' : vcf_files.collect { file -> file.getName() + ":md5," + path(file.toString()).vcf.variantsMD5 })
            }
        }

        return assertion
    }

    public static def get_test = { scenario ->
        return {
            if (scenario.stub && scenario.gpu) {
                options "-stub"
                tag "gpu_stub"
            } else if (scenario.stub && !scenario.gpu) {
                options "-stub"
                tag "cpu_stub"
            } else if (!scenario.stub && scenario.gpu) {
                tag "gpu"
            } else if (!scenario.stub && !scenario.gpu) {
                tag "cpu"
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
                        *UTILS.get_assertion(outdir: params.outdir, stub: scenario.stub, include_txt: scenario.include_txt, vcf_gzip_lines: scenario.vcf_gzip_lines)
                    ).match() }
                )
            }
        }
    }
}
