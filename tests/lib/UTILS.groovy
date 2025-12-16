// Helper functions for pipeline tests

class UTILS {

    public static def get_assertion = { Map args ->
        // Mandatory, as we always need an outdir
        def outdir = args.outdir

        // Get scenario and extract all properties dynamically
        def scenario = args.scenario ?: [:]

        // Pass down workflow for std capture
        def workflow = args.workflow

        // stable_name: All files + folders in ${outdir}/ with a stable name
        def stable_name = getAllFilesFromDir(outdir, relative: true, includeDir: true, ignore: ['pipeline_info/*.{html,json,txt}'])
        // stable_content: All files in ${outdir}/ with stable content
        def stable_content = getAllFilesFromDir(outdir, ignoreFile: 'tests/.nftignore', ignore: [scenario.ignoreFiles ])
        // bam_files: All bam files
        def bam_files = getAllFilesFromDir(outdir, include: ['**/*.bam'], ignore: [scenario.ignoreFiles ])
        // cram_files: All cram files
        def cram_files = getAllFilesFromDir(outdir, include: ['**/*.cram'], ignore: [scenario.ignoreFiles ])
        // Fasta file for cram verification with nft-bam
        def fasta_base = 'https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/'
        def fasta = fasta_base + 'genomics/homo_sapiens/genome/genome.fasta'
        // txt_files: MuSE txt files
        def txt_files = getAllFilesFromDir(outdir, include: ['**/*.MuSE.txt'])
        // vcf_files: All vcf files
        def vcf_files = getAllFilesFromDir(outdir, include: ['**/*.vcf{,.gz}'], ignore: [scenario.ignoreFiles ])
        // freebayes_unfiltered: vcf files from freebayes without quality filtering
        def freebayes_unfiltered = getAllFilesFromDir(outdir, include: ['**/*.freebayes.vcf.gz'])
        // varlociraptor vcf
        def varlociraptor_vcf = getAllFilesFromDir(outdir, include: ['**/*.varlociraptor.{vcf}{,.gz}'])

        def assertion = []

        if (!scenario.failure) {
            assertion.add(workflow.trace.succeeded().size())
            assertion.add(removeFromYamlMap("${outdir}/pipeline_info/nf_core_sarek_software_mqc_versions.yml", "Workflow"))
        }

        // At least always pipeline_info/ is created and stable
        assertion.add(stable_name)

        if (!scenario.stub) {
            assertion.add(stable_content.isEmpty() ? 'No stable content' : stable_content)
            assertion.add(bam_files.isEmpty() ? 'No BAM files' : bam_files.collect { file -> file.getName() + ":md5," + bam(file.toString()).readsMD5 })
            assertion.add(cram_files.isEmpty() ? 'No CRAM files' : cram_files.collect { file -> file.getName() + ":md5," + cram(file.toString(), fasta).readsMD5 })
            if (scenario.include_muse_txt) {
                // It will skip the first line of the txt file
                assertion.add(txt_files.isEmpty() ? 'No TXT files' : txt_files.collect{ file -> file.getName() + ":md5," + file.readLines()[2..-1].join('\n').md5() })
            }
            if (scenario.include_freebayes_unfiltered) {
                // It will only print the vcf summary to avoid differing md5sums because of small differences in QUAL score
                assertion.add(freebayes_unfiltered.isEmpty() ? 'No Freebayes unfiltered VCF files' : freebayes_unfiltered.collect { file -> [ file.getName(), path(file.toString()).vcf.summary ] })
            }
            if (scenario.no_vcf_md5sum) {
                // Will print the summary instead of the md5sum for vcf files
                assertion.add(vcf_files.isEmpty() ? 'No VCF files' : vcf_files.collect { file -> [ file.getName(), path(file.toString()).vcf.summary ] })
            } else {
                assertion.add(vcf_files.isEmpty() ? 'No VCF files' : vcf_files.collect { file -> file.getName() + ":md5," + path(file.toString()).vcf.variantsMD5 })
                if (scenario.include_varlociraptor_vcf) {
                    // It will use the summary method to extract the vcf file content
                    assertion.add(varlociraptor_vcf.isEmpty() ? 'No Varlociraptor VCF files' : varlociraptor_vcf.collect { file -> file.getName() + ":summary," + path(file.toString()).vcf.summary })
                }
            }
        }

        // Always capture stdout and stderr for any WARN message
        assertion.add(filterNextflowOutput(workflow.stderr + workflow.stdout, include: ["WARN"] ) ?: "No warnings")

        // Capture std for snapshot
        // Allow to capture either stderr, stdout or both
        // Additional possibilities to include and/or ignore some string
        if (scenario.snapshot) {
            def workflow_std = []

            scenario.snapshot.split(',').each { std ->
                if (std in ['stderr', 'stdout']) { workflow_std.add(workflow."$std") }
            }

            if (scenario.snapshot_include) {
                assertion.add(filterNextflowOutput(workflow_std.flatten(), ignore: [scenario.snapshot_ignore], include:[scenario.snapshot_include]))
            } else {
                assertion.add(filterNextflowOutput(workflow_std.flatten(), ignore: [scenario.snapshot_ignore]))
            }
        }

        return assertion
    }

    public static def get_test = { scenario ->
        // This function returns a closure that will be used to run the test and the assertion
        // It will create tags or options based on the scenario

        return {
            // If the test is for a gpu, we add the gpu tag
            // Otherwise, we add the cpu tag
            // If the tests has no conda incompatibilities
            // then we append "_conda" to the cpu/gpu tag
            // If the test is for a stub, we add options -stub
            // And we append "_stub" to the cpu/gpu tag

            // All options should be:
            // gpu (this is the default for gpu)
            // cpu (this is the default for tests without conda)
            // gpu_conda (this should never happen)
            // cpu_conda (this is the default for tests with conda compatibility)
            // gpu_stub
            // cpu_stub
            // gpu_conda_stub (this should never happen)
            // cpu_conda_stub

            if (scenario.stub) {
                options "-stub"
            }

            if (scenario.gpu) {
                tag "gpu${!scenario.no_conda ? '_conda' : ''}${scenario.stub ? '_stub' : ''}"
            }

            if (!scenario.gpu) {
                tag "cpu${!scenario.no_conda ? '_conda' : ''}${scenario.stub ? '_stub' : ''}"
            }

            // If a tag is provided, add it to the test
            if (scenario.tag) {
                tag scenario.tag
            }

            when {
                params {
                    // Mandatory, as we always need an outdir
                    outdir = "${outputDir}"
                    // Apply scenario-specific params
                    scenario.params.each { key, value ->
                        delegate."$key" = value
                    }
                }
            }

            then {
                // Assert failure/success, and fails early so we don't pollute console with massive diffs
                if (scenario.failure) {
                    assert workflow.failed
                } else {
                    assert workflow.success
                }
                assertAll(
                    { assert snapshot(
                        // All assertions based on the scenario
                        *UTILS.get_assertion(
                            outdir: params.outdir,
                            scenario: scenario,
                            workflow: workflow
                        )
                    ).match() }
                )
            }
            cleanup {
                if (System.getenv('NFT_CLEANUP')) {
                    println ""
                    println "CLEANUP"
                    println "Set NFT_CLEANUP to false to disable."
                    println "The following folders will be deleted:"
                    println "- ${workDir}"

                    new File("${workDir}").deleteDir()
                }
            }
        }
    }
}
