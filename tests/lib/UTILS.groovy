// Helper functions for pipeline tests

class UTILS {

    public static def get_assertion = { Map args ->
        // Mandatory, as we always need an outdir
        def outdir = args.outdir

        // Use this args to run the test with stub
        // It will disable all assertions but versions and stable_name
        def stub = args.stub

        // Use this args to include muse txt files in the assertion
        // It will skip the first line of the txt file
        def include_muse_txt = args.include_muse_txt

        // Use this args to include freebayes unfiltered vcf files in the assertion
        // It will only print the vcf summary to avoid differing md5sums because of small differences in QUAL score
        def include_freebayes_unfiltered = args.include_freebayes_unfiltered

        // Will print the summary instead of the md5sum for vcf files
        def no_vcf_md5sum = args.no_vcf_md5sum

        // Use this args to include varlociraptor vcf files in the assertion
        // It will use the summary method to extract the vcf file content
        def include_varlociraptor_vcf = args.include_varlociraptor_vcf

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
        def vcf_files = getAllFilesFromDir(outdir, include: ['**/*.vcf{,.gz}'], ignore: ['**/test{N,T}.germline.vcf{,.gz}', '**/*.freebayes.vcf{,.gz}',  '**/*.varlociraptor.{vcf}{,.gz}'])
        // freebayes_unfiltered: vcf files from freebayes without quality filtering
        def freebayes_unfiltered = getAllFilesFromDir(outdir, include: ['**/*.freebayes.vcf.gz'])
        // varlociraptor vcf
        def varlociraptor_vcf = getAllFilesFromDir(outdir, include: ['**/*.varlociraptor.{vcf}{,.gz}'] )


        def assertion = []

        assertion.add(removeFromYamlMap("${outdir}/pipeline_info/nf_core_sarek_software_mqc_versions.yml", "Workflow"))
        assertion.add(stable_name)

        if (!stub) {
            assertion.add(stable_content.isEmpty() ? 'No stable content' : stable_content)
            assertion.add(bam_files.isEmpty() ? 'No BAM files' : bam_files.collect { file -> file.getName() + ":md5," + bam(file.toString()).readsMD5 })
            assertion.add(cram_files.isEmpty() ? 'No CRAM files' : cram_files.collect { file -> file.getName() + ":md5," + cram(file.toString(), fasta).readsMD5 })
            if (include_muse_txt) {
                assertion.add(txt_files.isEmpty() ? 'No TXT files' : txt_files.collect{ file -> file.getName() + ":md5," + file.readLines()[2..-1].join('\n').md5() })
            }
            if (include_freebayes_unfiltered) {
                assertion.add(freebayes_unfiltered.isEmpty() ? 'No Freebayes unfiltered VCF files' : freebayes_unfiltered.collect { file -> [ file.getName(), path(file.toString()).vcf.summary ] })
            }
            if (no_vcf_md5sum) {
                assertion.add(vcf_files.isEmpty() ? 'No VCF files' : vcf_files.collect { file -> [ file.getName(), path(file.toString()).vcf.summary ] })
            } else {
                assertion.add(vcf_files.isEmpty() ? 'No VCF files' : vcf_files.collect { file -> file.getName() + ":md5," + path(file.toString()).vcf.variantsMD5 })
                if (include_varlociraptor_vcf) {
                    assertion.add(varlociraptor_vcf.isEmpty() ? 'No Varlociraptor VCF files' : varlociraptor_vcf.collect { file -> file.getName() + ":summary," + path(file.toString()).vcf.summary })
                }
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
            // If the test is for a stub, we add options -stub
            // And we append "_stub" to the cpu/gpu tag

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
                // Assert failure
                if (scenario.failure) {
                    // Early failure, so we don't pollute console with massive diffs
                    assert workflow.failed
                    // Check stdout if specified
                    if (scenario.stdout) {
                        assertAll(
                            { assert workflow.stdout.toString().contains(scenario.stdout) }
                        )
                    }
                    // Check stderr if specified
                    if (scenario.stderr) {
                        { assert snapshot(
                            workflow.stderr.toString().replaceAll(/\x1B\[[0-9;]*m/, '').replaceAll(/^\[/, '').replaceAll(/\]$/, '').replaceAll(/, /, ',').split(",").findAll { !it.matches(/.*Nextflow [0-9]+\.[0-9]+\.[0-9]+ is available.*/) }[scenario.stderr]
                        ).match() }
                    }
                // Assert success
                } else {
                    // Early failure, so we don't pollute console with massive diffs
                    assert workflow.success
                    assertAll(
                        { assert snapshot(
                            // Number of successful tasks
                            workflow.trace.succeeded().size(),
                            // All assertions based on the scenario
                            *UTILS.get_assertion(include_freebayes_unfiltered: scenario.include_freebayes_unfiltered, include_muse_txt: scenario.include_muse_txt, include_varlociraptor_vcf: scenario.include_varlociraptor_vcf, no_vcf_md5sum: scenario.no_vcf_md5sum, outdir: params.outdir, stub: scenario.stub)
                        ).match() }
                    )
                    // Check stdout if specified
                    if (scenario.stdout) {
                        assert workflow.stdout.toString().contains(scenario.stdout)
                    }
                }
            }
        }
    }
}
