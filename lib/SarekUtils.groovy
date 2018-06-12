class MyUtils {
  // Check if params is in this given list
  static def checkParams(it) {
    return it in [
      'annotate-tools',
      'annotate-VCF',
      'annotateTools',
      'annotateVCF',
      'build',
      'call-name',
      'callName',
      'contact-mail',
      'contactMail',
      'container-path',
      'containerPath',
      'containers',
      'docker',
      'download',
      'explicit-bqsr-needed',
      'explicitBqsrNeeded',
      'genome_base',
      'genome',
      'genomes',
      'help',
      'max_cpus',
      'max_memory',
      'max_time',
      'more',
      'nf-required-version',
      'nfRequiredVersion',
      'no-BAMQC',
      'no-GVCF',
      'no-reports',
      'noBAMQC',
      'noGVCF',
      'noReports',
      'nucleotides-per-second',
      'nucleotidesPerSecond',
      'only-QC',
      'onlyQC',
      'out-dir',
      'outDir',
      'params',
      'project',
      'push',
      'ref-dir',
      'refDir',
      'repository',
      'run-time',
      'runTime',
      'sample-dir',
      'sample',
      'sampleDir',
      'sequencing_center',
      'single-CPUMem',
      'singleCPUMem',
      'singularity',
      'step',
      'strelka-BP',
      'strelkaBP',
      'tag',
      'test',
      'tools',
      'total-memory',
      'totalMemory',
      'vcflist',
      'verbose',
      'version']
  }

  // Loop through all parameters to check their existence and spelling
  static def checkParameterList(list, realList) {
    return list.every{ checkParameterExistence(it, realList) }
  }

  // Check parameter existence
  static def checkParameterExistence(it, list) {
    if (!list.contains(it)) {
      println("Unknown parameter: ${it}")
      return false
    }
    return true
  }

  // Define map of directories
  static def defineDirectoryMap(outDir) {
    return [
      'nonRealigned'     : "${outDir}/Preprocessing/NonRealigned",
      'nonRecalibrated'  : "${outDir}/Preprocessing/NonRecalibrated",
      'recalibrated'     : "${outDir}/Preprocessing/Recalibrated",
      'ascat'            : "${outDir}/VariantCalling/Ascat",
      'freebayes'        : "${outDir}/VariantCalling/FreeBayes",
      'gvcf-hc'          : "${outDir}/VariantCalling/HaplotypeCallerGVCF",
      'haplotypecaller'  : "${outDir}/VariantCalling/HaplotypeCaller",
      'manta'            : "${outDir}/VariantCalling/Manta",
      'mutect1'          : "${outDir}/VariantCalling/MuTect1",
      'mutect2'          : "${outDir}/VariantCalling/MuTect2",
      'strelka'          : "${outDir}/VariantCalling/Strelka",
      'strelkabp'        : "${outDir}/VariantCalling/StrelkaBP",
      'snpeff'           : "${outDir}/Annotation/SnpEff",
      'vep'              : "${outDir}/Annotation/VEP",
      'bamQC'            : "${outDir}/Reports/bamQC",
      'bcftoolsStats'    : "${outDir}/Reports/BCFToolsStats",
      'fastQC'           : "${outDir}/Reports/FastQC",
      'markDuplicatesQC' : "${outDir}/Reports/MarkDuplicates",
      'multiQC'          : "${outDir}/Reports/MultiQC",
      'samtoolsStats'    : "${outDir}/Reports/SamToolsStats",
      'snpeffReports'    : "${outDir}/Reports/SnpEff",
      'vcftools'         : "${outDir}/Reports/VCFTools",
      'version'          : "${outDir}/Reports/ToolsVersion"
    ]
  }

  // Compare params to list of verified params
  static def isAllowedParams(params) {
    final test = true
    params.each{
      if (!checkParams(it.toString().split('=')[0])) {
        println "params ${it.toString().split('=')[0]} is unknown"
        test = false
      }
    }
    return test
  }
}
