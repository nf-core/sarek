class MyUtils {
  static def checkParams(it) {
    // Check if params is in this given list
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

  static def checkParameterList(list, realList) {
    // Loop through all parameters to check their existence and spelling
    return list.every{ checkParameterExistence(it, realList) }
  }

  static def checkParameterExistence(it, list) {
    // Check parameter existence
    if (!list.contains(it)) {
      println("Unknown parameter: ${it}")
      return false
    }
    return true
  }

  static def isAllowedParams(params) {
    // Compare params to list of verified params
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
