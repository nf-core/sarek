import static nextflow.Nextflow.file
import nextflow.Channel

class SarekUtils {
  // Check file extension
  static def checkFileExtension(it, extension) {
    if (!it.toString().toLowerCase().endsWith(extension.toLowerCase())) exit 1, "File: ${it} has the wrong extension: ${extension} see --help for more information"
  }

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

  static def checkParameterExistence(it, list) {
    // Check parameter existence
    if (!list.contains(it)) {
      println("Unknown parameter: ${it}")
      return false
    }
    return true
  }

  static def checkParameterList(list, realList) {
    // Loop through all parameters to check their existence and spelling
    return list.every{ checkParameterExistence(it, realList) }
  }

  // Loop through all the references files to check their existence
  static def checkReferenceMap(referenceMap) {
    referenceMap.every {
      referenceFile, fileToCheck ->
      SarekUtils.checkRefExistence(referenceFile, fileToCheck)
    }
  }

  static def checkRefExistence(referenceFile, fileToCheck) {
    if (fileToCheck instanceof List) return fileToCheck.every{ SarekUtils.checkRefExistence(referenceFile, it) }
    def f = file(fileToCheck)
    // this is an expanded wildcard: we can assume all files exist
    if (f instanceof List && f.size() > 0) return true
    else if (!f.exists()) {
      this.log.info  "Missing references: ${referenceFile} ${fileToCheck}"
      return false
    }
    return true
  }

  // Channeling the TSV file containing BAM.
  // Format is: "subject gender status sample bam bai"
  static def extractBams(tsvFile, mode) {
    Channel
      .from(tsvFile.readLines())
      .map{line ->
        def list      = SarekUtils.returnTSV(line.split(),6)
        def idPatient = list[0]
        def gender    = list[1]
        def status    = SarekUtils.returnStatus(list[2].toInteger())
        def idSample  = list[3]
        def bamFile   = SarekUtils.returnFile(list[4])
        def baiFile   = SarekUtils.returnFile(list[5])

        SarekUtils.checkFileExtension(bamFile,".bam")
        SarekUtils.checkFileExtension(baiFile,".bai")

        if (mode == "germline") return [ idPatient, status, idSample, bamFile, baiFile ]
        else return [ idPatient, gender, status, idSample, bamFile, baiFile ]
      }
  }

  static def extractGenders(channel) {
    def genders = [:]  // an empty map
    channel = channel.map{ it ->
      def idPatient = it[0]
      def gender = it[1]
      genders[idPatient] = gender
      [idPatient] + it[2..-1]
    }
    [genders, channel]
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

  // return file if it exists
  static def returnFile(it) {
    if (!file(it).exists()) exit 1, "Missing file in TSV file: ${it}, see --help for more information"
    return file(it)
  }

  // Return status [0,1]
  // 0 == Normal, 1 == Tumor
  static def returnStatus(it) {
    if (!(it in [0, 1])) exit 1, "Status is not recognized in TSV file: ${it}, see --help for more information"
    return it
  }

  // Sarek ascii art
  static def sarek_ascii() {
    println "    ____        _____                _     "
    println "  .' _  `.     / ____|              | |    "
    println " /  |\\`-_ \\   | (___   __ _ _ __ ___| | __ "
    println "|   | \\  `-|   \\___ \\ / _` | '__/ __| |/ / "
    println " \\ |   \\  /    ____) | (_| | | |  __|   <  "
    println "  `|____\\'    |_____/ \\__,_|_|  \\___|_|\\_\\ "
  }

  // return TSV if it has the correct number of items in row
  static def returnTSV(it, number) {
    if (it.size() != number) exit 1, "Malformed row in TSV file: ${it}, see --help for more information"
    return it
  }

}
