import static nextflow.Nextflow.file
import nextflow.Channel

class SarekUtils {

  // Check if a row has the expected number of item
  static def checkNumberOfItem(row, number) {
    if (row.size() != number) exit 1, "Malformed row in TSV file: ${row}, see --help for more information"
    return true
  }

  // Check parameter existence
  static def checkParameterExistence(it, list) {
    if (!list.contains(it)) {
      println("Unknown parameter: ${it}")
      return false
    }
    return true
  }

  // Compare each parameter with a list of parameters
  static def checkParameterList(list, realList) {
    return list.every{ checkParameterExistence(it, realList) }
  }

  // Return element in list of allowed params
  static def checkParams(it) {
    return it in [
      'ac-loci-GC',
      'ac-loci',
      'acLoci',
      'acLociGC',
      'annotate-tools',
      'annotate-VCF',
      'annotateTools',
      'annotateVCF',
      'annotation_cache',
      'awsqueue_tiny',
      'awsqueue',
      'build',
      'bwa-index',
      'bwaIndex',
      'call-name',
      'callName',
      'contact-mail',
      'contactMail',
      'container-path',
      'containerPath',
      'containers',
      'dbsnp-index',
      'dbsnp',
      'dbsnpIndex',
      'docker',
      'download',
      'explicit-bqsr-needed',
      'explicitBqsrNeeded',
      'genome_base',
      'genome-dict',
      'genome-file',
      'genome-index',
      'genome',
      'genomeDict',
      'genomeFile',
      'genomeIndex',
      'genomes',
      'help',
      'intervals',
      'known-indels-index',
      'known-indels',
      'knownIndels',
      'knownIndelsIndex',
      'local-report-dir',
      'localReportDir',
      'markdup_java_options',
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
      'publish-dir-mode',
      'publishDirMode',
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
      'snp-eff_cache',
      'snpEff_cache',
      'snpeff-db',
      'snpeffDb',
      'step',
      'strelka-BP',
      'strelkaBP',
      'tag',
      'target-BED',
      'targetBED',
      'test',
      'tools',
      'total-memory',
      'totalMemory',
      'vcflist',
      'vep_cache',
      'vep-cache-version',
      'vepCacheVersion',
      'verbose',
      'version']
  }

  // Loop through all the references files to check their existence
  static def checkReferenceMap(referenceMap) {
    referenceMap.every {
      referenceFile, fileToCheck ->
      SarekUtils.checkRefExistence(referenceFile, fileToCheck)
    }
  }

  // Loop through all the references files to check their existence

  static def checkRefExistence(referenceFile, fileToCheck) {
    if (fileToCheck instanceof List) return fileToCheck.every{ SarekUtils.checkRefExistence(referenceFile, it) }
    def f = file(fileToCheck)
    // this is an expanded wildcard: we can assume all files exist
    if (f instanceof List && f.size() > 0) return true
    else if (!f.exists()) {
			println  "Missing references: ${referenceFile} ${fileToCheck}"
      return false
    }
    return true
  }

  // Channeling the TSV file containing BAM.
  // Format is: "subject gender status sample bam bai"
  static def extractBams(tsvFile, mode) {
    Channel.from(tsvFile)
      .splitCsv(sep: '\t')
      .map { row ->
        SarekUtils.checkNumberOfItem(row, 6)
        def idPatient = row[0]
        def gender    = row[1]
        def status    = SarekUtils.returnStatus(row[2].toInteger())
        def idSample  = row[3]
        def bamFile   = SarekUtils.returnFile(row[4])
        def baiFile   = SarekUtils.returnFile(row[5])

        if (!SarekUtils.hasExtension(bamFile,".bam")) exit 1, "File: ${bamFile} has the wrong extension. See --help for more information"
        if (!SarekUtils.hasExtension(baiFile,".bai")) exit 1, "File: ${baiFile} has the wrong extension. See --help for more information"

        if (mode == "germline") return [ idPatient, status, idSample, bamFile, baiFile ]
        else return [ idPatient, gender, status, idSample, bamFile, baiFile ]
      }
  }

  // Extract gender from Channel as it's only used for CNVs
  static def extractGenders(channel) {
    def genders = [:]
    channel = channel.map{ it ->
      def idPatient = it[0]
      def gender = it[1]
      genders[idPatient] = gender
      [idPatient] + it[2..-1]
    }
    [genders, channel]
  }

  // Check file extension
  static def hasExtension(it, extension) {
    it.toString().toLowerCase().endsWith(extension.toLowerCase())
  }

  // Compare params to list of verified params
  static def isAllowedParams(params) {
    def test = true
    params.each{
      if (!checkParams(it.toString().split('=')[0])) {
        println "params ${it.toString().split('=')[0]} is unknown"
        test = false
      }
    }
    return test
  }

  // Return file if it exists
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
    println "    ____        _____               _    "
    println "  .' _  `.     / ____|             | |   "
    println " /  |\\`-_ \\   | (___  ___  _ __ __ | | __"
    println "|   | \\  `-|   \\___ \\/__ \\| ´__/ _\\| |/ /"
    println " \\ |   \\  /    ____) | __ | | |  __|   < "
    println "  `|____\\'    |_____/\\____|_|  \\__/|_|\\_\\"
  }

}
