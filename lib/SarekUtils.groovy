import static nextflow.Nextflow.file
import nextflow.Channel

class SarekUtils {

  // Check file extension
  static def checkFileExtension(it, extension) {
    if (!it.toString().toLowerCase().endsWith(extension.toLowerCase())) exit 1, "File: ${it} has the wrong extension: ${extension} see --help for more information"
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
      this.log.info  "Missing references: ${referenceFile} ${fileToCheck}"
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

  // Return TSV if it has the correct number of items in row
  static def returnTSV(it, number) {
    if (it.size() != number) exit 1, "Malformed row in TSV file: ${it}, see --help for more information"
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

}
