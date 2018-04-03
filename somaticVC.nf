#!/usr/bin/env nextflow

/*
vim: syntax=groovy
-*- mode: groovy;-*-
kate: syntax groovy; space-indent on; indent-width 2;
================================================================================
=                                 S  A  R  E  K                                =
================================================================================
 New Germline (+ Somatic) Analysis Workflow. Started March 2016.
--------------------------------------------------------------------------------
 @Authors
 Sebastian DiLorenzo <sebastian.dilorenzo@bils.se> [@Sebastian-D]
 Jesper Eisfeldt <jesper.eisfeldt@scilifelab.se> [@J35P312]
 Phil Ewels <phil.ewels@scilifelab.se> [@ewels]
 Maxime Garcia <maxime.garcia@scilifelab.se> [@MaxUlysse]
 Szilveszter Juhos <szilveszter.juhos@scilifelab.se> [@szilvajuhos]
 Max Käller <max.kaller@scilifelab.se> [@gulfshores]
 Malin Larsson <malin.larsson@scilifelab.se> [@malinlarsson]
 Marcel Martin <marcel.martin@scilifelab.se> [@marcelm]
 Björn Nystedt <bjorn.nystedt@scilifelab.se> [@bjornnystedt]
 Pall Olason <pall.olason@scilifelab.se> [@pallolason]
 Pelin Sahlén <pelin.akan@scilifelab.se> [@pelinakan]
--------------------------------------------------------------------------------
 @Homepage
 http://opensource.scilifelab.se/projects/sarek/
--------------------------------------------------------------------------------
 @Documentation
 https://github.com/SciLifeLab/Sarek/README.md
--------------------------------------------------------------------------------
 Processes overview
 - RunSamtoolsStats - Run Samtools stats on recalibrated BAM files
 - RunBamQC - Run qualimap BamQC on recalibrated BAM files
 - CreateIntervalBeds - Create and sort intervals into bed files
 - RunMutect1 - Run MuTect1 for Variant Calling (Parallelized processes)
 - RunMutect2 - Run MuTect2 for Variant Calling (Parallelized processes)
 - RunFreeBayes - Run FreeBayes for Variant Calling (Parallelized processes)
 - ConcatVCF - Merge results from MuTect1 and MuTect2
 - RunStrelka - Run Strelka for Variant Calling
 - RunManta - Run Manta for Structural Variant Calling
 - RunSingleManta - Run Manta for Single Structural Variant Calling
 - RunAlleleCount - Run AlleleCount to prepare for ASCAT
 - RunConvertAlleleCounts - Run convertAlleleCounts to prepare for ASCAT
 - RunAscat - Run ASCAT for CNV
 - RunBcftoolsStats - Run BCFTools stats on vcf files
================================================================================
=                           C O N F I G U R A T I O N                          =
================================================================================
*/

// Check that Nextflow version is up to date enough
// try / throw / catch works for NF versions < 0.25 when this was implemented
try {
    if( ! nextflow.version.matches(">= ${params.nfRequiredVersion}") ){
        throw GroovyException('Nextflow version too old')
    }
} catch (all) {
    log.error "====================================================\n" +
              "  Nextflow version ${params.nfRequiredVersion} required! You are running v${workflow.nextflow.version}.\n" +
              "  Pipeline execution will continue, but things may break.\n" +
              "  Please update Nextflow.\n" +
              "============================================================"
}

if (params.help) exit 0, helpMessage()
if (!SarekUtils.isAllowedParams(params)) exit 1, "params unknown, see --help for more information"
if (!checkUppmaxProject()) exit 1, "No UPPMAX project ID found! Use --project <UPPMAX Project ID>"

tools = params.tools ? params.tools.split(',').collect{it.trim().toLowerCase()} : []

directoryMap = defineDirectoryMap()
referenceMap = defineReferenceMap()
toolList = defineToolList()

if (!checkReferenceMap(referenceMap)) exit 1, 'Missing Reference file(s), see --help for more information'
if (!checkParameterList(tools,toolList)) exit 1, 'Unknown tool(s), see --help for more information'

if (params.test && params.genome in ['GRCh37', 'GRCh38']) {
  referenceMap.intervals = file("$workflow.projectDir/repeats/tiny_${params.genome}.list")
}

// TODO
// MuTect and Mutect2 could be run without a recalibrated BAM (they support
// the --BQSR option), but this is not implemented, yet.
// TODO
// FreeBayes does not need recalibrated BAMs, but we need to test whether
// the channels are set up correctly when we disable it

tsvPath = ''
if (params.sample) tsvPath = params.sample
else tsvPath = "${directoryMap.recalibrated}/recalibrated.tsv"

// Set up the bamFiles channel

bamFiles = Channel.empty()
if (tsvPath) {
  tsvFile = file(tsvPath)
  bamFiles = extractBams(tsvFile)
} else exit 1, 'No sample were defined, see --help'

(patientGenders, bamFiles) = extractGenders(bamFiles)

/*
================================================================================
=                               P R O C E S S E S                              =
================================================================================
*/

startMessage()

if (params.verbose) bamFiles = bamFiles.view {
  "BAMs to process:\n\
  ID    : ${it[0]}\tStatus: ${it[1]}\tSample: ${it[2]}\n\
  Files : [${it[3].fileName}, ${it[4].fileName}]"
}

// assume input is recalibrated, ignore explicitBqsrNeeded
(bamForBamQC, bamForSamToolsStats, recalibratedBam, recalTables) = bamFiles.into(4)

recalTables = recalTables.map{ it + [null] } // null recalibration table means: do not use --BQSR

recalTables = recalTables.map { [it[0]] + it[2..-1] } // remove status

if (params.verbose) recalibratedBam = recalibratedBam.view {
  "Recalibrated BAM for variant Calling:\n\
  ID    : ${it[0]}\tStatus: ${it[1]}\tSample: ${it[2]}\n\
  Files : [${it[3].fileName}, ${it[4].fileName}]"
}

process RunSamtoolsStats {
  tag {idPatient + "-" + idSample}

  publishDir directoryMap.samtoolsStats, mode: 'link'

  input:
    set idPatient, status, idSample, file(bam), file(bai) from bamForSamToolsStats

  output:
    file ("${bam}.samtools.stats.out") into samtoolsStatsReport

  when: !params.noReports

  script:
  """
  samtools stats ${bam} > ${bam}.samtools.stats.out
  """
}

if (params.verbose) samtoolsStatsReport = samtoolsStatsReport.view {
  "SAMTools stats report:\n\
  File  : [${it.fileName}]"
}

process RunBamQC {
  tag {idPatient + "-" + idSample}

  publishDir directoryMap.bamQC, mode: 'link'

  input:
    set idPatient, status, idSample, file(bam), file(bai) from bamForBamQC

  output:
    file("${idSample}") into bamQCreport

  when: !params.noReports && !params.noBAMQC

  script:
  """
  qualimap --java-mem-size=${task.memory.toGiga()}G \
  bamqc \
  -bam ${bam} \
  -outdir ${idSample} \
  -outformat HTML
  """
}

if (params.verbose) bamQCreport = bamQCreport.view {
  "BamQC report:\n\
  Dir   : [${it.fileName}]"
}

// Here we have a recalibrated bam set, but we need to separate the bam files based on patient status.
// The sample tsv config file which is formatted like: "subject status sample lane fastq1 fastq2"
// cf fastqFiles channel, I decided just to add _status to the sample name to have less changes to do.
// And so I'm sorting the channel if the sample match _0, then it's a normal sample, otherwise tumor.
// Then combine normal and tumor to get each possibilities
// ie. normal vs tumor1, normal vs tumor2, normal vs tumor3
// then copy this channel into channels for each variant calling
// I guess it will still work even if we have multiple normal samples

// separate recalibrateBams by status
bamsNormal = Channel.create()
bamsTumor = Channel.create()

recalibratedBam
  .choice(bamsTumor, bamsNormal) {it[1] == 0 ? 1 : 0}

bamsNormal = bamsNormal.ifEmpty{exit 1, "No normal sample defined, check TSV file: ${tsvFile}"}
bamsTumor = bamsTumor.ifEmpty{exit 1, "No tumor sample defined, check TSV file: ${tsvFile}"}

// Ascat, Strelka Germline & Manta Germline SV
bamsForAscat = Channel.create()
bamsForSingleManta = Channel.create()

(bamsTumorTemp, bamsTumor) = bamsTumor.into(2)
(bamsNormalTemp, bamsNormal) = bamsNormal.into(2)
(bamsForAscat, bamsForSingleManta) = bamsNormalTemp.mix(bamsTumorTemp).into(2)

// Removing status because not relevant anymore
bamsNormal = bamsNormal.map { idPatient, status, idSample, bam, bai -> [idPatient, idSample, bam, bai] }

bamsTumor = bamsTumor.map { idPatient, status, idSample, bam, bai -> [idPatient, idSample, bam, bai] }

// We know that MuTect2 (and other somatic callers) are notoriously slow.
// To speed them up we are chopping the reference into smaller pieces.
// (see repeats/centromeres.list).
// Do variant calling by this intervals, and re-merge the VCFs.
// Since we are on a cluster, this can parallelize the variant call processes.
// And push down the variant call wall clock time significanlty.

process CreateIntervalBeds {
  tag {intervals.fileName}

  input:
    file(intervals) from Channel.value(referenceMap.intervals)

  output:
    file '*.bed' into bedIntervals mode flatten

  script:
  // If the interval file is BED format, the fifth column is interpreted to
  // contain runtime estimates, which is then used to combine short-running jobs
  if (intervals.getName().endsWith('.bed'))
    """
    awk -vFS="\t" '{
      t = \$5  # runtime estimate
      if (t == "") {
        # no runtime estimate in this row, assume default value
        t = (\$3 - \$2) / ${params.nucleotidesPerSecond}
      }
      if (name == "" || (chunk > 600 && (chunk + t) > longest * 1.05)) {
        # start a new chunk
        name = sprintf("%s_%d-%d.bed", \$1, \$2+1, \$3)
        chunk = 0
        longest = 0
      }
      if (t > longest)
        longest = t
      chunk += t
      print \$0 > name
    }' ${intervals}
    """
  else
    """
    awk -vFS="[:-]" '{
      name = sprintf("%s_%d-%d", \$1, \$2, \$3);
      printf("%s\\t%d\\t%d\\n", \$1, \$2-1, \$3) > name ".bed"
    }' ${intervals}
    """
}

bedIntervals = bedIntervals
  .map { intervalFile ->
    final duration = 0.0
    for (line in intervalFile.readLines()) {
      final fields = line.split('\t')
      if (fields.size() >= 5) duration += fields[4].toFloat()
      else {
        start = fields[1].toInteger()
        end = fields[2].toInteger()
        duration += (end - start) / params.nucleotidesPerSecond
      }
    }
    [duration, intervalFile]
  }.toSortedList({ a, b -> b[0] <=> a[0] })
  .flatten().collate(2)
  .map{duration, intervalFile -> intervalFile}

if (params.verbose) bedIntervals = bedIntervals.view {
  "  Interv: ${it.baseName}"
}

(bamsNormalTemp, bamsNormal, bedIntervals) = generateIntervalsForVC(bamsNormal, bedIntervals)
(bamsTumorTemp, bamsTumor, bedIntervals) = generateIntervalsForVC(bamsTumor, bedIntervals)

bamsAll = bamsNormal.combine(bamsTumor)

// Since idPatientNormal and idPatientTumor are the same
// It's removed from bamsAll Channel (same for genderNormal)
// /!\ It is assumed that every sample are from the same patient
bamsAll = bamsAll.map {
  idPatientNormal, idSampleNormal, bamNormal, baiNormal, idPatientTumor, idSampleTumor, bamTumor, baiTumor ->
  [idPatientNormal, idSampleNormal, bamNormal, baiNormal, idSampleTumor, bamTumor, baiTumor]
}

// Manta and Strelka
(bamsForManta, bamsForStrelka, bamsAll) = bamsAll.into(3)

bamsTumorNormalIntervals = bamsAll.spread(bedIntervals)

// MuTect1, MuTect2, FreeBayes
(bamsFMT1, bamsFMT2, bamsFFB) = bamsTumorNormalIntervals.into(3)

process RunMutect1 {
  tag {idSampleTumor + "_vs_" + idSampleNormal + "-" + intervalBed.baseName}

  input:
    set idPatient, idSampleNormal, file(bamNormal), file(baiNormal), idSampleTumor, file(bamTumor), file(baiTumor), file(intervalBed) from bamsFMT1
    set file(genomeFile), file(genomeIndex), file(genomeDict), file(dbsnp), file(dbsnpIndex), file(cosmic), file(cosmicIndex) from Channel.value([
      referenceMap.genomeFile,
      referenceMap.genomeIndex,
      referenceMap.genomeDict,
      referenceMap.dbsnp,
      referenceMap.dbsnpIndex,
      referenceMap.cosmic,
      referenceMap.cosmicIndex
    ])

  output:
    set val("mutect1"), idPatient, idSampleNormal, idSampleTumor, file("${intervalBed.baseName}_${idSampleTumor}_vs_${idSampleNormal}.vcf") into mutect1Output

  when: 'mutect1' in tools && !params.onlyQC

  script:
  """
  java -Xmx${task.memory.toGiga()}g \
  -jar \$MUTECT_HOME/muTect.jar \
  -T MuTect \
  -R ${genomeFile} \
  --cosmic ${cosmic} \
  --dbsnp ${dbsnp} \
  -I:normal ${bamNormal} \
  -I:tumor ${bamTumor} \
  -L ${intervalBed} \
  --disable_auto_index_creation_and_locking_when_reading_rods \
  --out ${intervalBed.baseName}_${idSampleTumor}_vs_${idSampleNormal}.call_stats.out \
  --vcf ${intervalBed.baseName}_${idSampleTumor}_vs_${idSampleNormal}.vcf
  """
}

mutect1Output = mutect1Output.groupTuple(by:[0,1,2,3])

process RunMutect2 {
  tag {idSampleTumor + "_vs_" + idSampleNormal + "-" + intervalBed.baseName}

  input:
    set idPatient, idSampleNormal, file(bamNormal), file(baiNormal), idSampleTumor, file(bamTumor), file(baiTumor), file(intervalBed) from bamsFMT2
    set file(genomeFile), file(genomeIndex), file(genomeDict), file(dbsnp), file(dbsnpIndex), file(cosmic), file(cosmicIndex) from Channel.value([
      referenceMap.genomeFile,
      referenceMap.genomeIndex,
      referenceMap.genomeDict,
      referenceMap.dbsnp,
      referenceMap.dbsnpIndex,
      referenceMap.cosmic,
      referenceMap.cosmicIndex
    ])

  output:
    set val("mutect2"), idPatient, idSampleNormal, idSampleTumor, file("${intervalBed.baseName}_${idSampleTumor}_vs_${idSampleNormal}.vcf") into mutect2Output

  when: 'mutect2' in tools && !params.onlyQC

  script:
  """
  java -Xmx${task.memory.toGiga()}g \
  -jar \$GATK_HOME/GenomeAnalysisTK.jar \
  -T MuTect2 \
  -R ${genomeFile} \
  --cosmic ${cosmic} \
  --dbsnp ${dbsnp} \
  -I:normal ${bamNormal} \
  -I:tumor ${bamTumor} \
  --disable_auto_index_creation_and_locking_when_reading_rods \
  -L ${intervalBed} \
  -o ${intervalBed.baseName}_${idSampleTumor}_vs_${idSampleNormal}.vcf
  """
}

mutect2Output = mutect2Output.groupTuple(by:[0,1,2,3])

process RunFreeBayes {
  tag {idSampleTumor + "_vs_" + idSampleNormal + "-" + intervalBed.baseName}

  input:
    set idPatient, idSampleNormal, file(bamNormal), file(baiNormal), idSampleTumor, file(bamTumor), file(baiTumor), file(intervalBed) from bamsFFB
    file(genomeFile) from Channel.value(referenceMap.genomeFile)

  output:
    set val("freebayes"), idPatient, idSampleNormal, idSampleTumor, file("${intervalBed.baseName}_${idSampleTumor}_vs_${idSampleNormal}.vcf") into freebayesOutput

  when: 'freebayes' in tools && !params.onlyQC

  script:
  """
  freebayes \
    -f ${genomeFile} \
    --pooled-continuous \
    --pooled-discrete \
    --genotype-qualities \
    --report-genotype-likelihood-max \
    --allele-balance-priors-off \
    --min-alternate-fraction 0.03 \
    --min-repeat-entropy 1 \
    --min-alternate-count 2 \
    -t ${intervalBed} \
    ${bamTumor} \
    ${bamNormal} > ${intervalBed.baseName}_${idSampleTumor}_vs_${idSampleNormal}.vcf
  """
}

freebayesOutput = freebayesOutput.groupTuple(by:[0,1,2,3])

// we are merging the VCFs that are called separatelly for different intervals
// so we can have a single sorted VCF containing all the calls for a given caller

vcfsToMerge = mutect1Output.mix(mutect2Output, freebayesOutput)
if (params.verbose) vcfsToMerge = vcfsToMerge.view {
  "VCFs To be merged:\n\
  Tool  : ${it[0]}\tID    : ${it[1]}\tSample: [${it[3]}, ${it[2]}]\n\
  Files : ${it[4].fileName}"
}

process ConcatVCF {
  tag {variantCaller + "_" + idSampleTumor + "_vs_" + idSampleNormal}

  publishDir "${directoryMap."$variantCaller"}", mode: 'link'

  input:
    set variantCaller, idPatient, idSampleNormal, idSampleTumor, file(vcFiles) from vcfsToMerge
    file(genomeIndex) from Channel.value(referenceMap.genomeIndex)

  output:
    set variantCaller, idPatient, idSampleNormal, idSampleTumor, file("*.vcf.gz") into vcfConcatenated
    file("*.vcf.gz.tbi") into vcfConcatenatedTbi

  when: ('mutect1' in tools || 'mutect2' in tools || 'freebayes' in tools ) && !params.onlyQC

  script:
  outputFile = "${variantCaller}_${idSampleTumor}_vs_${idSampleNormal}.vcf"

  """
  # first make a header from one of the VCF intervals
  # get rid of interval information only from the GATK command-line, but leave the rest
  FIRSTVCF=\$(ls *.vcf | head -n 1)
  sed -n '/^[^#]/q;p' \$FIRSTVCF | \
  awk '!/GATKCommandLine/{print}/GATKCommandLine/{for(i=1;i<=NF;i++){if(\$i!~/intervals=/ && \$i !~ /out=/){printf("%s ",\$i)}}printf("\\n")}' \
  > header

  # Get list of contigs from the FASTA index (.fai). We cannot use the ##contig
  # header in the VCF as it is optional (FreeBayes does not save it, for example)
  CONTIGS=(\$(cut -f1 ${genomeIndex}))

  # concatenate VCFs in the correct order
  (
    cat header

    for chr in "\${CONTIGS[@]}"; do
      # Skip if globbing would not match any file to avoid errors such as
      # "ls: cannot access chr3_*.vcf: No such file or directory" when chr3
      # was not processed.
      pattern="\${chr}_*.vcf"
      if ! compgen -G "\${pattern}" > /dev/null; then continue; fi

      # ls -v sorts by numeric value ("version"), which means that chr1_100_
      # is sorted *after* chr1_99_.
      for vcf in \$(ls -v \${pattern}); do
        # Determine length of header.
        # The 'q' command makes sed exit when it sees the first non-header
        # line, which avoids reading in the entire file.
        L=\$(sed -n '/^[^#]/q;p' \${vcf} | wc -l)

        # Then print all non-header lines. Since tail is very fast (nearly as
        # fast as cat), this is way more efficient than using a single sed,
        # awk or grep command.
        tail -n +\$((L+1)) \${vcf}
      done
    done
  ) | bgzip > ${outputFile}.gz
  tabix ${outputFile}.gz
  """
}

if (params.verbose) vcfConcatenated = vcfConcatenated.view {
  "Variant Calling output:\n\
  Tool  : ${it[0]}\tID    : ${it[1]}\tSample: [${it[3]}, ${it[2]}]\n\
  File  : ${it[4].fileName}"
}

process RunStrelka {
  tag {idSampleTumor + "_vs_" + idSampleNormal}

  publishDir directoryMap.strelka, mode: 'link'

  input:
    set idPatient, idSampleNormal, file(bamNormal), file(baiNormal), idSampleTumor, file(bamTumor), file(baiTumor) from bamsForStrelka
    set file(genomeFile), file(genomeIndex), file(genomeDict) from Channel.value([
      referenceMap.genomeFile,
      referenceMap.genomeIndex,
      referenceMap.genomeDict
    ])

  output:
    set val("strelka"), idPatient, idSampleNormal, idSampleTumor, file("*.vcf.gz"), file("*.vcf.gz.tbi") into strelkaOutput

  when: 'strelka' in tools && !params.onlyQC

  script:
  """
  \$STRELKA_INSTALL_PATH/bin/configureStrelkaSomaticWorkflow.py \
  --tumor ${bamTumor} \
  --normal ${bamNormal} \
  --referenceFasta ${genomeFile} \
  --runDir Strelka

  python Strelka/runWorkflow.py -m local -j ${task.cpus}

  mv Strelka/results/variants/somatic.indels.vcf.gz \
    Strelka_${idSampleTumor}_vs_${idSampleNormal}_somatic_indels.vcf.gz
  mv Strelka/results/variants/somatic.indels.vcf.gz.tbi \
    Strelka_${idSampleTumor}_vs_${idSampleNormal}_somatic_indels.vcf.gz.tbi
  mv Strelka/results/variants/somatic.snvs.vcf.gz \
    Strelka_${idSampleTumor}_vs_${idSampleNormal}_somatic_snvs.vcf.gz
  mv Strelka/results/variants/somatic.snvs.vcf.gz.tbi \
    Strelka_${idSampleTumor}_vs_${idSampleNormal}_somatic_snvs.vcf.gz.tbi
  """
}

if (params.verbose) strelkaOutput = strelkaOutput.view {
  "Variant Calling output:\n\
  Tool  : ${it[0]}\tID    : ${it[1]}\tSample: [${it[3]}, ${it[2]}]\n\
  Files : ${it[4].fileName}\n\
  Index : ${it[5].fileName}"
}

process RunManta {
  tag {idSampleTumor + "_vs_" + idSampleNormal}

  publishDir directoryMap.manta, mode: 'link'

  input:
    set idPatient, idSampleNormal, file(bamNormal), file(baiNormal), idSampleTumor, file(bamTumor), file(baiTumor) from bamsForManta
    set file(genomeFile), file(genomeIndex) from Channel.value([
      referenceMap.genomeFile,
      referenceMap.genomeIndex
    ])

  output:
    set val("manta"), idPatient, idSampleNormal, idSampleTumor, file("*.vcf.gz"), file("*.vcf.gz.tbi") into mantaOutput
    set idPatient, idSampleNormal, file(bamNormal), file(baiNormal), idSampleTumor, file(bamTumor), file(baiTumor), file("*.candidateSmallIndels.vcf.gz"), file("*.candidateSmallIndels.vcf.gz.tbi") into mantaToStrelka

  when: 'manta' in tools && !params.onlyQC

  script:
  """
  \$MANTA_INSTALL_PATH/bin/configManta.py \
  --normalBam ${bamNormal} \
  --tumorBam ${bamTumor} \
  --reference ${genomeFile} \
  --runDir Manta

  python Manta/runWorkflow.py -m local -j ${task.cpus}

  mv Manta/results/variants/candidateSmallIndels.vcf.gz \
    Manta_${idSampleTumor}_vs_${idSampleNormal}.candidateSmallIndels.vcf.gz
  mv Manta/results/variants/candidateSmallIndels.vcf.gz.tbi \
    Manta_${idSampleTumor}_vs_${idSampleNormal}.candidateSmallIndels.vcf.gz.tbi
  mv Manta/results/variants/candidateSV.vcf.gz \
    Manta_${idSampleTumor}_vs_${idSampleNormal}.candidateSV.vcf.gz
  mv Manta/results/variants/candidateSV.vcf.gz.tbi \
    Manta_${idSampleTumor}_vs_${idSampleNormal}.candidateSV.vcf.gz.tbi
  mv Manta/results/variants/diploidSV.vcf.gz \
    Manta_${idSampleTumor}_vs_${idSampleNormal}.diploidSV.vcf.gz
  mv Manta/results/variants/diploidSV.vcf.gz.tbi \
    Manta_${idSampleTumor}_vs_${idSampleNormal}.diploidSV.vcf.gz.tbi
  mv Manta/results/variants/somaticSV.vcf.gz \
    Manta_${idSampleTumor}_vs_${idSampleNormal}.somaticSV.vcf.gz
  mv Manta/results/variants/somaticSV.vcf.gz.tbi \
    Manta_${idSampleTumor}_vs_${idSampleNormal}.somaticSV.vcf.gz.tbi
  """
}

if (params.verbose) mantaOutput = mantaOutput.view {
  "Variant Calling output:\n\
  Tool  : ${it[0]}\tID    : ${it[1]}\tSample: [${it[3]}, ${it[2]}]\n\
  Files : ${it[4].fileName}\n\
  Index : ${it[5].fileName}"
}

process RunSingleManta {
  tag {idSample + " - Tumor-Only"}

  publishDir directoryMap.manta, mode: 'link'

  input:
    set idPatient, status, idSample, file(bam), file(bai) from bamsForSingleManta
    set file(genomeFile), file(genomeIndex) from Channel.value([
      referenceMap.genomeFile,
      referenceMap.genomeIndex
    ])

  output:
    set val("singlemanta"), idPatient, idSample,  file("*.vcf.gz"), file("*.vcf.gz.tbi") into singleMantaOutput

  when: 'manta' in tools && status == 1 && !params.onlyQC

  script:
  """
  \$MANTA_INSTALL_PATH/bin/configManta.py \
  --tumorBam ${bam} \
  --reference ${genomeFile} \
  --runDir Manta

  python Manta/runWorkflow.py -m local -j ${task.cpus}

  mv Manta/results/variants/candidateSmallIndels.vcf.gz \
    Manta_${idSample}.candidateSmallIndels.vcf.gz
  mv Manta/results/variants/candidateSmallIndels.vcf.gz.tbi \
    Manta_${idSample}.candidateSmallIndels.vcf.gz.tbi
  mv Manta/results/variants/candidateSV.vcf.gz \
    Manta_${idSample}.candidateSV.vcf.gz
  mv Manta/results/variants/candidateSV.vcf.gz.tbi \
    Manta_${idSample}.candidateSV.vcf.gz.tbi
  mv Manta/results/variants/tumorSV.vcf.gz \
    Manta_${idSample}.tumorSV.vcf.gz
  mv Manta/results/variants/tumorSV.vcf.gz.tbi \
    Manta_${idSample}.tumorSV.vcf.gz.tbi
  """
}

if (params.verbose) singleMantaOutput = singleMantaOutput.view {
  "Variant Calling output:\n\
  Tool  : ${it[0]}\tID    : ${it[1]}\tSample: ${it[2]}\n\
  Files : ${it[3].fileName}\n\
  Index : ${it[4].fileName}"
}


process RunStrelkaBP {
  tag {idSampleTumor + "_vs_" + idSampleNormal}

  publishDir directoryMap.strelkabp, mode: 'link'

  input:
    set idPatient, idSampleNormal, file(bamNormal), file(baiNormal), idSampleTumor, file(bamTumor), file(baiTumor), file(mantaCSI), file(mantaCSIi) from mantaToStrelka
    set file(genomeFile), file(genomeIndex), file(genomeDict) from Channel.value([
      referenceMap.genomeFile,
      referenceMap.genomeIndex,
      referenceMap.genomeDict
    ])

  output:
    set val("strelkaBP"), idPatient, idSampleNormal, idSampleTumor, file("*.vcf.gz"), file("*.vcf.gz.tbi") into strelkaBPOutput

  when: 'strelka' in tools && 'manta' in tools && params.strelkaBP && !params.onlyQC

  script:
  """
  \$STRELKA_INSTALL_PATH/bin/configureStrelkaSomaticWorkflow.py \
  --tumor ${bamTumor} \
  --normal ${bamNormal} \
  --referenceFasta ${genomeFile} \
  --indelCandidates ${mantaCSI} \
  --runDir Strelka

  python Strelka/runWorkflow.py -m local -j ${task.cpus}

  mv Strelka/results/variants/somatic.indels.vcf.gz \
    Strelka_${idSampleTumor}_vs_${idSampleNormal}_somatic_indels.vcf.gz
  mv Strelka/results/variants/somatic.indels.vcf.gz.tbi \
    Strelka_${idSampleTumor}_vs_${idSampleNormal}_somatic_indels.vcf.gz.tbi
  mv Strelka/results/variants/somatic.snvs.vcf.gz \
    Strelka_${idSampleTumor}_vs_${idSampleNormal}_somatic_snvs.vcf.gz
  mv Strelka/results/variants/somatic.snvs.vcf.gz.tbi \
    Strelka_${idSampleTumor}_vs_${idSampleNormal}_somatic_snvs.vcf.gz.tbi
  """
}

if (params.verbose) strelkaBPOutput = strelkaBPOutput.view {
  "Variant Calling output:\n\
  Tool  : ${it[0]}\tID    : ${it[1]}\tSample: [${it[3]}, ${it[2]}]\n\
  Files : ${it[4].fileName}\n\
  Index : ${it[5].fileName}"
}

// Run commands and code from Malin Larsson
// Based on Jesper Eisfeldt's code
process RunAlleleCount {
  tag {idSample}

  input:
    set idPatient, status, idSample, file(bam), file(bai) from bamsForAscat
    set file(acLoci), file(genomeFile), file(genomeIndex), file(genomeDict) from Channel.value([
      referenceMap.acLoci,
      referenceMap.genomeFile,
      referenceMap.genomeIndex,
      referenceMap.genomeDict
    ])

  output:
    set idPatient, status, idSample, file("${idSample}.alleleCount") into alleleCountOutput

  when: 'ascat' in tools && !params.onlyQC

  script:
  """
  alleleCounter \
  -l ${acLoci} \
  -r ${genomeFile} \
  -b ${bam} \
  -o ${idSample}.alleleCount;
  """
}

alleleCountNormal = Channel.create()
alleleCountTumor = Channel.create()

alleleCountOutput
  .choice(alleleCountTumor, alleleCountNormal) {it[1] == 0 ? 1 : 0}

alleleCountOutput = alleleCountNormal.combine(alleleCountTumor)

alleleCountOutput = alleleCountOutput.map {
  idPatientNormal, statusNormal, idSampleNormal, alleleCountNormal,
  idPatientTumor,  statusTumor,  idSampleTumor,  alleleCountTumor ->
  [idPatientNormal, idSampleNormal, idSampleTumor, alleleCountNormal, alleleCountTumor]
}

// R script from Malin Larssons bitbucket repo:
// https://bitbucket.org/malinlarsson/somatic_wgs_pipeline
process RunConvertAlleleCounts {
  tag {idSampleTumor + "_vs_" + idSampleNormal}

  publishDir directoryMap.ascat, mode: 'link'

  input:
    set idPatient, idSampleNormal, idSampleTumor, file(alleleCountNormal), file(alleleCountTumor) from alleleCountOutput

  output:
    set idPatient, idSampleNormal, idSampleTumor, file("${idSampleNormal}.BAF"), file("${idSampleNormal}.LogR"), file("${idSampleTumor}.BAF"), file("${idSampleTumor}.LogR") into convertAlleleCountsOutput

  when: 'ascat' in tools && !params.onlyQC

  script:
  gender = patientGenders[idPatient]
  """
  convertAlleleCounts.r ${idSampleTumor} ${alleleCountTumor} ${idSampleNormal} ${alleleCountNormal} ${gender}
  """
}

// R scripts from Malin Larssons bitbucket repo:
// https://bitbucket.org/malinlarsson/somatic_wgs_pipeline
process RunAscat {
  tag {idSampleTumor + "_vs_" + idSampleNormal}

  publishDir directoryMap.ascat, mode: 'link'

  input:
    set idPatient, idSampleNormal, idSampleTumor, file(bafNormal), file(logrNormal), file(bafTumor), file(logrTumor) from convertAlleleCountsOutput

  output:
    set val("ascat"), idPatient, idSampleNormal, idSampleTumor, file("${idSampleTumor}.*.{png,txt}") into ascatOutput

  when: 'ascat' in tools && !params.onlyQC

  script:
  """
  # get rid of "chr" string if there is any
  for f in *BAF *LogR; do sed 's/chr//g' \$f > tmpFile; mv tmpFile \$f;done
  run_ascat.r ${bafTumor} ${logrTumor} ${bafNormal} ${logrNormal} ${idSampleTumor} ${baseDir}
  """
}

if (params.verbose) ascatOutput = ascatOutput.view {
  "Variant Calling output:\n\
  Tool  : ${it[0]}\tID    : ${it[1]}\tSample: [${it[3]}, ${it[2]}]\n\
  Files : [${it[4].fileName}]"
}

(strelkaIndels, strelkaSNVS) = strelkaOutput.into(2)
(mantaSomaticSV, mantaDiploidSV) = mantaOutput.into(2)

vcfForBCFtools = Channel.empty().mix(
  mantaDiploidSV.map {
    variantcaller, idPatient, idSampleNormal, idSampleTumor, vcf, tbi ->
    [variantcaller, vcf[2]]
  },
  mantaSomaticSV.map {
    variantcaller, idPatient, idSampleNormal, idSampleTumor, vcf, tbi ->
    [variantcaller, vcf[3]]
  },
  singleMantaOutput.map {
    variantcaller, idPatient, idSample, vcf, tbi ->
    [variantcaller, vcf[2]]
  },
  strelkaIndels.map {
    variantcaller, idPatient, idSampleNormal, idSampleTumor, vcf, tbi ->
    [variantcaller, vcf[0]]
  },
  strelkaSNVS.map {
    variantcaller, idPatient, idSampleNormal, idSampleTumor, vcf, tbi ->
    [variantcaller, vcf[1]]
  })

process RunBcftoolsStats {
  tag {vcf}

  publishDir directoryMap.bcftoolsStats, mode: 'link'

  input:
    set variantCaller, file(vcf) from vcfForBCFtools

  output:
    file ("${vcf.baseName}.bcf.tools.stats.out") into bcfReport

  when: !params.noReports

  script:
  """
  bcftools stats ${vcf} > ${vcf.baseName}.bcf.tools.stats.out
  """
}

if (params.verbose) bcfReport = bcfReport.view {
  "BCFTools stats report:\n\
  File  : [${it.fileName}]"
}

bcfReport.close()

/*
================================================================================
=                               F U N C T I O N S                              =
================================================================================
*/

def checkFileExtension(it, extension) {
  // Check file extension
  if (!it.toString().toLowerCase().endsWith(extension.toLowerCase())) exit 1, "File: ${it} has the wrong extension: ${extension} see --help for more information"
}

def checkParameterExistence(it, list) {
  // Check parameter existence
  if (!list.contains(it)) {
    println("Unknown parameter: ${it}")
    return false
  }
  return true
}

def checkParameterList(list, realList) {
  // Loop through all parameters to check their existence and spelling
  return list.every{ checkParameterExistence(it, realList) }
}

def checkParamReturnFile(item) {
  params."${item}" = params.genomes[params.genome]."${item}"
  return file(params."${item}")
}

def checkReferenceMap(referenceMap) {
  // Loop through all the references files to check their existence
  referenceMap.every {
    referenceFile, fileToCheck ->
    checkRefExistence(referenceFile, fileToCheck)
  }
}

def checkRefExistence(referenceFile, fileToCheck) {
  if (fileToCheck instanceof List) return fileToCheck.every{ checkRefExistence(referenceFile, it) }
  def f = file(fileToCheck)
  // this is an expanded wildcard: we can assume all files exist
  if (f instanceof List && f.size() > 0) return true
  else if (!f.exists()) {
    log.info  "Missing references: ${referenceFile} ${fileToCheck}"
    return false
  }
  return true
}

def checkUppmaxProject() {
  // check if UPPMAX project number is specified
  return !(workflow.profile == 'slurm' && !params.project)
}

def defineDirectoryMap() {
  return [
    'recalibrated'     : "${params.outDir}/Preprocessing/Recalibrated",
    'bamQC'            : "${params.outDir}/Reports/bamQC",
    'bcftoolsStats'    : "${params.outDir}/Reports/BCFToolsStats",
    'samtoolsStats'    : "${params.outDir}/Reports/SamToolsStats",
    'ascat'            : "${params.outDir}/VariantCalling/Ascat",
    'freebayes'        : "${params.outDir}/VariantCalling/FreeBayes",
    'manta'            : "${params.outDir}/VariantCalling/Manta",
    'mutect1'          : "${params.outDir}/VariantCalling/MuTect1",
    'mutect2'          : "${params.outDir}/VariantCalling/MuTect2",
    'strelka'          : "${params.outDir}/VariantCalling/Strelka",
    'strelkabp'        : "${params.outDir}/VariantCalling/StrelkaBP"
  ]
}

def defineReferenceMap() {
  if (!(params.genome in params.genomes)) exit 1, "Genome ${params.genome} not found in configuration"
  return [
    // loci file for ascat
    'acLoci'           : checkParamReturnFile("acLoci"),
    'dbsnp'            : checkParamReturnFile("dbsnp"),
    'dbsnpIndex'       : checkParamReturnFile("dbsnpIndex"),
    // cosmic VCF with VCF4.1 header
    'cosmic'           : checkParamReturnFile("cosmic"),
    'cosmicIndex'      : checkParamReturnFile("cosmicIndex"),
    // genome reference dictionary
    'genomeDict'       : checkParamReturnFile("genomeDict"),
    // FASTA genome reference
    'genomeFile'       : checkParamReturnFile("genomeFile"),
    // genome .fai file
    'genomeIndex'      : checkParamReturnFile("genomeIndex"),
    // intervals file for spread-and-gather processes
    'intervals'        : checkParamReturnFile("intervals")
  ]
}

def defineToolList() {
  return [
    'ascat',
    'freebayes',
    'haplotypecaller',
    'manta',
    'mutect1',
    'mutect2',
    'strelka'
  ]
}

def extractBams(tsvFile) {
  // Channeling the TSV file containing BAM.
  // Format is: "subject gender status sample bam bai"
  Channel
    .from(tsvFile.readLines())
    .map{line ->
      def list      = returnTSV(line.split(),6)
      def idPatient = list[0]
      def gender    = list[1]
      def status    = returnStatus(list[2].toInteger())
      def idSample  = list[3]
      def bamFile   = returnFile(list[4])
      def baiFile   = returnFile(list[5])

      checkFileExtension(bamFile,".bam")
      checkFileExtension(baiFile,".bai")

      [ idPatient, gender, status, idSample, bamFile, baiFile ]
    }
}

def extractGenders(channel) {
  def genders = [:]  // an empty map
  channel = channel.map{ it ->
    def idPatient = it[0]
    def gender = it[1]
    genders[idPatient] = gender

    [idPatient] + it[2..-1]
  }
  [genders, channel]
}

def generateIntervalsForVC(bams, intervals) {

  def (bamsNew, bamsForVC) = bams.into(2)
  def (intervalsNew, vcIntervals) = intervals.into(2)
  def bamsForVCNew = bamsForVC.combine(vcIntervals)
  return [bamsForVCNew, bamsNew, intervalsNew]
}

def grabRevision() {
  // Return the same string executed from github or not
  return workflow.revision ?: workflow.commitId ?: workflow.scriptId.substring(0,10)
}

def helpMessage() {
  // Display help message
  this.sarekMessage()
  log.info "    Usage:"
  log.info "       nextflow run somaticVC.nf --sample <file.tsv> [--tools TOOL[,TOOL]] --genome <Genome>"
  log.info "    --sample <file.tsv>"
  log.info "       Specify a TSV file containing paths to sample files."
  log.info "    --test"
  log.info "       Use a test sample."
  log.info "    --noReports"
  log.info "       Disable QC tools and MultiQC to generate a HTML report"
  log.info "    --tools"
  log.info "       Option to configure which tools to use in the workflow."
  log.info "         Different tools to be separated by commas."
  log.info "       Possible values are:"
  log.info "         mutect1 (use MuTect1 for VC)"
  log.info "         mutect2 (use MuTect2 for VC)"
  log.info "         freebayes (use FreeBayes for VC)"
  log.info "         strelka (use Strelka for VC)"
  log.info "         haplotypecaller (use HaplotypeCaller for normal bams VC)"
  log.info "         manta (use Manta for SV)"
  log.info "         ascat (use Ascat for CNV)"
  log.info "    --genome <Genome>"
  log.info "       Use a specific genome version."
  log.info "       Possible values are:"
  log.info "         GRCh37"
  log.info "         GRCh38 (Default)"
  log.info "         smallGRCh37 (Use a small reference (Tests only))"
  log.info "    --onlyQC"
  log.info "       Run only QC tools and gather reports"
  log.info "    --help"
  log.info "       you're reading it"
  log.info "    --verbose"
  log.info "       Adds more verbosity to workflow"
}

def minimalInformationMessage() {
  // Minimal information message
  log.info "Command Line: " + workflow.commandLine
  log.info "Profile     : " + workflow.profile
  log.info "Project Dir : " + workflow.projectDir
  log.info "Launch Dir  : " + workflow.launchDir
  log.info "Work Dir    : " + workflow.workDir
  log.info "Out Dir     : " + params.outDir
  log.info "TSV file    : ${tsvFile}"
  log.info "Genome      : " + params.genome
  log.info "Genome_base : " + params.genome_base
  log.info "Tools       : " + tools.join(', ')
  log.info "Containers  :"
  if (params.repository) log.info "  Repository   : ${params.repository}"
  else log.info "  ContainerPath: " + params.containerPath
  log.info "  Tag          : " + params.tag
  log.info "Reference files used:"
  log.info "  acLoci      :\n\t" + referenceMap.acLoci
  log.info "  cosmic      :\n\t" + referenceMap.cosmic
  log.info "\t" + referenceMap.cosmicIndex
  log.info "  dbsnp       :\n\t" + referenceMap.dbsnp
  log.info "\t" + referenceMap.dbsnpIndex
  log.info "  genome      :\n\t" + referenceMap.genomeFile
  log.info "\t" + referenceMap.genomeDict
  log.info "\t" + referenceMap.genomeIndex
  log.info "  intervals   :\n\t" + referenceMap.intervals
}

def nextflowMessage() {
  // Nextflow message (version + build)
  log.info "N E X T F L O W  ~  version ${workflow.nextflow.version} ${workflow.nextflow.build}"
}

def returnFile(it) {
  // return file if it exists
  if (!file(it).exists()) exit 1, "Missing file in TSV file: ${it}, see --help for more information"
  return file(it)
}

def returnStatus(it) {
  // Return status if it's correct
  // Status should be only 0 or 1
  // 0 being normal
  // 1 being tumor (or relapse or anything that is not normal...)
  if (!(it in [0, 1])) exit 1, "Status is not recognized in TSV file: ${it}, see --help for more information"
  return it
}

def returnTSV(it, number) {
  // return TSV if it has the correct number of items in row
  if (it.size() != number) exit 1, "Malformed row in TSV file: ${it}, see --help for more information"
  return it
}

def sarekMessage() {
  // Display Sarek message
  log.info "Sarek - Workflow For Somatic And Germline Variations ~ ${params.version} - " + this.grabRevision() + (workflow.commitId ? " [${workflow.commitId}]" : "")
}

def startMessage() {
  // Display start message
  this.sarekMessage()
  this.minimalInformationMessage()
}

workflow.onComplete {
  // Display complete message
  this.nextflowMessage()
  this.sarekMessage()
  this.minimalInformationMessage()
  log.info "Completed at: " + workflow.complete
  log.info "Duration    : " + workflow.duration
  log.info "Success     : " + workflow.success
  log.info "Exit status : " + workflow.exitStatus
  log.info "Error report: " + (workflow.errorReport ?: '-')
}

workflow.onError {
  // Display error message
  this.nextflowMessage()
  this.sarekMessage()
  log.info "Workflow execution stopped with the following message:"
  log.info "  " + workflow.errorMessage
}
