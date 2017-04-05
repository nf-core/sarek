#!/usr/bin/env nextflow

/*
vim: syntax=groovy
-*- mode: groovy;-*-
================================================================================
=               C A N C E R    A N A L Y S I S    W O R K F L O W              =
================================================================================
 New Cancer Analysis Workflow. Started March 2016.
--------------------------------------------------------------------------------
 @Authors
 Sebastian DiLorenzo <sebastian.dilorenzo@bils.se> [@Sebastian-D]
 Jesper Eisfeldt <jesper.eisfeldt@scilifelab.se> [@J35P312]
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
 http://opensource.scilifelab.se/projects/caw/
--------------------------------------------------------------------------------
*/

/*
================================================================================
=                               P R O C E S S E S                              =
================================================================================
*/
if (params.genome == "smallGRCh37") {
  referencesToDownload =
    Channel.from([
      '1000G_phase1.indels.b37.small.vcf.gz',
      '1000G_phase3_20130502_SNP_maf0.3.small.loci',
      'b37_cosmic_v54_120711.small.vcf.gz',
      'b37_cosmic_v74.noCHR.sort.4.1.small.vcf.gz',
      'dbsnp_138.b37.small.vcf.gz',
      'human_g1k_v37_decoy.small.fasta.gz',
      'Mills_and_1000G_gold_standard.indels.b37.small.vcf.gz',
      'small.intervals'
    ])
} else {exit 1, "Can't build this reference genome"}

process DownloadReference {
  tag {reference}

  input:
    val(reference) from referencesToDownload

  output:
    file(reference) into downloadedFiles

  script:

  """
  set -euo pipefail
  wget https://github.com/MaxUlysse/smallRef/raw/master/$reference
  """
}

gzFiles = Channel.create()
notGZfiles = Channel.create()

downloadedFiles
  .choice(gzFiles, notGZfiles) {it =~ ".gz" ? 0 : 1}

notGZfiles.collectFile(storeDir: "References/" + params.genome)

process DecompressFile {
  tag {reference}

  input:
    file(reference) from gzFiles

  output:
    file("*.{vcf,fasta}") into decompressedFiles

  script:

  """
  set -euo pipefail
  realReference=`readlink $reference`
  gzip -d -c \$realReference > $reference.baseName
  """
}

fastaFile = Channel.create()
vcfFiles = Channel.create()
decompressedFiles
  .choice(fastaFile, vcfFiles) {it =~ ".fasta" ? 0 : 1}

fastaForBWA = Channel.create()
fastaForPicard = Channel.create()
fastaForSAMTools = Channel.create()

fastaFile.into(fastaForBWA,fastaForPicard,fastaForSAMTools)

process BuildBWAindexes {
  tag {reference}

  publishDir "References/" + params.genome, mode: 'copy'

  input:
    file(reference) from fastaForBWA

  output:
    set file(reference), file("*.{amb,ann,bwt,pac,sa}") into bwaIndexes

  script:

  """
  set -euo pipefail
  bwa index $reference
  """
}

process BuildPicardIndex {
  tag {reference}

  publishDir "References/" + params.genome, mode: 'copy'

  input:
    file(reference) from fastaForPicard

  output:
    file("*.dict") into picardIndex

  script:
  """
  set -euo pipefail
  java -Xmx${task.memory.toGiga()}g \
  -jar \$PICARD_HOME/picard.jar \
  CreateSequenceDictionary \
  REFERENCE=$reference \
  OUTPUT=${reference.baseName}.dict
  """
}

process BuildSAMToolsIndex {
  tag {reference}

  publishDir "References/" + params.genome, mode: 'copy'

  input:
    file(reference) from fastaForSAMTools

  output:
    file("*.fai") into samtoolsIndex

  script:
  """
  set -euo pipefail
  samtools faidx $reference
  """
}

process BuildVCFIndex {
  tag {reference}

  publishDir "References/" + params.genome, mode: 'copy'

  input:
    file(reference) from vcfFiles

  output:
    set file(reference), file("*.idx") into vcfIndexed

  script:
  """
  set -euo pipefail
  \$IGVTOOLS_HOME/igvtools index $reference
  """
}
