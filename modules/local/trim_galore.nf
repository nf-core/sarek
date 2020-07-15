process TRIM_GALORE {
    label 'TrimGalore'

    tag "${idPatient}-${idRun}"

    publishDir "${params.outdir}/Reports/${idSample}/TrimGalore/${idSample}_${idRun}", mode: params.publish_dir_mode,
      saveAs: {filename ->
        if (filename.indexOf("_fastqc") > 0) "FastQC/$filename"
        else if (filename.indexOf("trimming_report.txt") > 0) "logs/$filename"
        else if (params.save_trimmed) filename
        else null
      }

    input:
        tuple val(idPatient), val(idSample), val(idRun), file("${idSample}_${idRun}_R1.fastq.gz"), file("${idSample}_${idRun}_R2.fastq.gz")

    output:
        path "*.{html,zip,txt}",  emit: report
        tuple idPatient, idSample, idRun, file("${idSample}_${idRun}_R1_val_1.fq.gz"), file("${idSample}_${idRun}_R2_val_2.fq.gz") , emit: trimmed_reads


    script:
    // Calculate number of --cores for TrimGalore based on value of task.cpus
    // See: https://github.com/FelixKrueger/TrimGalore/blob/master/Changelog.md#version-060-release-on-1-mar-2019
    // See: https://github.com/nf-core/atacseq/pull/65
    def cores = 1
    if (task.cpus) {
      cores = (task.cpus as int) - 4
      if (cores < 1) cores = 1
      if (cores > 4) cores = 4
      }
    c_r1 = params.clip_r1 > 0 ? "--clip_r1 ${params.clip_r1}" : ''
    c_r2 = params.clip_r2 > 0 ? "--clip_r2 ${params.clip_r2}" : ''
    tpc_r1 = params.three_prime_clip_r1 > 0 ? "--three_prime_clip_r1 ${params.three_prime_clip_r1}" : ''
    tpc_r2 = params.three_prime_clip_r2 > 0 ? "--three_prime_clip_r2 ${params.three_prime_clip_r2}" : ''
    nextseq = params.trim_nextseq > 0 ? "--nextseq ${params.trim_nextseq}" : ''
    """
    trim_galore \
        --cores ${cores} \
        --paired \
        --fastqc \
        --gzip \
        ${c_r1} ${c_r2} \
        ${tpc_r1} ${tpc_r2} \
        ${nextseq} \
        ${idSample}_${idRun}_R1.fastq.gz ${idSample}_${idRun}_R2.fastq.gz

    mv *val_1_fastqc.html "${idSample}_${idRun}_R1.trimmed_fastqc.html"
    mv *val_2_fastqc.html "${idSample}_${idRun}_R2.trimmed_fastqc.html"
    mv *val_1_fastqc.zip "${idSample}_${idRun}_R1.trimmed_fastqc.zip"
    mv *val_2_fastqc.zip "${idSample}_${idRun}_R2.trimmed_fastqc.zip"
    """
}
