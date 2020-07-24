// process BASE_RECALIBRATION {
//     label 'cpus_1'

//     tag "${idPatient}-${idSample}-${intervalBed.baseName}"

//     input:
//         tuple idPatient, idSample, file(bam), file(bai), file(intervalBed) //from bamBaseRecalibrator
//         path dbsnp //from dbsnp
//         path dbsnpIndex// from dbsnp_tbi
//         path fasta //from fasta
//         path dict // from dict
//         path fastaFai // from fai
//         path knownIndels // from known_indels
//         path knownIndelsIndex // from known_indels_tbi

//     output:
//         tuple idPatient, idSample, file "${prefix}${idSample}.recal.table", emit: tableGatherBQSRReports
//         tuple idPatient, idSample, emit: recalTableTSVnoInt

//     //when: params.known_indels

//     script:
//     dbsnpOptions = params.dbsnp ? "--known-sites ${dbsnp}" : ""
//     knownOptions = params.known_indels ? knownIndels.collect{"--known-sites ${it}"}.join(' ') : ""
//     prefix = params.no_intervals ? "" : "${intervalBed.baseName}_"
//     intervalsOptions = params.no_intervals ? "" : "-L ${intervalBed}"
//     // TODO: --use-original-qualities ???
//     """
//     gatk --java-options -Xmx${task.memory.toGiga()}g \
//         BaseRecalibrator \
//         -I ${bam} \
//         -O ${prefix}${idSample}.recal.table \
//         --tmp-dir . \
//         -R ${fasta} \
//         ${intervalsOptions} \
//         ${dbsnpOptions} \
//         ${knownOptions} \
//         --verbosity INFO
//     """
// }