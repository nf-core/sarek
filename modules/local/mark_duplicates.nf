process MARK_DUPLICATES {
     label 'cpus_16'
     tag "${idPatient}-${idSample}"
     publishDir params.outdir, mode: params.publish_dir_mode,
         saveAs: {
             if (it == "${idSample}.bam.metrics") "Reports/${idSample}/MarkDuplicates/${it}"
             else "Preprocessing/${idSample}/DuplicatesMarked/${it}"
         }
     input:
         tuple idPatient, idSample, path("${idSample}.bam") 
     output:
         tuple idPatient, idSample, path("${idSample}.md.bam"), path("${idSample}.md.bam.bai"), emit:  bam_duplicates_marked
         tuple idPatient, idSample, emit: tsv_bam_duplicates_marked
         path "${idSample}.bam.metrics", emit: duplicates_marked_report //is optional , applies when skip_qc is used(not implemented yet)
     
     //when: !(params.skip_markduplicates)
     
     script:
     markdup_java_options = task.memory.toGiga() > 8 ? params.markdup_java_options : "\"-Xms" +  (task.memory.toGiga() / 2).trunc() + "g -Xmx" + (task.memory.toGiga() - 1) + "g\""
     //metrics = 'markduplicates' in skip_qc ? '' : "-M ${idSample}.bam.metrics"
     metrics = "-M ${idSample}.bam.metrics"
     if (params.no_gatk_spark)
     """
     gatk --java-options ${markdup_java_options} \
         MarkDuplicates \
         --MAX_RECORDS_IN_RAM 50000 \
         --INPUT ${idSample}.bam \
         --METRICS_FILE ${idSample}.bam.metrics \
         --TMP_DIR . \
         --ASSUME_SORT_ORDER coordinate \
         --CREATE_INDEX true \
         --OUTPUT ${idSample}.md.bam
     mv ${idSample}.md.bai ${idSample}.md.bam.bai
     """
     else
     """
     gatk --java-options ${markdup_java_options} \
         MarkDuplicatesSpark \
         -I ${idSample}.bam \
         -O ${idSample}.md.bam \
         ${metrics} \
         --tmp-dir . \
         --create-output-bam-index true \
         --spark-master local[${task.cpus}]
     """
 }