process CREATE_INTERVALS_BED {
    tag "$intervals"

    conda (params.enable_conda ? "anaconda::gawk=5.1.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gawk:5.1.0' :
        'quay.io/biocontainers/gawk:5.1.0' }"

    input:
    path(intervals)

    output:
    path("*.bed")       , emit: bed
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // If intervals file is in BED format,
    // Fifth column is interpreted to contain runtime estimates
    // Which is then used to combine short-running jobs
    if (intervals.toString().toLowerCase().endsWith("bed")) {
        """
        awk -vFS="\t" '{
            t = \$5  # runtime estimate
            if (t == "") {
                # no runtime estimate in this row, assume default value
                t = (\$3 - \$2) / ${params.nucleotides_per_second}
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

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            gawk: \$(awk -Wversion | sed '1!d; s/.*Awk //; s/,.*//')
        END_VERSIONS
        """
    } else if (intervals.toString().toLowerCase().endsWith("interval_list")) {
        """
        grep -v '^@' ${intervals} | awk -vFS="\t" '{
            name = sprintf("%s_%d-%d", \$1, \$2, \$3);
            printf("%s\\t%d\\t%d\\n", \$1, \$2-1, \$3) > name ".bed"
        }'

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            gawk: \$(awk -Wversion | sed '1!d; s/.*Awk //; s/,.*//')
        END_VERSIONS
        """
    } else {
        """
        awk -vFS="[:-]" '{
            name = sprintf("%s_%d-%d", \$1, \$2, \$3);
            printf("%s\\t%d\\t%d\\n", \$1, \$2-1, \$3) > name ".bed"
        }' ${intervals}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            gawk: \$(awk -Wversion | sed '1!d; s/.*Awk //; s/,.*//')
        END_VERSIONS
        """
    }
}
