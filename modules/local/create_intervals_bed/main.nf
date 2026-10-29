process CREATE_INTERVALS_BED {
    tag "${intervals}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/gawk:5.3.0'
        : 'biocontainers/gawk:5.3.0'}"

    input:
    path intervals
    tuple val(meta2), path(fasta_fai)
    val nucleotides_per_second

    output:
    path ("*.bed"),        emit: bed
    path ("sorted/*.bed"), emit: sorted_bed, optional: true
    tuple val("${task.process}"), val('gawk'), eval("awk --version | head -n1 | sed 's/GNU Awk //; s/, .*//'"), emit: versions_gawk, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // If intervals file is in BED format,
    // Fifth column is interpreted to contain runtime estimates
    // Which is then used to combine short-running jobs
    // Rows are first sorted in reference order: tabix rejects unsorted positions, and
    // GatherPileupSummaries orders chunks by their first row only (nf-core/sarek#2281)
    if (intervals.toString().toLowerCase().endsWith("bed")) {
        """
        mkdir sorted

        # Sort rows by fasta index order, then by start, keeping header lines on top
        awk -vFS="\t" '
        NR == FNR {
            rank[\$1] = FNR
            next
        }
        /^(#|track|browser)/ {
            print
            next
        }
        {
            # contigs missing from the fasta index go last, in order of appearance
            if (!(\$1 in rank))
                rank[\$1] = length(rank) + 1
            # zero-padded so that string order equals numeric order, FNR keeps the sort stable
            row[sprintf("%012d %012d %012d", rank[\$1], \$2, FNR)] = \$0
        }
        END {
            PROCINFO["sorted_in"] = "@ind_str_asc"
            for (key in row)
                print row[key]
        }' ${fasta_fai} ${intervals} > sorted/${intervals.name}

        awk -vFS="\t" '{
            t = \$5  # runtime estimate
            if (t == "") {
                # no runtime estimate in this row, assume default value
                t = (\$3 - \$2) / ${nucleotides_per_second}
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
        }' sorted/${intervals.name}
        """
    }
    else if (intervals.toString().toLowerCase().endsWith("interval_list")) {
        """
        grep -v '^@' ${intervals} | awk -vFS="\t" '{
            name = sprintf("%s_%d-%d", \$1, \$2, \$3);
            printf("%s\\t%d\\t%d\\n", \$1, \$2-1, \$3) > name ".bed"
        }'
        """
    }
    else {
        """
        awk -vFS="[:-]" '{
            name = sprintf("%s_%d-%d", \$1, \$2, \$3);
            printf("%s\\t%d\\t%d\\n", \$1, \$2-1, \$3) > name ".bed"
        }' ${intervals}
        """
    }

    stub:
    def prefix = task.ext.prefix ?: "${intervals.baseName}"
    def sorted = intervals.toString().toLowerCase().endsWith("bed") ? "mkdir sorted && touch sorted/${intervals.name}" : ""
    """
    touch ${prefix}.stub.bed
    ${sorted}
    """
}
