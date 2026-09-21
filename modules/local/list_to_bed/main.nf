process LIST_TO_BED {
    tag "${meta.id}"
    label 'process_single'

    input:
    tuple val(meta), path(intervals)

    output:
    tuple val(meta), path("${prefix}.bed"), emit: bed

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"

    // If intervals file is in list format,
    // Match from end of line (int or ',', int or ',', ':')
    // .list is 1-based inclusive, .bed is 0-based exclusive
    """
    awk 'BEGIN{OFS="\\t"}
      /^#/ || NF==0 {next}
      {
        i = match(\$0, /:[0-9,]+(-[0-9,]+)?\$/)
        if (i == 0) { print "No coordinates: " \$0 > "/dev/stderr"; exit 1 }
        chr = substr(\$0, 1, i-1)
        rng = substr(\$0, i+1); gsub(/,/, "", rng)
        n = split(rng, a, "-")
        print chr, a[1]-1, (n == 2 ? a[2] : a[1])
      }' ${intervals} > ${prefix}.bed
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.bed
    """
}
