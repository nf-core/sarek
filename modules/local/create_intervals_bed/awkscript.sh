awk -vFS="\t" '{
            t = $5  # runtime estimate
            if (t == "") {
                # no runtime estimate in this row, assume default value
                t = ($3 - $2) / 20
            }
            if (name == "" || (chunk > 600 && (chunk + t) > longest * 1.05)) {
                # start a new chunk
                name = sprintf("%s_%d-%d.bed", $1, $2+1, longest)
                chunk = 0
                longest = 0
            }
            if (t > longest)
                longest = t
            chunk += t
            print \$0 > name
        }'
