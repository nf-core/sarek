#  This awk script includes a variant caller flag into the INFO fields
#  you should call it like
# awk -f addCallerINFO.awk -v caller=VCCALLER variants.vcf > new_file.vcf
BEGIN {
    VCSincluded = 0;
    OFS="\t"
}
# metainfo lines
/^##/ {
    if(/^##INFO/) {
        if(0 == VCSincluded) {
            # include the Variant Calling Software info (VCS)
            print "##INFO=<ID="caller",Number=0,Type=Flag,Description=\"Called by "caller">"
            print $0        # print the first INFO line
            VCSincluded = -1
        } else {
            print
        }
    } else {
        print
    }
}
# header line
/^#CHR/{ print }
# remainders are the call lines - we are allowing only PASS-ed entries
!/^#/ && $7 ~/PASS/{
    # we are adding the caller info to the INFO field that is the 8th field
    $8 = caller";"$8
    print
}
