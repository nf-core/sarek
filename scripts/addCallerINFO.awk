# header lines
BEGIN {
    VCSincluded = 0;
    OFS="\t"
}
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
/^#CHR/{ print }
# remainders are the call lines 
!/^#/ {
    # we are adding the caller info to the INFO field that is the 8th field
    $8 = caller";"$8
    print
}
