BEGIN { 
    OFS="\t"
}
# if it is a header print it out, but get rid of HPV and oter decoy parts
/^@/ && (!/NC_007605/ && !/hs37d5/){ 
    print
}
# if it is a read, print only if the other part of the read is not in a decoy region 
# I am not entirely sure it is the right thing, but will ask Jesper
!/^@/ &&($7 !~/hs37d5/ && $7 !~/NC_007605/ && $3!~/hs37d5/ && $3 !~/NC_007605/) {
    print
}
