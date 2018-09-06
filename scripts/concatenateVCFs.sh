#!/usr/bin/env bash
# this script concatenates all VCFs that are in the local directory: the 
# purpose is to make a single VCF from all the VCFs that were created from different intervals

set -euo pipefail

genomeIndex=$1
cpus=$2
outputFile=$3
targetBED=$4
# first make a header from one of the VCF intervals
# get rid of interval information only from the GATK command-line, but leave the rest
FIRSTVCF=$(ls *.vcf | head -n 1)
sed -n '/^[^#]/q;p' $FIRSTVCF | \
awk '!/GATKCommandLine/{print}/GATKCommandLine/{for(i=1;i<=NF;i++){if($i!~/intervals=/ && $i !~ /out=/){printf("%s ",$i)}}printf("\\n")}' \
> header

# Get list of contigs from the FASTA index (.fai). We cannot use the ##contig
# header in the VCF as it is optional (FreeBayes does not save it, for example)
CONTIGS=($(cut -f1 ${genomeIndex}))

# concatenate VCFs in the correct order
(
  cat header

  for chr in "${CONTIGS[@]}"; do
    # Skip if globbing would not match any file to avoid errors such as
    # "ls: cannot access chr3_*.vcf: No such file or directory" when chr3
    # was not processed.
    pattern="${chr}_*.vcf"
    if ! compgen -G "${pattern}" > /dev/null; then continue; fi

    # ls -v sorts by numeric value ("version"), which means that chr1_100_
    # is sorted *after* chr1_99_.
    for vcf in $(ls -v ${pattern}); do
      # Determine length of header.
      # The 'q' command makes sed exit when it sees the first non-header
      # line, which avoids reading in the entire file.
      L=$(sed -n '/^[^#]/q;p' ${vcf} | wc -l)

      # Then print all non-header lines. Since tail is very fast (nearly as
      # fast as cat), this is way more efficient than using a single sed,
      # awk or grep command.
      tail -n +$((L+1)) ${vcf}
    done
  done
) | bgzip -@${cpus} > rawcalls.vcf.gz
tabix rawcalls.vcf.gz

# now we have the concatenated VCF file, check for WES/panel targets, and generate a subset if there is a BED provided
if [ -s "${targetBED}" ]; then
	bcftools isec --targets-file ${targetBED} rawcalls.vcf.gz | bgzip -@${cpus} > ${outputFile}.gz
	tabix ${outputFile}.gz
else
	# simply rename the raw calls as WGS results
	for f in rawcalls*; do mv -v $f ${outputFile}${f#rawcalls}; done
fi

