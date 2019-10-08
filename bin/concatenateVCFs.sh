#!/usr/bin/env bash
set -euo pipefail

# This script concatenates all VCF files that are in the local directory,
# that were created from different intervals to make a single final VCF

usage() { echo "Usage: $0 [-i genome_index_file] [-o output.file.no.gz.extension] <-t target.bed> <-c cpus>" 1>&2; exit 1; }

while [[ $# -gt 0 ]]
do
  key=$1
  case $key in
		-i)
			genomeIndex=$2
			shift # past argument
	    shift # past value
			;;
		-c)
			cpus=$2
			shift # past argument
	    shift # past value
			;;
		-o)
			outputFile=$2
			shift # past argument
	    shift # past value
			;;
		-t)
			targetBED=$2
			shift # past argument
	    shift # past value
			;;
		*)
			usage
			shift # past argument
			;;
	esac
done

if [ -z ${genomeIndex} ]; then echo "Missing index file "; usage; fi
if [ -z ${cpus} ]; then echo "No CPUs defined: setting to 1"; cpus=1; fi
if [ -z ${outputFile} ]; then echo "Missing output file name"; usage; fi

# First make a header from one of the VCF
# Remove interval information from the GATK command-line, but leave the rest
FIRSTVCF=$(ls *.vcf | head -n 1)
sed -n '/^[^#]/q;p' $FIRSTVCF | \
awk '!/GATKCommandLine/{print}/GATKCommandLine/{for(i=1;i<=NF;i++){if($i!~/intervals=/ && $i !~ /out=/){printf("%s ",$i)}}printf("\n")}' \
> header

# Get list of contigs from the FASTA index (.fai)
# ##contig header in the VCF cannot be used as it is optional (FreeBayes does not save it, for example)

CONTIGS=($(cut -f1 ${genomeIndex}))

# Concatenate VCFs in the correct order
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

set +u

# Now we have the concatenated VCF file, check for WES/panel targets, and generate a subset if there is a BED provided
if [ ! -z ${targetBED+x} ]; then
	echo "Target is $targetBED - Selecting subset..."
	bcftools isec --targets-file ${targetBED} rawcalls.vcf.gz | bgzip -@${cpus} > ${outputFile}.gz
	tabix ${outputFile}.gz
else
	# Rename the raw calls as WGS results
	for f in rawcalls*; do mv -v $f ${outputFile}${f#rawcalls.vcf}; done
fi
