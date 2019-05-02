#!/usr/bin/env bash
# this script concatenates all VCFs that are in the local directory: the 
# purpose is to make a single VCF from all the VCFs that were created from different intervals

usage() { echo "Usage: $0 [-i genome_index_file] [-o output.file.no.gz.extension] <-t target.bed> <-c cpus>" 1>&2; exit 1; }

while getopts "i:c:o:t:" p; do
	case "${p}" in
		i)
			genomeIndex=${OPTARG}
			;;
		c)
			cpus=${OPTARG}
			;;
		o)
			outputFile=${OPTARG}
			;;
		t)
			targetBED=${OPTARG}
			;;
		*)
			usage
			;;
	esac
done
shift $((OPTIND-1))

if [ -z ${genomeIndex} ]; then echo "Missing index file "; usage; fi
if [ -z ${cpus} ]; then echo "No CPUs defined: setting to 1"; cpus=1; fi
if [ -z ${outputFile} ]; then echo "Missing output file name"; usage; fi

set -euo pipefail

# first make a header from one of the VCF intervals
# get rid of interval information only from the GATK command-line, but leave the rest
FIRSTVCF=$(ls *.vcf | head -n 1)
sed -n '/^[^#]/q;p' $FIRSTVCF | \
awk '!/GATKCommandLine/{print}/GATKCommandLine/{for(i=1;i<=NF;i++){if($i!~/intervals=/ && $i !~ /out=/){printf("%s ",$i)}}printf("\n")}' \
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

set +u

# now we have the concatenated VCF file, check for WES/panel targets, and generate a subset if there is a BED provided
echo "target is $targetBED"
if [ ! -z ${targetBED+x} ]; then
	echo "Selecting subset..."
	bcftools isec --targets-file ${targetBED} rawcalls.vcf.gz | bgzip -@${cpus} > ${outputFile}.gz
	tabix ${outputFile}.gz
else
	# simply rename the raw calls as WGS results
	for f in rawcalls*; do mv -v $f ${outputFile}${f#rawcalls.vcf}; done
fi

