# Genomes

CAW currently uses GRCh37 by default. Support for GRCh38 is not fully working!


## GRCh37

...


## GRCh38

Use `--genome=GRCh38` to map against GRCh38. Before doing so and if you
are not on Uppmax, you need to adjust the settings in `genomes.config` to your
needs.

To get the needed files, download the GATK bundle for GRCh38 from
<ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/>.

The MD5SUM of `Homo_sapiens_assembly38.fasta` included in that file is
7ff134953dcca8c8997453bbb80b6b5e.

From the `beta/` directory, which seems to be an older version of the bundle,
only Homo_sapiens_assembly38.known_indels.vcf* is needed. Also, you can omit
dbsnp_138* and dbsnp_144 files as we use dbsnp_146. The old ones also use the
wrong chromosome naming convention.

Afterwards, the following needs to be done:

    gunzip Homo_sapiens_assembly38.fasta.gz
    bwa index -6 Homo_sapiens_assembly38.fasta
    awk '!/^@/{printf("%s:%d-%d\n", $1, $2, $3)}' wgs_calling_regions.hg38.interval_list > wgs_calling_regions.hg38.list
