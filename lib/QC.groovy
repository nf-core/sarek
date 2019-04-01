class QC {
// Run bcftools on vcf file
  static def bcftools(vcf) {
    """
    bcftools stats ${vcf} > ${vcf.fileName.toString().minus(".ann").minus(".vcf").minus(".gz")}.bcf.tools.stats.out
    """
  }

// Run samtools stats on bam file
  static def samtoolsStats(bam) {
    """
    samtools stats ${bam} > ${bam}.samtools.stats.out
    """
  }

// Run vcftools on vcf file
  static def vcftools(vcf) {
    """
    vcftools \
    --gzvcf ${vcf} \
    --relatedness2 \
    --out ${vcf.fileName.toString().minus(".ann").minus(".vcf").minus(".gz")}

    vcftools \
    --gzvcf ${vcf} \
    --TsTv-by-count \
    --out ${vcf.fileName.toString().minus(".ann").minus(".vcf").minus(".gz")}

    vcftools \
    --gzvcf ${vcf} \
    --TsTv-by-qual \
    --out ${vcf.fileName.toString().minus(".ann").minus(".vcf").minus(".gz")}

    vcftools \
    --gzvcf ${vcf} \
    --FILTER-summary \
    --out ${vcf.fileName.toString().minus(".ann").minus(".vcf").minus(".gz")}
    """
  }
}
