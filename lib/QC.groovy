class QC {
// Run bcftools on vcf file
  static def bcftools(vcf) {
    """
    bcftools stats ${vcf} > ${vcf.minus(".ann").minus(".vcf").minus(".gz")[0]}.bcf.tools.stats.out
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
    --out ${vcf.minus(".ann").minus(".vcf").minus(".gz")[0]}

    vcftools \
    --gzvcf ${vcf} \
    --TsTv-by-count \
    --out ${vcf.minus(".ann").minus(".vcf").minus(".gz")[0]}

    vcftools \
    --gzvcf ${vcf} \
    --TsTv-by-qual \
    --out ${vcf.minus(".ann").minus(".vcf").minus(".gz")[0]}

    vcftools \
    --gzvcf ${vcf} \
    --FILTER-summary \
    --out ${vcf.minus(".ann").minus(".vcf").minus(".gz")[0]}
    """
  }
}
