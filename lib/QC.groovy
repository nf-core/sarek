class QC {
// Run bcftools on vcf file
  static def bcftools(vcf) {
    """
    bcftools stats ${vcf} > ${vcf.simpleName}.bcf.tools.stats.out
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
    --out ${vcf.simpleName}

    vcftools \
    --gzvcf ${vcf} \
    --TsTv-by-count \
    --out ${vcf.simpleName}

    vcftools \
    --gzvcf ${vcf} \
    --TsTv-by-qual \
    --out ${vcf.simpleName}

    vcftools \
    --gzvcf ${vcf} \
    --FILTER-summary \
    --out ${vcf.simpleName}
    """
  }


// Get SnpEFF version
  static def getVersionSnpEFF() {
    """
    echo "SNPEFF version"\$(java -jar \$SNPEFF_HOME/snpEff.jar -h 2>&1) > v_snpeff.txt
    """
  }

// Get VEP version
  static def getVersionVEP() {
    """
    /opt/vep/src/ensembl-vep/vep --help > v_vep.txt
    """
  }
}
