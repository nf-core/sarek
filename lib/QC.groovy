class QC {
// Run bamQC on vcf file
  static def bamQC(bam, idSample, mem) {
    """
    qualimap --java-mem-size=${mem.toGiga()}G \
    bamqc \
    -bam ${bam} \
    -outdir ${idSample} \
    -outformat HTML
    """
  }

// Run bcftools on vcf file
  static def bcftools(vcf) {
    """
    bcftools stats ${vcf} > ${vcf.baseName}.bcf.tools.stats.out
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
    --out ${vcf.baseName}

    vcftools \
    --gzvcf ${vcf} \
    --TsTv-by-count \
    --out ${vcf.baseName}

    vcftools \
    --gzvcf ${vcf} \
    --TsTv-by-qual \
    --out ${vcf.baseName}

    vcftools \
    --gzvcf ${vcf} \
    --FILTER-summary \
    --out ${vcf.baseName}
    """
  }

// Get BCFtools version
  static def getVersionBCFtools() {
    """
    bcftools version > v_bcftools.txt
    """
  }

// Get GATK version
  static def getVersionGATK() {
    """
    echo "GATK version"\$(java -jar \$GATK_HOME/GenomeAnalysisTK.jar --version 2>&1) > v_gatk.txt
    """
  }

// Get Manta version
  static def getVersionManta() {
    """
    cat \$MANTA_INSTALL_PATH/lib/python/configBuildTimeInfo.py | grep workflowVersion > v_manta.txt
    """
  }

// Get SnpEFF version
  static def getVersionSnpEFF() {
    """
    echo "SNPEFF version"\$(java -jar \$SNPEFF_HOME/snpEff.jar -h 2>&1) > v_snpeff.txt
    """
  }

// Get Strelka version
  static def getVersionStrelka() {
    """
    cat \$STRELKA_INSTALL_PATH/lib/python/configBuildTimeInfo.py | grep workflowVersion > v_strelka.txt
    """
  }

// Get VCFtools version
  static def getVersionVCFtools() {
    """
    vcftools --version > v_vcftools.txt
    """
  }

// Get VEP version
  static def getVersionVEP() {
    """
    vep --help > v_vep.txt
    """
  }
}
