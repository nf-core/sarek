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

// Get BCFtools version
  static def getVersionBCFtools() {
    """
    bcftools version > v_bcftools.txt
    """
  }

// Get GATK version
  static def getVersionGATK() {
    """
		gatk ApplyBQSR --help 2>&1| awk -F/ '/java/{for(i=1;i<=NF;i++){if(\$i~/gatk4/){sub("gatk4-","",\$i);print \$i>"v_gatk.txt"}}}'
    """
  }

// Get Manta version
  static def getVersionManta() {
    """
		configManta.py --version > v_manta.txt
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
		configureStrelkaGermlineWorkflow.py --version > v_strelka.txt
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
