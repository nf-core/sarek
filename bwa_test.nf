//params.reads = "data/tiny/tiny_normal_L001_R{1,2}.fastq.gz"
params.reads = "data/tiny/tiny_normal_L00?_R{1,2}.fastq.gz"
params.genome = "/sw/data/uppnex/ToolBox/ReferenceAssemblies/hg38make/bundle/2.8/b37/human_g1k_v37_decoy.fasta"
genome_file = file(params.genome)
Channel
    .fromFilePairs( params.reads )                                             
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }  
    .set { read_pairs }

process align {
    module "bwa/0.7.13:samtools/1.3"

    input:
    set pair_id, file(reads) from read_pairs
  
    output:
    file "*.bam" into bam
 
    script: 
    """
    bwa mem -R "@RG\tID:normal\tSM:normal\tLB:normal\tPL:illumina" -B 3 -t 1 -M ${params.genome} ${reads} |samtools view -bS -t /sw/data/uppnex/ToolBox/ReferenceAssemblies/hg38make/bundle/2.8/b37/human_g1k_v37_decoy.fasta - |   samtools sort - > test.bam
    """
}
