#!/usr/bin/env nextflow
 
/*
 * Defines pipeline parameters in order to specify the refence genomes
 * and read pairs by using the command line options
 */

params.genome = "/sw/data/uppnex/ToolBox/ReferenceAssemblies/hg38make/bundle/2.8/b37/human_g1k_v37_decoy.fasta"
params.genomeidx = "${params.genome}.fai"
params.genomedict = "/sw/data/uppnex/ToolBox/ReferenceAssemblies/hg38make/bundle/2.8/b37/human_g1k_v37_decoy.dict"
params.out = "$PWD"
params.kgindels = "/sw/data/uppnex/ToolBox/ReferenceAssemblies/hg38make/bundle/2.8/b37/1000G_phase1.indels.b37.vcf"
params.kgidx = "/sw/data/uppnex/ToolBox/ReferenceAssemblies/hg38make/bundle/2.8/b37/1000G_phase1.indels.b37.vcf.idx"
params.dbsnp = "/sw/data/uppnex/ToolBox/ReferenceAssemblies/hg38make/bundle/2.8/b37/dbsnp_138.b37.vcf"
params.dbsnpidx = "/sw/data/uppnex/ToolBox/ReferenceAssemblies/hg38make/bundle/2.8/b37/dbsnp_138.b37.vcf.idx"
params.millsindels = "/sw/data/uppnex/ToolBox/ReferenceAssemblies/hg38make/bundle/2.8/b37/Mills_and_1000G_gold_standard.indels.b37.vcf"
params.millsidx = "/sw/data/uppnex/ToolBox/ReferenceAssemblies/hg38make/bundle/2.8/b37/Mills_and_1000G_gold_standard.indels.b37.vcf.idx"
params.sample = "sample.config.tsv" // override on cl

if (!params.sample) {
  exit 1, "Please specify the sample config file"
}


/*
 * validate input and params
 */


genome_file = file(params.genome)
genome_index = file(params.genomeidx)
genome_dict = file(params.genomedict)
kgindels = file(params.kgindels)
kgidx = file(params.kgidx)
dbsnp = file(params.dbsnp)
dbsnpidx = file(params.dbsnpidx)
millsindels = file(params.millsindels)
millsidx = file(params.millsidx)
sconfig = file(params.sample)


if( !genome_file.exists() ) exit 1, "Missing reference: ${genome_file}"
if( !genome_dict.exists() ) exit 1, "Missing index: ${genome_dict}"
if( !genome_index.exists() ) exit 1, "Missing index: ${genome_index}"
if( !kgindels.exists() ) exit 1, "Missing vcf: ${kgindels}"
if( !dbsnp.exists() ) exit 1, "Missing vcf: ${dbsnp}"
if( !millsindels.exists() ) exit 1, "Missing vcf: ${millsindels}"


/*
 * Read config file, lets presume its "subject sample fastq1 fastq2"
 * for now and channel this out for mapping
 *
 */


fastqs = Channel
.from(sconfig.readLines())
.map { line ->
  list = line.split()
  mergeId = list[0]
  id = list[1]
  idRun = list[2]
  fq1path = file(list[3])
  fq2path = file(list[4])
  [ mergeId, id, idRun, fq1path, fq2path ]
}


/*
 * processes
 *
 */	


process MappingBwa {

	module 'bioinfo-tools'
	module 'bwa'
	module 'samtools/1.3'

	cpus 1

	input:
	file genome_file
	set mergeId, id, idRun, file(fq1), file(fq2) from fastqs

	output:
//	file "${name}.bam" into mapped_bam
	set mergeId, id, idRun, file("${idRun}.bam") into bams

//	script:
//	lanePattern = ~/_L00[1-8]/
//	id = (name - lanePattern)



// here I use params.genome for bwa ref so I dont have to link to all bwa index files


	script:
	rgString="\"@RG\\tID:${idRun}\\tSM:${id}\\tLB:${id}\\tPL:illumina\""

	"""
	bwa mem -R ${rgString} -B 3 -t ${task.cpus} -M ${params.genome} ${fq1} ${fq2} | samtools view -bS -t ${genome_index} - | samtools sort - > ${idRun}.bam
	"""

}


/*
 * Borrowed code from chip.nf
 * 
 * Now, we decide whether bam is standalone or should be merged by sample (id (column 1) from channel bams)
 *  
 */

// Merge or rename bam
singleBam = Channel.create()
groupedBam = Channel.create()

//bams.groupTuple(by: [0,3,4])
bams.groupTuple(by: [1])
.choice(singleBam, groupedBam) {
  it[3].size() > 1 ? 1 : 0
}

process MergeBam {
    
    module 'bioinfo-tools'
    module 'samtools/1.3'

    input:
    set mergeId, id, idRun, file(bam) from groupedBam

    output:
    set mergeId, id, idRun, file("${id}.bam") into mergedBam
//    set mergeId, id, file("${id}.bam") into mergedBam

    script:
//    cpus = task.cpus
    idRun = idRun.sort().join(':')

    """
    echo ${mergeId} ${id} ${idRun} ${bam} > ble
    samtools merge ${id}.bam ${bam}
    """

//    echo ${idRun} > id.txt
//    (
//      samtools view -H ${bam} | egrep -v '@RG|@PG';
//      for f in ${bam}; do 
//        samtools view -H \$f | grep '@RG';
//      done
//    ) > header.txt && \
//    samtools merge -h header.txt ${id}.bam ${bam}

}

/*
 * merge all merged and single bams to a single channel
 */


bamList = Channel.create()
bamList = singleBam
.mix(mergedBam)
.map { mergeId, id, idRun, bam -> [mergeId[0], id, bam].flatten()
//.map { mergeId, id, idRun, bam -> [mergeId, id, bam].flatten()
}


/*
 *  mark duplicates all bams
 */


process MarkDuplicates {

//	module 'bioinfo-tools'
//	module 'picard'


	input:
//	set mergeId, id, idRun, file(mBam) from bamList
	set mergeId, id, file(mBam) from bamList
//
//	wtf is coming through the channel??
//

	
	output:
//	set mergeId, idRun, file("${idRun}.md.bam"), file("${idRun}.md.bai") into markdupBamInts
	set mergeId, idRun, file("${id}.md.bam"), file("${id}.md.bai") into markdupBamInts
//	set mergeId, idRun, file("${idRun}.md.bam"), file("${idRun}.md.bai") into markdupBam
	set mergeId, idRun, file("${id}.md.bam"), file("${id}.md.bai") into markdupBam

//	echo "${mergeId} : ${id} : ${idRun} : ${mBam}" > ble


	"""
	echo "${mergeId} : ${id} : ${mBam}" > ble
	java -Xmx7g -jar /sw/apps/bioinfo/picard/1.118/milou/MarkDuplicates.jar \
	INPUT=${mBam} \
	METRICS_FILE=${mBam}.metrics \
	TMP_DIR=. \
	ASSUME_SORTED=true \
	VALIDATION_STRINGENCY=LENIENT \
	CREATE_INDEX=TRUE \
	OUTPUT=${id}.md.bam
	"""

// 	debug
//	echo ${mergeId} ${id} ${idRun} ${mergedBam} > ble

}

/*
 * create realign intervals, use both tumor+normal as input
 */

// group the bams by overall subject/patient id (mergeId)
mdbi=Channel.create()
mdbi=markdupBamInts
.groupTuple()
//. map {$mergeId, $idRun, $bam, $bai -> $mergeId, $idRun "-I $bam", $bai}

// group the bams by overall subject/patient id (mergeId)
mdb=Channel.create()
mdb=markdupBam
.groupTuple()
//. map {$mergeId, $idRun, $bam, $bai -> $mergeId, $idRun "-I $bam", $bai}


process CreateIntervals {

	cpus 4
	
	input:
	set mergeId, id, file(mdBam), file(mdBai) from mdbi
	file gf from genome_file 
	file gi from genome_index
	file gd from genome_dict
	file ki from kgindels
	file mi from millsindels

	output:
	file("${mergeId}.intervals") into intervals


	script:
//	ifile="${id}.md.bam"
//	input="-I $mdBam.join(' -I ')"
//	input = GATKInput.spread(mdBam)
	input = mdBam.collect{"-I $it"}.join(' ')


	"""
	echo ${mergeId} ${id} ${mdBam} > ble
	java -Xmx7g -jar /sw/apps/bioinfo/GATK/3.3.0/GenomeAnalysisTK.jar \
	-T RealignerTargetCreator \
	$input \
	-R $gf \
	-known $ki \
	-known $mi \
	-nt ${task.cpus} \
	-o ${mergeId}.intervals

 	"""	

}


/*
 * realign, use nWayOut to split into tumor/normal again
 */


process realign {


	input:
	set mergeId, id, file(mdBam), file(mdBai) from mdb
	file gf from genome_file
	file gi from genome_index
	file gd from genome_dict
	file ki from kgindels
	file mi from millsindels
	file intervals from intervals

	output:
	set mergeId, id, file("${id}.real.bam") into realBams


	script:
	input = mdBam.collect{"-I $it"}.join(' ')

	"""
	java -Xmx7g -jar /sw/apps/bioinfo/GATK/3.3.0/GenomeAnalysisTK.jar \
	-T IndelRealigner \
	$input \
	-R $gf \
	-targetIntervals $intervals \
	-known $ki \
	-known $mi \
	-nWayOut '.real.bam'
	"""

}



/*
process create_recal_table_tumor {

	cpus 2

	input:
	file tumor_real_bam_table
	file tumor_real_bai_table
	file genome_file
	file genome_dict
	file genome_index
	file dbsnp
	file dbsnpidx
	file kgindels
	file kgidx
	file millsindels
	file millsidx

	output:
	file '*.tumor.recal.table' into tumor_recal_table


	"""
	java -Xmx7g \
	-jar /sw/apps/bioinfo/GATK/3.3.0/GenomeAnalysisTK.jar \
	-T BaseRecalibrator \
	-l INFO -R $genome_file \
	-I $tumor_real_bam_table \
	-knownSites $dbsnp \
	-knownSites $kgindels \
	-knownSites $millsindels \
	-nct ${task.cpus} \
	-o ${params.sample}.tumor.recal.table
	"""
}
*/

/*
process recalibrate_bam_tumor {


	
	input:
	file tumor_real_bam_recal
	file tumor_real_bai_recal
	file genome_file
	file genome_dict
	file genome_index
	file dbsnp
	file dbsnpidx
	file kgindels
	file kgidx
	file millsindels
	file millsidx
	file tumor_recal_table

	output:
	file '*.tumor.recal.bam' into tumor_recal_bam
	file '*.tumor.recal.bai' into tumor_recal_bai

	"""
	java -Xmx7g -Djava.io.tmpdir=\$SNIC_TMP \
	-jar /sw/apps/bioinfo/GATK/3.3.0/GenomeAnalysisTK.jar \
	-R $genome_file \
	-I $tumor_real_bam_recal \
	-T PrintReads \
	--BQSR $tumor_recal_table \
	-o ${params.sample}.tumor.recal.bam
	"""


}
*/

/*
process genotype_gvcf_tumor {

	cpus 2
	
	input:
	file tumor_recal_bam
	file tumor_recal_bai
	file genome_file
	file genome_dict
	file genome_index
	file dbsnp
	file dbsnpidx
	file kgindels
	file kgidx
	file millsindels
	file millsidx	

	output:
	file '*.tumor.g.vcf.gz' into 'tumor_gvcf'

	"""
        java -Xmx7g \
	-jar /sw/apps/bioinfo/GATK/3.3.0/GenomeAnalysisTK.jar \
	-R $genome_file \
	-T HaplotypeCaller \
	-I $tumor_recal_bam \
	--emitRefConfidence GVCF \
	--variant_index_type LINEAR \
	--dbsnp $dbsnp \
	--variant_index_parameter 128000 \
	-nct ${task.cpus} \
	-o ${params.sample}.tumor.g.vcf.gz
	"""

}
*/

/*
process create_recal_table_normal {

	cpus 2

	input:
	file normal_real_bam_table
	file normal_real_bai_table
	file genome_file
	file genome_dict
	file genome_index
	file dbsnp
	file dbsnpidx
	file kgindels
	file kgidx
	file millsindels
	file millsidx

	output:
	file '*.normal.recal.table' into normal_recal_table


	"""
	java -Xmx7g \
	-jar /sw/apps/bioinfo/GATK/3.3.0/GenomeAnalysisTK.jar \
	-T BaseRecalibrator \
	-l INFO -R $genome_file \
	-I $normal_real_bam_table \
	-knownSites $dbsnp \
	-knownSites $kgindels \
	-knownSites $millsindels \
	-nct ${task.cpus} \
	-o ${params.sample}.normal.recal.table
	"""
}
*/

/*
process recalibrate_bam_normal {


	
	input:
	file normal_real_bam_recal
	file normal_real_bai_recal
	file genome_file
	file genome_dict
	file genome_index
	file dbsnp
	file dbsnpidx
	file kgindels
	file kgidx
	file millsindels
	file millsidx
	file normal_recal_table

	output:
	file '*.normal.recal.bam' into normal_recal_bam
	file '*.normal.recal.bai' into normal_recal_bai

	"""
	java -Xmx7g -Djava.io.tmpdir=\$SNIC_TMP \
	-jar /sw/apps/bioinfo/GATK/3.3.0/GenomeAnalysisTK.jar \
	-R $genome_file \
	-I $normal_real_bam_recal \
	-T PrintReads \
	--BQSR $normal_recal_table \
	-o ${params.sample}.normal.recal.bam
	"""


}
*/


/*
process genotype_gvcf_normal {

	cpus 2
	
	input:
	file normal_recal_bam
	file normal_recal_bai
	file genome_file
	file genome_dict
	file genome_index
	file dbsnp
	file dbsnpidx
	file kgindels
	file kgidx
	file millsindels
	file millsidx	

	output:
	file '*.normal.g.vcf.gz' into 'normal_gvcf'

	"""
        java -Xmx7g \
	-jar /sw/apps/bioinfo/GATK/3.3.0/GenomeAnalysisTK.jar \
	-R $genome_file \
	-T HaplotypeCaller \
	-I $normal_recal_bam \
	--emitRefConfidence GVCF \
	--variant_index_type LINEAR \
	--dbsnp $dbsnp \
	--variant_index_parameter 128000 \
	-nct ${task.cpus} \
	-o ${params.sample}.normal.g.vcf.gz
	"""

}

*/


// ############################### FUNCTIONS

 def readPrefix( Path actual, template ) {

    final fileName = actual.getFileName().toString()

    def filePattern = template.toString()
    int p = filePattern.lastIndexOf('/')
    if( p != -1 ) filePattern = filePattern.substring(p+1)
    if( !filePattern.contains('*') && !filePattern.contains('?') ) 
        filePattern = '*' + filePattern 
  
    def regex = filePattern
                    .replace('.','\\.')
                    .replace('*','(.*)')
                    .replace('?','(.?)')
                    .replace('{','(?:')
                    .replace('}',')')
                    .replace(',','|')

    def matcher = (fileName =~ /$regex/)
    if( matcher.matches() ) {  
        def end = matcher.end(matcher.groupCount() )      
        def prefix = fileName.substring(0,end)
        while(prefix.endsWith('-') || prefix.endsWith('_') || prefix.endsWith('.') ) 
          prefix=prefix[0..-2]
          
        return prefix
    }
    
    return fileName
}


/* unused crap: ######################################################################


process mark_duplicates_normal {

	module 'bioinfo-tools'
	module 'picard'


	input:
	file normal_bam
	
	output:
	file '*.normal.md.bam' into normal_md_bam_intervals, normal_md_bam_real
	file '*.normal.md.bai' into normal_md_bai_intervals, normal_md_bai_real

	"""
	java -Xmx7g -jar /sw/apps/bioinfo/picard/1.118/milou/MarkDuplicates.jar \
	INPUT=${normal_bam} \
	METRICS_FILE=${normal_bam}.metrics \
	TMP_DIR=. \
	ASSUME_SORTED=true \
	VALIDATION_STRINGENCY=LENIENT \
	CREATE_INDEX=TRUE \
	OUTPUT=${params.sample}.normal.md.bam	
	"""

}


process merge_bam {

	module 'bioinfo-tools'
	module 'picard'

	input:
	file (bams:'*') from mapped_bam.toList() 

	output:
	file ble


	
	
//	java -jar picard.jar MergeSamFiles I=input_1.bam I=input_2.bam O=merged_files.bam

	"""
	ls ${bams} > ble
        """
}


 

process mapping_tumor_bwa {

	module 'bioinfo-tools'
	module 'bwa'
	module 'samtools/1.3'


	cpus 1

	input:
	file genome_file
	
	file tp1
	file tp2

	output:
	file '*.tumor.bam' into tumor_bam

	"""
	bwa mem -R "@RG\\tID:${params.sample}.tumor\\tSM:${params.sample}\\tLB:${params.sample}.tumor\\tPL:illumina" \
	-B 3 -t ${task.cpus} \
	-M ${params.genome} ${tp1} ${tp2} \
	| samtools view -bS -t ${genome_index} - \
	| samtools sort - > ${params.sample}.tumor.bam
	"""	

}

process mapping_normal_bwa {

	module 'bioinfo-tools'
	module 'bwa'
	module 'samtools/1.3'

	cpus 1


	input:
	file genome_file
	file np1
	file np2

	output:
	file '*.normal.bam' into normal_bam 

	"""
	bwa mem -R "@RG\\tID:${params.sample}.normal\\tSM:${params.sample}\\tLB:${params.sample}.normal\\tPL:illumina" \
	-B 3 -t ${task.cpus} \
	-M ${params.genome} ${np1} ${np2} \
	| samtools view -bS -t ${genome_index} - \
	| samtools sort - > ${params.sample}.normal.bam
	"""	

}




if (!params.index) {
  exit 1, "Please specify the input table file"

//params.sample = "tcga.cl" // override with --sample <SAMPLE>
//params.tpair1 = "data/${params.sample}.tumor_R1.fastq.gz"
//params.npair1 = "data/${params.sample}.normal_R1.fastq.gz"
//params.tpair2 = "data/${params.sample}.tumor_R2.fastq.gz"
//params.npair2 = "data/${params.sample}.normal_R2.fastq.gz"


//tp1 = file(params.tpair1)
//tp2 = file(params.tpair2)
//np1 = file(params.npair1)
//np2 = file(params.npair2)


//if( !tp1.exists() ) exit 2, "Missing read ${tp1}"
//if( !tp2.exists() ) exit 2, "Missing read ${tp2}"
//if( !np1.exists() ) exit 2, "Missing read ${np1}"
//if( !np2.exists() ) exit 2, "Missing read ${np2}"


// Pattern to grab all fastq pairs:
params.r = "data/*_R[12].fastq.gz"
Channel
    .fromPath( params.r )
    .ifEmpty { error "Cannot find any reads matching: ${params.r}" }
    .map { path -> 
       def prefix = readPrefix(path, params.r)
       tuple(prefix, path) 
    }
    .groupTuple(sort: true)
    .set { read } 


// read .subscribe {println it}

//      ln -s ${name}.bam ../../../${id}
//	mkdir -f bam_${name}
//	mkdir -p ../../../${id}
//       	bwa mem -R "@RG\\tID:${params.sample}\\tSM:${params.sample}\\tLB:${params.sample}\\tPL:illumina" -B 3 -t ${task.cpus} -M ${params.genome} ${reads} | samtools view -bS -t ${genome_index} - | samtools sort - > ${name}.bam


*/