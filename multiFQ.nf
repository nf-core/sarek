#!/usr/bin/env nextflow
/*
    Cancer Analysis Workflow for tumor/normal samples. Usage:
    
    $ nextflow run simplified.nf -c <config> --sample <sample.tsv>

    Parameters (like location of the genome reference, dbSNP location) have to be given at the command-line.
    We are not giving paramters here, but checking the existence of files.

    You also must give the sample .tsv file that is defining the tumor/normal pairs. See explanation below.
 */

//
// we want to list all the missing parameters at the beginning, hence checking only this value
// and exiting only at the very end if some of the parameters failed
//
allParamsDefined = true    

/*
    Use this closure to loop through all the parameters.
    We can get an AssertionError exception from the file() method as well.
 */
checkExistence = {
    paramFile,paramToCheck -> 
    try {
        paramFile= file(paramToCheck)
        assert paramFile.exists()
    } catch(AssertionError ae) {
        println("Missing file: ${paramFile} ${paramToCheck}")
        allParamsDefined = false;
    }
}

refs = [ "genome_file": params.genome,      // genome reference
    "genome_index" : params.genomeIndex,    // genome reference index
    "genome_dict": params.genomeDict,       // genome reference dictionary
    "kgindels": params.kgIndels,            // 1000 Genomes SNPs 
    "kgidx": params.kgIndex,                // 1000 Genomes SNPs index
    "dbsnp": params.dbsnp,                  // dbSNP
    "dbsnpidx": params.dbsnpIndex,          // dbSNP index
    "millsindels": params.millsIndels,      // Mill's Golden set of SNPs
    "millsidx": params.millsIndex,          // Mill's Golden set index
    "sample": params.sample                 // the sample sheet (multilane data refrence table, see below)
  ]
refs.each(checkExistence)

if(!allParamsDefined) {
    exit 1, "Missing file or parameter: please review your config file. Use the -c <config> command line option and --sample <sample.tsv> for data."
}

/*
    Time to check the sample file. Its format is like: "subject sample lane fastq1 fastq2": 

tcga.cl	tcga.cl.normal	tcga.cl.normal_1	data/tcga.cl.normal_L001_R1.fastq.gz	data/tcga.cl.normal_L001_R2.fastq.gz
tcga.cl	tcga.cl.tumor	tcga.cl.tumor_1	data/tcga.cl.tumor_L001_R1.fastq.gz	data/tcga.cl.tumor_L001_R2.fastq.gz
tcga.cl tcga.cl.tumor	tcga.cl.tumor_2	data/tcga.cl.tumor_L002_R1.fastq.gz	data/tcga.cl.tumor_L002_R2.fastq.gz
 */
sconfig = file(params.sample)
if (!params.sample) {
  exit 1, "Please specify the sample config file"
}

/*
 * Read config file, lets presume its "subject sample fastq1 fastq2"
 * for now and channel this out for mapping
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

fastqs = logChannelContent("FASTQ files and IDs to process: ",fastqs)

/*
 * processes
 */	
process MappingBwa {

	module 'bioinfo-tools'
	module 'bwa'
	module 'samtools/1.3'

	cpus 1

	input:
	file refs["genome_file"]
	set mergeId, id, idRun, file(fq1), file(fq2) from fastqs

	output:
	set mergeId, id, idRun, file("${idRun}.bam") into bams

// here I use params.genome for bwa ref so I dont have to link to all bwa index files


	script:
	rgString="\"@RG\\tID:${idRun}\\tSM:${id}\\tLB:${id}\\tPL:illumina\""

	"""
	bwa mem -R ${rgString} -B 3 -t ${task.cpus} -M ${params.genome} ${fq1} ${fq2} | samtools view -bS -t ${refs["genomeIndex"]} - | samtools sort - > ${idRun}.bam
	"""
}

/*
 * Borrowed code from chip.nf
 * 
 * Now, we decide whether bam is standalone or should be merged by sample (id (column 1) from channel bams)
 * http://www.nextflow.io/docs/latest/operator.html?highlight=grouptuple#grouptuple 
 * 
 */

// Merge or rename bam
singleBam = Channel.create()
groupedBam = Channel.create()

bams.groupTuple(by: [1])
    .choice(singleBam, groupedBam) {
        it[3].size() > 1 ? 1 : 0
    }

singleBam = logChannelContent("Single BAMs before merge:",singleBam )
groupedBam = logChannelContent("Grouped BAMs before merge:",groupedBam)

process MergeBam {
    
    module 'bioinfo-tools'
    module 'samtools/1.3'

    input:
    set mergeId, id, idRun, file(bam) from groupedBam

    output:
    set mergeId, id, idRun, file("${id}.bam") into mergedBam

    script:
    idRun = idRun.sort().join(':')

    """
    echo ${mergeId} ${id} ${idRun} ${bam} > ble
    samtools merge ${id}.bam ${bam}
    """
}

/*
 * merge all merged and single bams to a single channel
 */
bamList = Channel.create()
bamList = singleBam
    .mix(mergedBam)
    .map { mergeId, id, idRun, bam -> [mergeId[0], id, bam].flatten()
}

bamList = logChannelContent("BAM list for MarkDuplicates",bamList)

/*
 *  mark duplicates all bams
 */
process MarkDuplicates {

	module 'bioinfo-tools'
	module 'picard'

	input:
	set mergeId, id, idRun, file(mBam) from bamList
	set mergeId, id, file(mBam) from bamList
//
//	Channel content should be in the log before
//
	output:
	set mergeId, idRun, file("${id}.md.bam"), file("${id}.md.bai") into markdupBamInts
	set mergeId, idRun, file("${id}.md.bam"), file("${id}.md.bai") into markdupBam

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
}

/*
 * create realign intervals, use both tumor+normal as input
 */

// group the bam itervals by overall subject/patient id (mergeId)
mdbi=Channel.create()
mdbi=markdupBamInts
.groupTuple()

mdbi = logChannelContent("Grouped BAM intervals by overall subject/patient ID ",  mdbi)

// group the bams by overall subject/patient id (mergeId)
mdb=Channel.create()
mdb=markdupBam
.groupTuple()

mdb = logChannelContent("Grouped BAMs by overall subject/patient ID ",  mdb)

/*
    Creating target intervals for indel realigner.
    Though VCF indexes are not needed explicitly, we are adding them so they will be linked, and not re-created on the fly.
 */
process CreateIntervals {

	cpus 4
	
	input:
	set mergeId, id, file(mdBam), file(mdBai) from mdbi
	file gf from file(refs["genome_file"]) 
	file gi from file(refs["genome_index"])
	file gd from file(refs["genome_dict"])
	file ki from file(refs["kgindels"])
    file kix from file(refs["kgidx"])
	file mi from file(refs["millsindels"])
	file mix from file(refs["millsidx"])

	output:
	file("${mergeId}.intervals") into intervals

	script:
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


///*
// * realign, use nWayOut to split into tumor/normal again
// */
//process realign {
//
//	input:
//	set mergeId, id, file(mdBam), file(mdBai) from mdb
//	file gf from file(refs["genome_file"])
//	file gi from file(refs["genome_index"])
//	file gd from file(refs["genome_dict"])
//	file ki from file(refs["kgindels"])
//    file kix from file(refs["kgidx"])
//	file mi from file(refs["millsindels"])
//	file mix from file(refs["millsidx"])
//	file intervals from intervals
//
//	output:
//	//set mergeId, id, file("${id}.real.bam") into realBams
//	set mergeId, id, file("${mergeId}.real.bam") into realBams
//
//	script:
//	input = mdBam.collect{"-I $it"}.join(' ')
//
//	"""
//	java -Xmx7g -jar /sw/apps/bioinfo/GATK/3.3.0/GenomeAnalysisTK.jar \
//	-T IndelRealigner \
//	$input \
//	-R $gf \
//	-targetIntervals $intervals \
//	-known $ki \
//	-known $mi \
//	-nWayOut '.real.bam'
//	"""
//
//}
//
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

/*
    Paolo Di Tommaso told that it should be solved by the Channel view() method, but frankly that was not working
    as expected. This simple function makes two channels from one, and logs the content for one (thus consuming that channel)
    and returns with the intact copy for further processing.
 */
def logChannelContent(aMessage, aChannel) {
   
    resChannel = Channel.create()
    logChannel = Channel.create()
    Channel
        .from aChannel
        .separate(resChannel,logChannel) {a -> [a,a] }
    logChannel.subscribe { log.info aMessage + " -- $it" }
    return resChannel
}
