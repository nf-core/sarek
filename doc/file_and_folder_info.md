# CAW
Cancer Analysis Workflow Prototype at SciLifeLab

## Project folder structure

The workflow is started for a set of samples taken from the same person (i.e. a cancer patient), identified by "ID".

ID = Individual

Every sample taken from the same ID has a unique identifier. So for one ID there are several Samples

Sample = "Normal", "Tumor 1", "Tumor 2" etc corresponding to all physical samples

The file structure for one ID is shown below. 


![folder_structure.jpg]


## Input Fastq file name conventions

The input fastq files for one ID (Individual) should be organized into one subfolder for every sample, where all fastq files for that sample should be collected.

ID  
    sample1
        sample1_lib_flowcell-index_lane_R1_1000.fastq.gz  
        sample1_lib_flowcell-index_lane_R2_1000.fastq.gz  
        sample1_lib_flowcell-index_lane_R1_1000.fastq.gz  
        sample1_lib_flowcell-index_lane_R2_1000.fastq.gz  
    sample2  
        sample2_lib_flowcell-index_lane_R1_1000.fastq.gz  
        sample2_lib_flowcell-index_lane_R2_1000.fastq.gz  
    sample3  
        sample3_lib_flowcell-index_lane_R1_1000.fastq.gz  
        sample3_lib_flowcell-index_lane_R2_1000.fastq.gz  
        sample3_lib_flowcell-index_lane_R1_1000.fastq.gz  
        sample3_lib_flowcell-index_lane_R2_1000.fastq.gz  



Fastq filename structure:  
sample_lib_flowcell-index_lane_R1_1000.fastq.gz  
and   
sample_lib_flowcell-index_lane_R2_1000.fastq.gz  

Where  
sample = sample id  
lib = indentifier of libaray preparation  
flowcell = identifyer of flow cell for the sequencing run  
lane = identifier of the lane of the sequencing run  
  
Read group information will be parsed from fastq file names according to this:  
RGID = “sample_lib_flowcell-index_lane"  
RGPL = “Illumina"  
PU = sample  
RGLB = lib  