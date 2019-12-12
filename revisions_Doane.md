# Notes
    - I made changes to main.nf in a private branch to deal with issues below, as my solutions may break on other systems. I added comments to main.nf, but otherwise this PR should only be adding the cluster config files, conf/panda.config and conf/nygc.config.  See comments below for changes I made on my private branch.
    
# Suggested revisions

    - Some processes that use pipes are using more than the requested resources.  See comments in main.nf for `process MapReads`.  Commands that are piped together will run concurrently, and the resources used will be the sum of the resources for each piped command. So `bwa mem -t 2 ... | samtools sort -@2 ...` will use 4 cpus, not 2.  The sample applies for memory.  I made a private branch with specific cpu and memory values as needed for my use case (WGS), but you likely need a more generalizeable solution. 

    - Feature request: option to set temp directory path or environment variable. Several processes make use of a temp directory using `/tmp`.  The compute nodes at Weill Cornell and New York Genome Center (NYGC) have limited space in `/tmp` and use the env variable `$TMPDIR`, which specifies a local scratch disc.
    
    - BaseRecalibrator: reduce memory resource request.  The bam files are spllit up for parallel processing, so using max memory was way overkill on my system and caused jobs to delay in queue. Something like  10G is probably more than enough.  In a private branch, I changed it to "memory_singleCPU_2_task" and this worked well for WGS data."
    
    - Exit code 140 for Univa grid engine could be added to base.config to resubmit job.
    
    - Issue: some ApplyBQSR jobs are failing with exit code 141.  Cause: In `base.config`, setting pipefail option causes jobs to sometimes fail with exit code 141.  Known issue with pipefail and exit 141.  See https://gitter.im/nextflow-io/nextflow/archives/2016/04/12 
    Solution for now, is changing `process.shell = ['/bin/bash', '-euo', 'pipefail']` to `process.shell = ['/bin/bash', '-eu']`
