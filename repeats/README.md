Intervals in tiny_GRCh37.list were lifted over from GRCh37 to GRCh38 in the following
way:

```
module load bioinfo-tools samtools/1.4 bwa
samtools faidx /sw/data/uppnex/ToolBox/ReferenceAssemblies/hg38make/bundle/2.8/b37/human_g1k_v37_decoy.fasta $(< tiny.list) > intervals.fasta
bwa mem /sw/data/uppnex/ToolBox/hg38bundle/Homo_sapiens_assembly38.fasta intervals.fasta > intervals.sam
grep -v '^@' intervals.sam | awk '{printf("%s:%d-%d\n", $3, $4, $4+$6-1)}' > tiny-GRCh38.list
```
