
from dataclasses import dataclass
import typing
import typing_extensions

from flytekit.core.annotation import FlyteAnnotation

from latch.types.metadata import NextflowParameter
from latch.types.file import LatchFile
from latch.types.directory import LatchDir, LatchOutputDir

# Import these into your `__init__.py` file:
#
# from .parameters import generated_parameters

@dataclass(frozen=True)
class Sample:
    patient: str
    lane: int
    status: int
    sample: str
    sex: str
    fastq_1: LatchFile
    fastq_2: LatchFile


generated_parameters = {
    'input': NextflowParameter(
        type=typing.List[Sample],
        samplesheet=True,
        samplesheet_type='csv',
        section_title='Input/output options',
        description='Path to comma-separated file containing information about the samples in the experiment.',
    ),
    'step': NextflowParameter(
        type=str,
        default='mapping',
        section_title=None,
        description='Starting step',
    ),
    'outdir': NextflowParameter(
        type=typing_extensions.Annotated[LatchDir, FlyteAnnotation({'output': True})],
        default=None,
        section_title=None,
        description='The output directory where the results will be saved. You have to use absolute paths to storage on Cloud infrastructure.',
    ),
    'split_fastq': NextflowParameter(
        type=typing.Optional[int],
        default=50000000,
        section_title='Main options',
        description='Specify how many reads each split of a FastQ file contains. Set 0 to turn off splitting at all.',
    ),
    'wes': NextflowParameter(
        type=typing.Optional[bool],
        default=None,
        section_title=None,
        description='Enable when exome or panel data is provided.',
    ),
    'intervals': NextflowParameter(
        type=typing.Optional[str],
        default=None,
        section_title=None,
        description='Path to target bed file in case of whole exome or targeted sequencing or intervals file.',
    ),
    'nucleotides_per_second': NextflowParameter(
        type=typing.Optional[int],
        default=200000,
        section_title=None,
        description='Estimate interval size.',
    ),
    'no_intervals': NextflowParameter(
        type=typing.Optional[bool],
        default=None,
        section_title=None,
        description='Disable usage of intervals.',
    ),
    'tools': NextflowParameter(
        type=typing.Optional[str],
        default=None,
        section_title=None,
        description='Tools to use for duplicate marking, variant calling and/or for annotation.',
    ),
    'skip_tools': NextflowParameter(
        type=typing.Optional[str],
        default=None,
        section_title=None,
        description='Disable specified tools.',
    ),
    'trim_fastq': NextflowParameter(
        type=typing.Optional[bool],
        default=None,
        section_title='FASTQ Preprocessing',
        description='Run FastP for read trimming',
    ),
    'umi_read_structure': NextflowParameter(
        type=typing.Optional[str],
        default=None,
        section_title=None,
        description='Specify UMI read structure',
    ),
    'aligner': NextflowParameter(
        type=typing.Optional[str],
        default='bwa-mem',
        section_title='Preprocessing',
        description='Specify aligner to be used to map reads to reference genome.',
    ),
    'save_mapped': NextflowParameter(
        type=typing.Optional[bool],
        default=None,
        section_title=None,
        description='Save mapped files.',
    ),
    'save_output_as_bam': NextflowParameter(
        type=typing.Optional[bool],
        default=None,
        section_title=None,
        description='Saves output from mapping (if `--save_mapped`), Markduplicates & Baserecalibration as BAM file instead of CRAM',
    ),
    'use_gatk_spark': NextflowParameter(
        type=typing.Optional[str],
        default=None,
        section_title=None,
        description='Enable usage of GATK Spark implementation for duplicate marking and/or base quality score recalibration',
    ),
    'concatenate_vcfs': NextflowParameter(
        type=typing.Optional[bool],
        default=None,
        section_title='Variant Calling',
        description='Option for concatenating germline vcf-files.',
    ),
    'only_paired_variant_calling': NextflowParameter(
        type=typing.Optional[bool],
        default=None,
        section_title=None,
        description='If true, skips germline variant calling for matched normal to tumor sample. Normal samples without matched tumor will still be processed through germline variant calling tools.',
    ),
    'joint_germline': NextflowParameter(
        type=typing.Optional[bool],
        default=None,
        section_title=None,
        description='Turn on the joint germline variant calling for GATK haplotypecaller',
    ),
    'joint_mutect2': NextflowParameter(
        type=typing.Optional[bool],
        default=None,
        section_title=None,
        description='Runs Mutect2 in joint (multi-sample) mode for better concordance among variant calls of tumor samples from the same patient. Mutect2 outputs will be stored in a subfolder named with patient ID under `variant_calling/mutect2/` folder. Only a single normal sample per patient is allowed. Tumor-only mode is also supported.',
    ),
    'vep_custom_args': NextflowParameter(
        type=typing.Optional[str],
        default='--everything --filter_common --per_gene --total_length --offline --format vcf',
        section_title='Annotation',
        description='Add an extra custom argument to VEP.',
    ),
    'vep_version': NextflowParameter(
        type=typing.Optional[str],
        default='111.0-0',
        section_title=None,
        description='Should reflect the VEP version used in the container.',
    ),
    'bcftools_annotations': NextflowParameter(
        type=typing.Optional[str],
        default=None,
        section_title=None,
        description='A vcf file containing custom annotations to be used with bcftools annotate. Needs to be bgzipped.',
    ),
    'bcftools_annotations_tbi': NextflowParameter(
        type=typing.Optional[str],
        default=None,
        section_title=None,
        description='Index file for `bcftools_annotations`',
    ),
    'bcftools_header_lines': NextflowParameter(
        type=typing.Optional[str],
        default=None,
        section_title=None,
        description='Text file with the header lines of `bcftools_annotations`',
    ),
    'genome': NextflowParameter(
        type=typing.Optional[str],
        default='GATK.GRCh38',
        section_title='Reference genome options',
        description='Name of iGenomes reference.',
    ),
    'dbsnp_vqsr': NextflowParameter(
        type=typing.Optional[str],
        default=None,
        section_title=None,
        description='label string for VariantRecalibration (haplotypecaller joint variant calling)',
    ),
    'fasta': NextflowParameter(
        type=typing.Optional[LatchFile],
        default=None,
        section_title=None,
        description='Path to FASTA genome file.',
    ),
    'fasta_fai': NextflowParameter(
        type=typing.Optional[str],
        default=None,
        section_title=None,
        description='Path to FASTA reference index.',
    ),
    'known_indels_vqsr': NextflowParameter(
        type=typing.Optional[str],
        default=None,
        section_title=None,
        description='If you use AWS iGenomes, this has already been set for you appropriately.\n\n1st label string for VariantRecalibration (haplotypecaller joint variant calling)',
    ),
    'known_snps': NextflowParameter(
        type=typing.Optional[str],
        default=None,
        section_title=None,
        description='If you use AWS iGenomes, this has already been set for you appropriately.\n\nPath to known snps file.',
    ),
    'known_snps_tbi': NextflowParameter(
        type=typing.Optional[str],
        default=None,
        section_title=None,
        description='Path to known snps file snps.',
    ),
    'known_snps_vqsr': NextflowParameter(
        type=typing.Optional[str],
        default=None,
        section_title=None,
        description='If you use AWS iGenomes, this has already been set for you appropriately.\n\nlabel string for VariantRecalibration (haplotypecaller joint variant calling)',
    ),
    'ngscheckmate_bed': NextflowParameter(
        type=typing.Optional[str],
        default=None,
        section_title=None,
        description='Path to SNP bed file for sample checking with NGSCheckMate',
    ),
    'snpeff_db': NextflowParameter(
        type=typing.Optional[str],
        default=None,
        section_title=None,
        description='snpEff DB version.',
    ),
    'snpeff_genome': NextflowParameter(
        type=typing.Optional[str],
        default=None,
        section_title=None,
        description='snpEff genome.',
    ),
    'vep_genome': NextflowParameter(
        type=typing.Optional[str],
        default=None,
        section_title=None,
        description='VEP genome.',
    ),
    'vep_species': NextflowParameter(
        type=typing.Optional[str],
        default=None,
        section_title=None,
        description='VEP species.',
    ),
    'vep_cache_version': NextflowParameter(
        type=typing.Optional[str],
        default=None,
        section_title=None,
        description='VEP cache version.',
    ),
    'save_reference': NextflowParameter(
        type=typing.Optional[bool],
        default=None,
        section_title=None,
        description='Save built references.',
    ),
    'build_only_index': NextflowParameter(
        type=typing.Optional[bool],
        default=None,
        section_title=None,
        description='Only built references.',
    ),
    'download_cache': NextflowParameter(
        type=typing.Optional[bool],
        default=None,
        section_title=None,
        description='Download annotation cache.',
    ),
    'igenomes_base': NextflowParameter(
        type=typing.Optional[LatchDir],
        default=None,
        section_title=None,
        description='Directory / URL base for iGenomes references.',
    ),
    'igenomes_ignore': NextflowParameter(
        type=typing.Optional[bool],
        default=None,
        section_title=None,
        description='Do not load the iGenomes reference config.',
    ),
    'vep_cache': NextflowParameter(
        type=typing.Optional[LatchDir],
        default=None,
        section_title=None,
        description='Path to VEP cache.',
    ),
    'snpeff_cache': NextflowParameter(
        type=typing.Optional[LatchDir],
        default=None,
        section_title=None,
        description='Path to snpEff cache.',
    ),
    'email': NextflowParameter(
        type=typing.Optional[str],
        default=None,
        section_title='Generic options',
        description='Email address for completion summary.',
    ),
    'multiqc_title': NextflowParameter(
        type=typing.Optional[str],
        default=None,
        section_title=None,
        description='MultiQC report title. Printed as page header, used for filename if not otherwise specified.',
    ),
    'multiqc_methods_description': NextflowParameter(
        type=typing.Optional[str],
        default=None,
        section_title=None,
        description='Custom MultiQC yaml file containing HTML including a methods description.',
    ),
}

