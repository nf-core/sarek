import os
import shutil
import subprocess
import typing
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Any, List, Optional

import requests
import typing_extensions
from flytekit.core.annotation import FlyteAnnotation
from latch.ldata.path import LPath
from latch.resources.tasks import custom_task, nextflow_runtime_task
from latch.resources.workflow import workflow
from latch.types import metadata
from latch.types.directory import LatchDir, LatchOutputDir
from latch.types.file import LatchFile
from latch_cli.nextflow.utils import _get_execution_name
from latch_cli.nextflow.workflow import get_flag
from latch_cli.services.register.utils import import_module_by_path
from latch_cli.utils import urljoins

meta = Path("latch_metadata") / "__init__.py"
import_module_by_path(meta)
import latch_metadata


@dataclass(frozen=True)
class Sample:
    patient: str
    status: int
    sample: str
    lane: str
    fastq_1: LatchFile
    fastq_2: typing.Optional[LatchFile]

@custom_task(cpu=0.25, memory=0.5, storage_gib=1)
def initialize() -> str:
    token = os.environ.get("FLYTE_INTERNAL_EXECUTION_ID")
    if token is None:
        raise RuntimeError("failed to get execution token")

    headers = {"Authorization": f"Latch-Execution-Token {token}"}

    print("Provisioning shared storage volume... ", end="")
    resp = requests.post(
        "http://nf-dispatcher-service.flyte.svc.cluster.local/provision-storage",
        headers=headers,
        json={
            "storage_gib": 100,
        },
    )
    resp.raise_for_status()
    print("Done.")

    return resp.json()["name"]

def get_flag_defaults(name: str, val: Any, default_val: Optional[Any]):
    if val == default_val or val is None:
        return ""
    else:
        return get_flag(name=name, val=val)

@nextflow_runtime_task(cpu=4, memory=8, storage_gib=100)
def nextflow_runtime(
    pvc_name: str,
    input: typing.List[Sample],
    outdir: typing_extensions.Annotated[LatchDir, FlyteAnnotation({"output": True})],
    wes: typing.Optional[bool],
    intervals: typing.Optional[str],
    no_intervals: typing.Optional[bool],
    tools: typing.Optional[str],
    skip_tools: typing.Optional[str],
    trim_fastq: typing.Optional[bool],
    umi_read_structure: typing.Optional[str],
    save_mapped: typing.Optional[bool],
    save_output_as_bam: typing.Optional[bool],
    use_gatk_spark: typing.Optional[str],
    concatenate_vcfs: typing.Optional[bool],
    only_paired_variant_calling: typing.Optional[bool],
    joint_germline: typing.Optional[bool],
    joint_mutect2: typing.Optional[bool],
    bcftools_annotations: typing.Optional[str],
    bcftools_annotations_tbi: typing.Optional[str],
    bcftools_header_lines: typing.Optional[str],
    dbsnp_vqsr: typing.Optional[str],
    fasta: typing.Optional[LatchFile],
    fasta_fai: typing.Optional[str],
    known_indels_vqsr: typing.Optional[str],
    known_snps: typing.Optional[str],
    known_snps_tbi: typing.Optional[str],
    known_snps_vqsr: typing.Optional[str],
    ngscheckmate_bed: typing.Optional[str],
    snpeff_db: typing.Optional[str],
    snpeff_genome: typing.Optional[str],
    vep_genome: typing.Optional[str],
    vep_species: typing.Optional[str],
    vep_cache_version: typing.Optional[str],
    save_reference: typing.Optional[bool],
    build_only_index: typing.Optional[bool],
    download_cache: typing.Optional[bool],
    igenomes_base: typing.Optional[LatchDir],
    igenomes_ignore: typing.Optional[bool],
    vep_cache: typing.Optional[LatchDir],
    snpeff_cache: typing.Optional[LatchDir],
    email: typing.Optional[str],
    multiqc_title: typing.Optional[str],
    multiqc_methods_description: typing.Optional[str],
    step: str,
    split_fastq: typing.Optional[int],
    nucleotides_per_second: typing.Optional[int],
    aligner: typing.Optional[str],
    vep_custom_args: typing.Optional[str],
    vep_version: typing.Optional[str],
    genome: typing.Optional[str],
) -> None:
    try:
        shared_dir = Path("/nf-workdir")

        ignore_list = [
            "latch",
            ".latch",
            "nextflow",
            ".nextflow",
            "work",
            "results",
            "miniconda",
            "anaconda3",
            "mambaforge",
        ]

        shutil.copytree(
            Path("/root"),
            shared_dir,
            ignore=lambda src, names: ignore_list,
            ignore_dangling_symlinks=True,
            dirs_exist_ok=True,
        )

        cmd = [
            "/root/nextflow",
            "run",
            str(shared_dir / "main.nf"),
            "-work-dir",
            str(shared_dir),
            "-profile",
            "docker",
            "-c",
            "latch.config",


            *get_flag("input", input),
            *get_flag("step", step),
            *get_flag("outdir", outdir),
            *get_flag("split_fastq", split_fastq),
            *get_flag_defaults("wes", wes, False),
            *get_flag("intervals", intervals),
            *get_flag("nucleotides_per_second", nucleotides_per_second),
            *get_flag_defaults("no_intervals", no_intervals, False),
            *get_flag("tools", tools),
            *get_flag("skip_tools", skip_tools),
            *get_flag("trim_fastq", trim_fastq),
            *get_flag("umi_read_structure", umi_read_structure),
            *get_flag("aligner", aligner),
            *get_flag_defaults("save_mapped", save_mapped, False),
            *get_flag_defaults("save_output_as_bam", save_output_as_bam, False),
            *get_flag("use_gatk_spark", use_gatk_spark),
            *get_flag_defaults("concatenate_vcfs", concatenate_vcfs, False),
            *get_flag_defaults("only_paired_variant_calling", only_paired_variant_calling, False),
            *get_flag_defaults("joint_germline", joint_germline, False),
            *get_flag_defaults("joint_mutect2", joint_mutect2, False),
            *get_flag("vep_custom_args", vep_custom_args),
            *get_flag("vep_version", vep_version),
            *get_flag("bcftools_annotations", bcftools_annotations),
            *get_flag("bcftools_annotations_tbi", bcftools_annotations_tbi),
            *get_flag("bcftools_header_lines", bcftools_header_lines),
            *get_flag("genome", genome),
            *get_flag("dbsnp_vqsr", dbsnp_vqsr),
            *get_flag("fasta", fasta),
            *get_flag("fasta_fai", fasta_fai),
            *get_flag("known_indels_vqsr", known_indels_vqsr),
            *get_flag("known_snps", known_snps),
            *get_flag("known_snps_tbi", known_snps_tbi),
            *get_flag("known_snps_vqsr", known_snps_vqsr),
            *get_flag("ngscheckmate_bed", ngscheckmate_bed),
            *get_flag("snpeff_db", snpeff_db),
            *get_flag("snpeff_genome", snpeff_genome),
            *get_flag("vep_genome", vep_genome),
            *get_flag("vep_species", vep_species),
            *get_flag("vep_cache_version", vep_cache_version),
            *get_flag_defaults("save_reference", save_reference, False),
            *get_flag_defaults("build_only_index", build_only_index, False),
            *get_flag("download_cache", download_cache),
            *get_flag("igenomes_ignore", igenomes_ignore),
            *get_flag("igenomes_base", igenomes_base),
            *get_flag("vep_cache", vep_cache),
            *get_flag("snpeff_cache", snpeff_cache),
            *get_flag("email", email),
            *get_flag("multiqc_title", multiqc_title),
            *get_flag("multiqc_methods_description", multiqc_methods_description),
        ]

        print("Launching Nextflow Runtime")
        print(" ".join(cmd))
        print(flush=True)

        env = {
            **os.environ,
            "NXF_HOME": "/root/.nextflow",
            "NXF_OPTS": "-Xms2048M -Xmx8G -XX:ActiveProcessorCount=4",
            "K8S_STORAGE_CLAIM_NAME": pvc_name,
            "NXF_DISABLE_CHECK_LATEST": "true",
        }
        subprocess.run(
            cmd,
            env=env,
            check=True,
            cwd=str(shared_dir),
        )
    finally:
        print()

        nextflow_log = shared_dir / ".nextflow.log"
        if nextflow_log.exists():
            name = _get_execution_name()
            if name is None:
                print("Skipping logs upload, failed to get execution name")
            else:
                remote = LPath(urljoins("latch:///your_log_dir/nf_nf_core_sarek", name, "nextflow.log"))
                print(f"Uploading .nextflow.log to {remote.path}")
                remote.upload_from(nextflow_log)


@workflow(metadata._nextflow_metadata)
def nf_nf_core_sarek(
    input: typing.List[Sample],
    outdir: typing_extensions.Annotated[LatchDir, FlyteAnnotation({"output": True})],
    intervals: typing.Optional[str],
    tools: typing.Optional[str],
    skip_tools: typing.Optional[str],
    umi_read_structure: typing.Optional[str],
    use_gatk_spark: typing.Optional[str],
    bcftools_annotations: typing.Optional[str],
    bcftools_annotations_tbi: typing.Optional[str],
    bcftools_header_lines: typing.Optional[str],
    dbsnp_vqsr: typing.Optional[str],
    fasta: typing.Optional[LatchFile],
    fasta_fai: typing.Optional[str],
    known_indels_vqsr: typing.Optional[str],
    known_snps: typing.Optional[str],
    known_snps_tbi: typing.Optional[str],
    known_snps_vqsr: typing.Optional[str],
    ngscheckmate_bed: typing.Optional[str],
    snpeff_db: typing.Optional[str],
    snpeff_genome: typing.Optional[str],
    vep_genome: typing.Optional[str],
    vep_species: typing.Optional[str],
    vep_cache_version: typing.Optional[str],
    igenomes_base: typing.Optional[LatchDir],
    vep_cache: typing.Optional[LatchDir],
    snpeff_cache: typing.Optional[LatchDir],
    email: typing.Optional[str],
    multiqc_title: typing.Optional[str],
    multiqc_methods_description: typing.Optional[str],
    wes: bool = False,
    no_intervals: bool = False,
    trim_fastq: bool = False,
    save_mapped: bool = False,
    save_output_as_bam: bool = False,
    concatenate_vcfs: bool = False,
    only_paired_variant_calling: bool = False,
    joint_germline: bool = False,
    joint_mutect2:  bool = False,
    save_reference:  bool = False,
    build_only_index: bool = False,
    download_cache: bool = False,
    igenomes_ignore: bool = False,
    step: str = "mapping",
    split_fastq: typing.Optional[int] = 50000000,
    nucleotides_per_second: typing.Optional[int] = 200000,
    aligner: typing.Optional[str] = "bwa-mem",
    vep_custom_args: typing.Optional[
        str
    ] = "--everything --filter_common --per_gene --total_length --offline --format vcf",
    vep_version: typing.Optional[str] = "111.0-0",
    genome: typing.Optional[str] = "GATK.GRCh38",
) -> None:
    """
    nf-core/sarek

    Sample Description
    """

    pvc_name: str = initialize()
    nextflow_runtime(
        pvc_name=pvc_name,
        input=input,
        step=step,
        outdir=outdir,
        split_fastq=split_fastq,
        wes=wes,
        intervals=intervals,
        nucleotides_per_second=nucleotides_per_second,
        no_intervals=no_intervals,
        tools=tools,
        skip_tools=skip_tools,
        trim_fastq=trim_fastq,
        umi_read_structure=umi_read_structure,
        aligner=aligner,
        save_mapped=save_mapped,
        save_output_as_bam=save_output_as_bam,
        use_gatk_spark=use_gatk_spark,
        concatenate_vcfs=concatenate_vcfs,
        only_paired_variant_calling=only_paired_variant_calling,
        joint_germline=joint_germline,
        joint_mutect2=joint_mutect2,
        vep_custom_args=vep_custom_args,
        vep_version=vep_version,
        bcftools_annotations=bcftools_annotations,
        bcftools_annotations_tbi=bcftools_annotations_tbi,
        bcftools_header_lines=bcftools_header_lines,
        genome=genome,
        dbsnp_vqsr=dbsnp_vqsr,
        fasta=fasta,
        fasta_fai=fasta_fai,
        known_indels_vqsr=known_indels_vqsr,
        known_snps=known_snps,
        known_snps_tbi=known_snps_tbi,
        known_snps_vqsr=known_snps_vqsr,
        ngscheckmate_bed=ngscheckmate_bed,
        snpeff_db=snpeff_db,
        snpeff_genome=snpeff_genome,
        vep_genome=vep_genome,
        vep_species=vep_species,
        vep_cache_version=vep_cache_version,
        save_reference=save_reference,
        build_only_index=build_only_index,
        download_cache=download_cache,
        igenomes_base=igenomes_base,
        igenomes_ignore=igenomes_ignore,
        vep_cache=vep_cache,
        snpeff_cache=snpeff_cache,
        email=email,
        multiqc_title=multiqc_title,
        multiqc_methods_description=multiqc_methods_description,
    )
