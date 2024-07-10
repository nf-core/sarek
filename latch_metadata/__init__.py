
from latch.types.metadata import (
    NextflowMetadata,
    LatchAuthor,
    NextflowRuntimeResources
)
from latch.types.directory import LatchDir

from .parameters import generated_parameters

NextflowMetadata(
    display_name='nf-core/sarek',
    author=LatchAuthor(
        name="LatchBio",
    ),
    parameters=generated_parameters,
    runtime_resources=NextflowRuntimeResources(
        cpus=4,
        memory=8,
        storage_gib=500,
    ),
    log_dir=LatchDir("latch:///nextflow_sarek_logs"),
)

