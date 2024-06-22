from latch.types.directory import LatchDir
from latch.types.metadata import (LatchAuthor, NextflowMetadata,
                                  NextflowRuntimeResources)

from .parameters import generated_parameters

# from .parameters import flow, generated_parameters


NextflowMetadata(
    display_name="nf-core/sarek",
    author=LatchAuthor(
        name="nf-core",
    ),
    parameters=generated_parameters,
    # flow=flow,
    runtime_resources=NextflowRuntimeResources(
        cpus=4,
        memory=8,
        storage_gib=100,
    ),
    log_dir=LatchDir("latch:///your_log_dir"),
)
