#!/bin/bash

nextflow run . -profile docker_skadi --test --steps preprocessing,MultiQC
nextflow run . -profile docker_skadi --test --steps realign,MultiQC
nextflow run . -profile docker_skadi --test --steps recalibrate,MultiQC
nextflow run . -profile docker_skadi --test --steps skipPreprocessing,HaplotypeCaller,MuTect1,MuTect2
