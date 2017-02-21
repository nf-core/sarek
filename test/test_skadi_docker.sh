#!/bin/bash

nextflow run . -profile docker_skadi --test
nextflow run . -profile docker_skadi --testRealign
nextflow run . -profile docker_skadi --testRecalibrate
nextflow run . -profile docker_skadi --testCoreVC
