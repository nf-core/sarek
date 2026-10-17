process GPU_SMOKE {
    accelerator 1
    arch 'linux/amd64'
    container 'docker.io/nvidia/cuda:12.4.1-base-ubuntu22.04'
    publishDir params.outdir, mode: 'copy'

    output:
    path 'gpu.txt'

    script:
    '''
    nvidia-smi -L > gpu.txt
    grep -q '^GPU [0-9]' gpu.txt
    nvidia-smi --query-gpu=name,uuid,driver_version --format=csv,noheader >> gpu.txt
    uname -m >> gpu.txt
    cat gpu.txt
    '''
}

process ARM_SMOKE {
    arch 'linux/arm64'
    container 'docker.io/library/ubuntu:24.04'
    publishDir params.outdir, mode: 'copy'

    output:
    path 'arm.txt'

    script:
    '''
    uname -m > arm.txt
    test "$(cat arm.txt)" = aarch64
    echo "arithmetic=$((19 * 23))" >> arm.txt
    grep -qx 'arithmetic=437' arm.txt
    cat arm.txt
    '''
}

workflow {
    if (!params.outdir) {
        error('Specify --outdir for smoke evidence')
    }
    if (params.target == 'gpu') {
        GPU_SMOKE()
    }
    else if (params.target == 'arm') {
        ARM_SMOKE()
    }
    else {
        error('Specify --target gpu or --target arm')
    }
}
