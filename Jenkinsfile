pipeline {
    agent any

    environment {
        JENKINS_API = credentials('api')
    }

    stages {
        stage('Docker setup') {
            steps {
                 sh "docker pull nfcore/sarek:dev"
                 sh "docker tag nfcore/sarek:dev nfcore/sarek:2.5.2"
                 sh "docker pull nfcore/sareksnpeff:dev.GRCh37"
                 sh "docker tag nfcore/sareksnpeff:dev.GRCh37 nfcore/sareksnpeff:2.5.2.GRCh37"
                 sh "docker pull nfcore/sarekvep:dev.GRCh37"
                 sh "docker tag nfcore/sarekvep:dev.GRCh37 nfcore/sarekvep:2.5.2.GRCh37"
            }
        }
        stage('Annotation') {
            steps {
                sh "nextflow run . -profile test_annotation,kraken --verbose --tools snpeff,vep,merge"
            }
        }
        stage('Germline') {
            steps {
                sh "rm -rf data/"
                sh "git clone --single-branch --branch sarek https://github.com/nf-core/test-datasets.git data"
                sh "nextflow run . -profile test,kraken --input data/testdata/tiny/normal"
                sh "nextflow run . -profile test,kraken --input=false --step recalibrate -resume"
                sh "nextflow run . -profile test,kraken --input=false --step variantCalling"
                sh "rm -rf data/"
            }
        }
        stage('Minimal') {
            steps {
                sh "nextflow run . -profile test,kraken --skipQC all --verbose --genome smallerGRCh37 --no_intervals --tools Manta,mpileup,Strelka"
                sh "nextflow run . -profile test,kraken --skipQC all --verbose --genome smallerGRCh37 --tools Manta,mpileup,Strelka"
                sh "nextflow run . -profile test,kraken --skipQC all --verbose --genome minimalGRCh37 --no_intervals --tools Manta,mpileup,Strelka"
                sh "nextflow run . -profile test,kraken --skipQC all --verbose --genome minimalGRCh37 --tools Manta,mpileup,Strelka"
            }
        }
        stage('Profile') {
            steps {
                sh "nextflow run . -profile test_splitfastq,kraken --verbose"
                sh "nextflow run . -profile test_targeted,kraken --verbose"
            }
        }
        stage('Tools') {
            steps {
                sh "nextflow run . -profile test_tool,kraken --verbose --tools Haplotypecaller,Freebayes,Manta,mpileup,Mutect2,Strelka"
            }
        }
    }

    post {
        failure {
            script {
                def response = sh(script: "curl -u ${JENKINS_API_USR}:${JENKINS_API_PSW} ${BUILD_URL}/consoleText", returnStdout: true).trim().replace('\n', '<br>')
                def comment = pullRequest.comment("## :rotating_light: Buil log output:<br><summary><details>${response}</details></summary>")
            }
        }
    }
}
