pipeline {
    agent any

    environment {
        JENKINS_API = credentials('api')
    }

    stages {
        stage('Setup environment') {
            steps {
                sh "docker pull nfcore/sarek:dev"
                sh "docker tag nfcore/sarek:dev"
                sh "docker pull nfcore/sareksnpeff:dev.GRCh37"
                sh "docker tag nfcore/sareksnpeff:dev.GRCh37 nfcore/sareksnpeff:dev.smallGRCh37"
                sh "docker pull nfcore/sarekvep:dev.GRCh37"
                sh "docker tag nfcore/sarekvep:dev.GRCh37 nfcore/sarekvep:dev.smallGRCh37"
            }
        }
        stage('Build') {
            steps {
                sh "rm -rf data"
                sh "git clone --single-branch --branch sarek https://github.com/nf-core/test-datasets.git data"
                sh "nextflow run build.nf -profile docker --genome smallGRCh37 --refdir data/reference --outdir references --publishDirMode link -ansi-log false"
                sh "rm -rf work/ references/pipeline_info .nextflow*"
            }
        }
        stage('Germline from directory') {
            steps {
                sh "nextflow run main.nf -profile docker --sample data/testdata/tiny/normal --tools HaplotypeCaller,Manta,Strelka,snpEff,VEP,merge --genome smallGRCh37 --igenomes_base references --publishDirMode link -ansi-log false"
                sh "rm -rf work/ .nextflow* results/"
            }
        }
        stage('Somatic multiple') {
            steps {
                sh "nextflow run main.nf -profile docker --sample data/testdata/tsv/tiny-multiple.tsv --tools HaploTypeCaller,Manta,Strelka,MuTecT2,FreeBayes,snpEff,VEP,merge --genome smallGRCh37 --igenomes_base references --publishDirMode link -ansi-log false"
                sh "rm -rf work/ .nextflow* results/"
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
