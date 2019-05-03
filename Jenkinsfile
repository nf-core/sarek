pipeline {
    agent any

    environment {
        JENKINS_API = credentials('api')
    }

    stages {
        stage('Setup environment') {
            steps {
                sh "docker pull nfcore/sarek:dev"
            }
        }
        stage('Build') {
            steps {
              sh "git clone --single-branch --branch sarek https://github.com/nf-core/test-datasets.git data"
              sh "nextflow run build.nf -profile docker --genome smallGRCh37 --refdir data/reference --outdir references --publishDirMode link -ansi-log false"
              sh "rm -rf work/ references/pipeline_info .nextflow*"
            }
        }
        stage('SampleDir') {
            steps {
                sh "nextflow run main.nf -profile docker --sampleDir data/testdata/tiny/normal --tools HaplotypeCaller,Manta,Strelka --genome smallGRCh37 --igenomes_base references --publishDirMode link -ansi-log false"
                sh "rm -rf work/ .nextflow* results/"
            }
        }
        stage('Multiple') {
            steps {
                sh "nextflow run main.nf -profile docker --sample data/testdata/tsv/tiny-multiple.tsv --tools HaplotypeCaller,Manta,Strelka --genome smallGRCh37 --igenomes_base references --publishDirMode link -ansi-log false"
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
