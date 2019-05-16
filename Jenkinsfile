pipeline {
    agent any

    environment {
        JENKINS_API = credentials('api')
    }

    stages {
        stage('Setup environment') {
            steps {
                sh "./bin/download_docker.sh -t ALL"
            }
        }
        stage('Build') {
            steps {
                sh "rm -rf references/"
                sh "./bin/build_reference.sh --test ALL --build"
                sh "rm -rf .nextflow* references/pipeline_info work/"
            }
        }
        stage('Germline') {
            steps {
                sh "rm -rf data/"
                sh "git clone --single-branch --branch sarek https://github.com/nf-core/test-datasets.git data"
                sh "./bin/run_tests.sh --test GERMLINE"
                sh "rm -rf .nextflow* data/ results/ work/"
            }
        }
        stage('Somatic') {
            steps {
                sh "./bin/run_tests.sh --test SOMATIC"
                sh "rm -rf .nextflow* results/ work/"
            }
        }
        stage('targeted') {
            steps {
                sh "./bin/run_tests.sh --test TARGETED"
                sh "rm -rf .nextflow* results/ work/"
            }
        }
        stage('Annotation') {
            steps {
                sh "./bin/run_tests.sh --test ANNOTATEALL"
                sh "rm -rf .nextflow* results/ work/"
            }
        }
        stage('Multiple') {
            steps {
                sh "./bin/run_tests.sh --test MULTIPLE"
                sh "rm -rf .nextflow* results/ work/"
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
