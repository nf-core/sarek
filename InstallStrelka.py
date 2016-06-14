#!/usr/bin/env python

from subprocess import call

from os import chdir
from os import mkdir
from os import path
import argparse

parser = argparse.ArgumentParser(description='Install Strelka.')
parser.add_argument('pathway', help='pathway to install Strelka')
args = parser.parse_args()
pathway = args.pathway

version='0.1'
name='InstallStrelka'

print "Installing dependancies for CAW"
print "Installing Strelka"

call(["wget", "ftp://strelka:%27%27@ftp.illumina.com/v1-branch/v1.0.14/strelka_workflow-1.0.14.tar.gz"])

call(["tar", "-xzf", "strelka_workflow-1.0.14.tar.gz"])

chdir("strelka_workflow-1.0.14")


prefix = "--prefix="+pathway

call(["./configure", prefix])

call(["make"])

demo = pathway+"/bin/demo/run_demo.bash"
call([demo])
