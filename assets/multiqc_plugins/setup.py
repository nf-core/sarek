#!/usr/bin/env python
"""
Setup code for Genomics pipeline MultiQC plugin.
For more information about MultiQC, see http://multiqc.info
"""

from setuptools import setup, find_packages

version = '1.0.0'

setup(
    name = 'multiqc_zymo',
    version = version,
    description = "MultiQC plugins for genomics pipeline",
    packages = find_packages(),
    include_package_data = True,
    install_requires = ['multiqc==1.21'],
    entry_points = {
        'multiqc.templates.v1': [
            'zymo = multiqc_zymo.templates.zymo'
        ],
        'multiqc.hooks.v1': [
            'before_config = multiqc_zymo.utils:plugin_before_config'
        ]
    }
)
