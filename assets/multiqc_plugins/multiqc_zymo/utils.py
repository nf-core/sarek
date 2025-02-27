#!/usr/bin/env python
"""
We can add any custom Python functions here and call them
using the setuptools plugin hooks.
"""

from __future__ import print_function
from pkg_resources import get_distribution
import logging

from multiqc.utils import report, util_functions, config

# Initialise the main MultiQC logger
log = logging.getLogger('multiqc')

# Save this plugin's version number (defined in setup.py) to the MultiQC config
config.multiqc_zymo_version = get_distribution("multiqc_zymo").version

# Add default config options that can be overriden by user config
def plugin_before_config():
    
    # Use the zymo template by default
    config.template = 'zymo'
