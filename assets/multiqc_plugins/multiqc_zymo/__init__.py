#! /usr/bin/env python

from pkg_resources import get_distribution
from multiqc.utils import config

__version__ = get_distribution("multiqc_zymo").version
config.multiqc_zymo_version = __version__
