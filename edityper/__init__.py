#!/usr/bin/env python

"""
EdiTyper
===

A Python program for aligning CRISPR-edited RNA-seq data

Available submodules
--------------------
alignment
    Functions for aligning reads to a reference
analysis
    Tools for analysing read alignment
plots
    Utilities for plotting analysed reads and alignment quality scores
quality_control
    Tools for ensuring reference and template are good to use
recnw
    The aligner itself
sam
    Support for SAM/BAM output of aligned reads
toolkit
    Miscellaneous functions that serve as helpers throughout EdiTyper
"""

try:
    __EDITYPER_SETUP__
except NameError:
    from . import alignment
    from . import analysis
    from . import plots
    from . import quality_control
    from . import sam
    from . import toolkit


__all__ = ['alignment', 'analysis', 'quality_control', 'toolkit', 'sam', 'recnw']
__version__ = '0.2.0'
