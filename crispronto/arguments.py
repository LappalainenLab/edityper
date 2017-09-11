#!/usr/bin/env python

'''Set the arguments for the CRISPR alingment and analysis pipeline'''

import os
import sys
import argparse

_OUTDIR_DEFAULT = os.path.join(os.getcwd(), 'output')
_MODE_DEFAULT = 'SNP'
_PROJECT_DEFAULT = 'crispronto'
_VERBOSITY_DEFAULT = 'info'
_VERBOSITY_LEVELS = (
    'debug',
    'info',
    'warning',
    'error',
    'critical'
)
_VALID_PLATFORMS = (
    'CAPILLARY',
    'LS454',
    'ILLUMINA',
    'SOLID',
    'HELICOS',
    'IONTORRENT',
    'ONT',
    'PACBIO'
)

def make_argument_parser():
    """Create an argument parser"""
    parser = argparse.ArgumentParser(
        add_help=False
    )
    #   Arguments for verbosity and logging
    parser.add_argument( # Verbosity
        '-v',
        '--verbosity',
        dest='verbosity',
        type=str.lower,
        choices=_VERBOSITY_LEVELS,
        default=_VERBOSITY_DEFAULT,
        required=False,
        metavar='VERBOSITY',
        help="Set the verbosity level, choose from '%s'; defaults to '%s'" % ("', '".join(_VERBOSITY_LEVELS), _VERBOSITY_DEFAULT)
    )
    parser.add_argument( # Log file
        '-f',
        '--logfile',
        dest='logfile',
        type=argparse.FileType(mode='w'),
        default=sys.stderr,
        required=False,
        metavar='LOG FILE',
        help="Specify a file for the log messages, defaults to stderr"
    )
    align_opts = parser.add_argument_group( # Alignment options
        title='alignment arguments',
        description='Set parameters for alignment'
    )
    align_opts.add_argument( # Analysis mode
        '-m',
        '--analysis-mode',
        dest='analysis_mode',
        type=str,
        choices=('SNP', 'SNP+PAM'),
        default=_MODE_DEFAULT,
        required=False,
        metavar='ANALYSIS MODE',
        # help="Set the analysis mode, choose from 'SNP' or 'SNP+PAM', defaults to '%s'" % _MODE_DEFAULT
        help=argparse.SUPPRESS
    )
    align_opts.add_argument( # p-value threshold
        '-p',
        '--pvalue-threshold',
        dest='pvalue_threshold',
        type=float,
        default=1e-3,
        required=False,
        metavar='PVALUE THRESHOLD',
        help="Set the p-value threshold, defaults to '1e-3'"
    )
    align_opts.add_argument( # Gap opening penaltiy
        '-g',
        '--gap-opening',
        dest='gap_opening',
        type=int,
        default=8,
        required=False,
        metavar='GAP OPENING PENALTY',
        help="Set the gap opeing penalty, defaults to '8'"
    )
    align_opts.add_argument( # Gap extension penalty
        '-e',
        '--gap-extension',
        dest='gap_extension',
        type=int,
        default=1,
        required=False,
        metavar='GAP EXTENSION PENALTY',
        help="Set the gap extension penalty, defaults to '1'"
    )
    align_opts.add_argument(
        '-n',
        '--num-cores',
        dest='num_cores',
        type=int,
        default=None,
        required=False,
        metavar='NUM CORES',
        help="How many cores should we use for multiprocessing? Defaults to the number of FASTQ files provided or the number of cores available on the system, whichever is lower"
    )
    ref_opts = parser.add_argument_group( # Reference and template sequences
        title='reference arguments',
        description='Provide FASTA files for the reference and template sequences'
    )
    ref_opts.add_argument( # Reference sequence
        '-r',
        '--reference-sequence',
        dest='reference',
        type=str,
        default=None,
        required=True,
        metavar='REFERENCE SEQUENCE',
        help="Choose a reference FASTA file"
    )
    ref_opts.add_argument( # Template sequence
        '-t',
        '--template-sequence',
        dest='template',
        type=str,
        default=None,
        required=True,
        metavar='TEMPLATE SEQUENCE',
        help="Choose a template FASTA file"
    )
    in_opts = parser.add_argument_group( # Input FASTQ options
        title='input arguments',
        description='Provide either a single FASTQ file or a list of FASTQ files. Note: we currently do NOT support paired-end FASTQ files'
    )
    in_files = in_opts.add_mutually_exclusive_group(required=True)
    in_files.add_argument( # Are we providing an input file?
        '-i',
        '--input-file',
        dest='input_file',
        type=str,
        default=None,
        metavar='INPUT FILE',
        help="Provide a single input file, mutually exclusive with '-l | --sample-list'"
    )
    in_files.add_argument( # Are we using a sample list?
        '-l',
        '--sample-list',
        dest='sample_list',
        type=str,
        default=None,
        metavar='SAMPLE LIST',
        help="Provdide a sample list, with each sample on its own line, mutually exclusive with '-i | --input-file'"
    )
    out_opts = parser.add_argument_group( # Output options
        title='output arguments',
        description='Provide an output directory and project name. All files, including the configuration file, will be placed in the output directory with a basename of the project name.'
    )
    out_opts.add_argument( # Output directory
        '-d',
        '--output-directory',
        dest='outdirectory',
        type=str,
        default=_OUTDIR_DEFAULT,
        required=False,
        metavar='OUTPUT DIRECTORY',
        help="Choose where all output files are to be stored; defaults to '%s'" % _OUTDIR_DEFAULT
    )
    out_opts.add_argument( # Project name
        '-j',
        '--project',
        dest='project',
        type=str,
        default=_PROJECT_DEFAULT,
        required=False,
        metavar='PROJECT',
        help="Provide a project name to be used as the basename for all output files; defaults to '%s'" % _PROJECT_DEFAULT
    )
    rg_opts = parser.add_argument_group( # Arguments for read groups
        title='read group arguments',
        description="Provide extra read group information. Note: this information will be applied to ALL read groups"
    )
    rg_opts.add_argument( # Read center
        '-rc',
        '--read-center',
        dest='read_center',
        type=str,
        default=None,
        required=False,
        metavar='READ CENTER',
        help="Name of the sequencing center (@RG CN)"
    )
    rg_opts.add_argument( # Read library
        '-rl',
        '--read-library',
        dest='read_library',
        type=str,
        default=None,
        required=False,
        metavar='READ LIBRARY',
        help="Library (@RG LB)"
    )
    rg_opts.add_argument( # Read platform
        '-rp',
        '--read-platform',
        dest='read_platform',
        type=str,
        default=None,
        required=False,
        choices=_VALID_PLATFORMS,
        metavar='READ PLATFORM',
        help="Platform technology for sequencing, choose from '%s' (@RG PL)" % "', '".join(_VALID_PLATFORMS)
    )
    rg_opts.add_argument( # Read sample
        '-rs',
        '--read-sample',
        dest='read_sample',
        type=str,
        default=None,
        required=False,
        metavar='READ SAMPLE',
        help="Sample, use pool name when a pool is being sequenced (@RG SM)"
    )
    suppression = parser.add_argument_group( # Arguments for suppressing outputs
        title='suppression arguments',
        description="Choose to suppress some or all of the output files"
    )
    suppression.add_argument( # Suppress SAM output
        '--suppress-sam',
        dest='suppress_sam',
        action='store_true',
        required=False,
        help="Suppress SAM output"
    )
    suppression.add_argument( # Suppress .events table
        '--suppress-events',
        dest='suppress_events',
        action='store_true',
        required=False,
        help="Suppress events output"
    )
    suppression.add_argument( # Suppress .classifications table
        '--suppress-classification',
        dest='suppress_classification',
        action='store_true',
        required=False,
        help="Suppress read classification"
    )
    suppression.add_argument( # Suppress both .events and .classifications tables
        '--suppress-tables',
        dest='suppress_tables',
        action='store_true',
        required=False,
        help="Suppress both events and read classification tables"
    )
    suppression.add_argument( # Suppress plots
        '--suppress-plots',
        dest='suppress_plots',
        action='store_true',
        required=False,
        help="Suppress plots"
    )
    parser.add_argument(
        '--profile',
        dest='profile',
        action='store_true',
        required=False,
        # help="Profile CRISPRonto"
        help=argparse.SUPPRESS
    )
    parser.add_argument(
        '--xkcd',
        dest='xkcd',
        action='store_true',
        required=False,
        help=argparse.SUPPRESS
    )
    parser.add_argument(
        '-h',
        '--help',
        action='store_true',
        required=False,
        help=argparse.SUPPRESS
    )
    return parser
