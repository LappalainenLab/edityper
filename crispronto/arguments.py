#!/usr/bin/env python

# TODO: Add argument for genomic location of reference
#   Use in SAM output
#   BED input?

'''Set the arguments for the CRISPR alingment and analysis pipeline'''

import os
import sys
import argparse
import multiprocessing

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

def _num_cores(value):
    try:
        value = int(value)
    except ValueError:
        raise argparse.ArgumentTypeError("Must pass an integer value")
    try:
        assert value < multiprocessing.cpu_count()
    except AssertionError:
        raise argparse.ArgumentTypeError("Cannot have more jobs (%s) than cores available (%s)" % (value, multiprocessing.cpu_count()))
    except NotImplementedError:
        pass
    return value


def _validate_bam_index(value):
    if value is False or value in ('bai', 'csi'):
        return value
    else:
        raise argparse.ArgumentTypeError("Can only pass '--bam' or '--bam csi', not '--bam %s'" % value)


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
        metavar='verbosity',
        help="Set the verbosity level, choose from '%s'; defaults to '%s'" % ("', '".join(_VERBOSITY_LEVELS), _VERBOSITY_DEFAULT)
    )
    parser.add_argument( # Log file
        '-f',
        '--logfile',
        dest='logfile',
        type=argparse.FileType(mode='w'),
        default=sys.stderr,
        required=False,
        metavar='log file',
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
        metavar='analysis mode',
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
        metavar='pvalue threshold',
        help="Set the p-value threshold, defaults to '1e-3'"
    )
    align_opts.add_argument( # Gap opening penaltiy
        '-g',
        '--gap-opening',
        dest='gap_opening',
        type=int,
        default=8,
        required=False,
        metavar='gap opening penalty',
        help="Set the gap opeing penalty, defaults to '8'"
    )
    align_opts.add_argument( # Gap extension penalty
        '-e',
        '--gap-extension',
        dest='gap_extension',
        type=int,
        default=1,
        required=False,
        metavar='gap extension penalty',
        help="Set the gap extension penalty, defaults to '1'"
    )
    align_opts.add_argument( # Number of cores
        '--parallel',
        dest='num_cores',
        # type=int,
        type=_num_cores,
        const=None,
        default=1,
        nargs='?',
        required=False,
        # metavar='num cores',
        metavar='num jobs',
        help="Run %(prog)s in parallel. If passed, can optionally specify the number of jobs to run at once."
        # help="How many cores should we use for multiprocessing? Defaults to the number of FASTQ files provided or the number of cores available on the system, whichever is lower"
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
        metavar='reference sequence',
        help="Choose a reference FASTA file"
    )
    ref_opts.add_argument( # Template sequence
        '-t',
        '--template-sequence',
        dest='template',
        type=str,
        default=None,
        required=True,
        metavar='template sequence',
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
        nargs='+',
        metavar='input file',
        help="Provide one or more FASTQ files, mutually exclusive with '-l | --sample-list' and '-d | --directory'"
    )
    in_files.add_argument( # Are we using a sample list?
        '-l',
        '--sample-list',
        dest='sample_list',
        type=str,
        default=None,
        metavar='sample list',
        help="Provdide a sample list, with each sample on its own line, mutually exclusive with '-i | --input-file' and '-d | --directory'"
    )
    in_files.add_argument( # Are we loading all FASTQ files from a directory?
        '-d',
        '--fastq-directory',
        dest='fastq_directory',
        type=str,
        default=None,
        metavar='FASTQ directory',
        help="Provide a directory with FASTQ files; note that each file must end in .fq or .fastq (optional .gz), mutually exclusive with '-i | --input-file' and '-l | --sample-list'"
    )
    out_opts = parser.add_argument_group( # Output options
        title='output arguments',
        description='Provide an output directory and project name. All files, including the configuration file, will be placed in the output directory with a basename of the project name.'
    )
    out_opts.add_argument( # Output directory
        '-o',
        '--output-directory',
        dest='outdirectory',
        type=str,
        default=_OUTDIR_DEFAULT,
        required=False,
        metavar='output directory',
        help="Choose where all output files are to be stored; defaults to '%s'" % _OUTDIR_DEFAULT
    )
    out_opts.add_argument( # Project name
        '-j',
        '--project',
        dest='project',
        type=str,
        default=_PROJECT_DEFAULT,
        required=False,
        metavar='project',
        help="Provide a project name to be used as the basename for all output files; defaults to '%s'" % _PROJECT_DEFAULT
    )
    out_opts.add_argument( # Output BAM instead of SAM
        '--bam',
        dest='bam',
        type=_validate_bam_index,
        const='bai',
        default=False,
        nargs='?',
        required=False,
        metavar='csi',
        help="Output BAM format instead of SAM, ignored if '--suppress-sam' is passed; for CSI indices, pass '--bam csi', otherwise uses BAI indices"
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
        metavar='read center',
        help="Name of the sequencing center (@RG CN)"
    )
    rg_opts.add_argument( # Read library
        '-rl',
        '--read-library',
        dest='read_library',
        type=str,
        default=None,
        required=False,
        metavar='read library',
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
        metavar='read platform',
        help="Platform technology for sequencing, choose from '%s' (@RG PL)" % "', '".join(_VALID_PLATFORMS)
    )
    rg_opts.add_argument( # Read sample
        '-rs',
        '--read-sample',
        dest='read_sample',
        type=str,
        default=None,
        required=False,
        metavar='read sample',
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
    parser.add_argument( # Enable profiling
        '--profile',
        dest='profile',
        action='store_true',
        required=False,
        help=argparse.SUPPRESS
    )
    parser.add_argument( # Enable XKCD
        '--xkcd',
        dest='xkcd',
        action='store_true',
        required=False,
        help=argparse.SUPPRESS
    )
    return parser
