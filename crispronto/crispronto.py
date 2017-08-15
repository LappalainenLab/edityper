#!/usr/bin/env python

"""CRISPRonto"""

from __future__ import print_function
from __future__ import division

import sys
PYTHON_VERSION = sys.version_info.major

import os
import time
import logging
import signal
import itertools
from multiprocessing import Lock
from collections import Counter, defaultdict

try:
    if PYTHON_VERSION is 3:
        from crispronto import toolkit
        from crispronto import nw_align
        from crispronto import alignment
        from crispronto import quality_control
        from crispronto.arguments import make_argument_parser
        from multiprocessing import pool as Pool
    elif PYTHON_VERSION is 2:
        import toolkit
        import nw_align
        import alignment
        import quality_control
        from arguments import make_argument_parser
        from multiprocessing import Pool
        from itertools import izip as zip
        from itertools import imap as map
        range = xrange
    else:
        sys.exit('WTF MATE')
except ImportError as error:
    sys.exit("FAIL")


def _set_verbosity(level): # type: (str) -> int
    level = level.upper()
    if level == 'DEBUG':
        log_level = logging.DEBUG # type: int
    elif level == 'INFO':
        log_level = logging.INFO # type: int
    elif level == 'WARNING':
        log_level = logging.WARNING # type: int
    elif level == 'ERROR':
        log_level = logging.ERROR # type: int
    elif level == 'CRITICAL':
        log_level = logging.CRITICAL # type: int
    else:
        raise ValueError("'level' must be one of 'DEBUG', 'INFO', 'WARNING', 'ERROR', or 'CRITICAL'")
    return log_level


def _check_suppressions(suppressions): # type: (Dict[str, bool]) -> bool
    suppress_tables = suppressions['suppress_events'] and suppressions['suppress_classification']
    return suppressions['suppress_plots'] and suppressions['suppress_sam'] and (suppressions['suppress_tables'] or suppress_tables)


def crispr_analysis(tup_args):
    """Do the CRISPR analysis"""
    fastq_file, aligned_reference, args_dict = tup_args
    #   Make a name for our FASTQ file
    fastq_name = os.path.basename(fastq_file) # type: str
    if fastq_name.count('.') == 2:
        fastq_name = fastq_name.split('.')[0]
    else:
        fastq_name = os.path.splitext(fastq_name)[0]
    logging.info("FASTQ %s: Starting analysis", fastq_name)
    analysis_start = time.time()
    #   Load the reads
    reads = toolkit.load_fastq(fastq_file=fastq_file) # type: Tuple[toolpack.Read]
    #   Determine the alignment direction for our reads
    unique_reads = Counter(map(lambda read: read.seq, reads)) # type: Counter
    do_reverse, score_threshold = quality_control.determine_alignment_direction(
        fastq_name=fastq_name,
        unique_reads=unique_reads.keys(),
        reference=aligned_reference.sequence,
        gap_open=args_dict['gap_opening'],
        gap_extension=args_dict['gap_extension'],
        pvalue_threshold=args_dict['pvalue_threshold']
    )
    if do_reverse:
        unique_reads = {toolkit.reverse_complement(sequence=read): count for read, count in unique_reads.items()}
    #   Sort the reads by length
    reads_by_length = alignment.sort_reads_by_length(reads=unique_reads.keys(), fastq_name=fastq_name)
    alignments = alignment.align_recurse( # type: Dict[int, List[alignment.Alignment]]
        fastq_name=fastq_name,
        reads_by_length=reads_by_length,
        reference=aligned_reference.sequence,
        gap_open=args_dict['gap_opening'],
        gap_extension=args_dict['gap_extension']
    )
    logging.info("Collecting alignment scores")
    logging.debug("FASTQ %s: Analysis took %s seconds", fastq_name, round(time.time() - analysis_start, 3))


def main():
    """CRISPronto"""
    #   Setup CRISPRonto
    #   Parse arguments
    parser = make_argument_parser() # type: argparse.ArgumentParser
    if not sys.argv[1:]:
        sys.exit(parser.print_help())
    args = {key: value for key, value in vars(parser.parse_args()).items() if value is not None} # type: Dict[str, Any]
    #   Setup logger
    logging.basicConfig(
        format='%(asctime)s %(levelname)s:\t%(message)s',
        stream=args['logfile'],
        level=_set_verbosity(level=args['verbosity']),
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    #   Check suppression values
    suppressions = {key: value for key, value in args.items() if 'suppress' in key} # type: Dict[str, bool]
    if _check_suppressions(suppressions=suppressions): # All output suppressed? Error
        sys.exit(logging.critical("All output suppressed, not running"))
    if suppressions['suppress_sam']: # Suppressed SAM output?
        logging.warning("SAM output suppressed, not writing SAM file")
    if suppressions['suppress_events'] or suppressions['suppress_tables']: # Suppressed events table?
        logging.warning("Events output suppressed, not writing events table")
    if suppressions['suppress_classification'] or suppressions['suppress_tables']: # Suppressed classification table?
        logging.warning("Read classification suppressed, not writing classification table")
    if suppressions['suppress_plots']: # Suppressed plots?
        logging.warning("Plots suppressed, not creating plots")
    if args['xkcd']:
        # plots._XKCD = True
        pass
    #   Tell worker processes to ignore SIGINT (^C)
    #   by turning INTERUPT signals into IGNORED signals
    sigint_handler = signal.signal(signal.SIGINT, signal.SIG_IGN)
    signal.signal(signal.SIGINT, sigint_handler)
    #   Read in reference and template sequences
    logging.info("Quality control")
    qc_start = time.time()
    reference = toolkit.load_seq(args['reference'])
    template = toolkit.load_seq(args['template'])
    #   Align template and reference sequences to determine alignment direction
    al_ref_seq, al_temp_seq = quality_control.align_reference(
        reference=reference.sequence,
        template=template.sequence,
        gap_penalty=args['gap_opening']
    )
    aligned_reference = toolkit.NamedSequence(name=reference.name, sequence=al_ref_seq)
    aligned_template = toolkit.NamedSequence(name=template.name, sequence=al_temp_seq)
    #   QC the alignments
    logging.info("Validating reference/template alignment")
    alignment_validation = time.time()
    if '-' in set(aligned_reference.sequence):
        raise ValueError(logging.error("Cannot have insertions in the reference"))
    if '-' in set(toolkit.side_trimmer(seq=aligned_template.sequence)):
        raise ValueError(logging.error("Cannot have deletions in the template sequence"))
    template_reference_mismatch = toolkit.get_mismatch(seq_a=aligned_reference.sequence, seq_b=aligned_template.sequence)
    if not template_reference_mismatch:
        raise ValueError(logging.error("There must be at least one mismatch between the template and reference sequences"))
    if len(template_reference_mismatch) > 1 and args['analysis_mode'] == 'SNP':
        raise ValueError(logging.error("There can only be one mismatch between the template and reference sequences in SNP mode"))
    logging.debug("Reference/template aligmnent validation took %s seconds", round(time.time() - alignment_validation, 3))
    #   Get SNP information
    snp_index, reference_state, target_snp = quality_control.get_snp_states(
        reference=aligned_reference.sequence,
        template=aligned_template.sequence,
        mismatch=template_reference_mismatch,
        mode=args['analysis_mode']
    )
    logging.debug("Quality control took %s seconds", round(time.time() - qc_start, 3))
    crispr_analysis(tup_args=(args['input_file'], aligned_reference, args))
    #   Setup our multiprocessing pool
    #   Allow the user to specify the number of jobs to run at once
    #   If not specified, let multiprocessing figure it out
    try:
        pool = Pool(processes=args['num_cores']) # type: multiprocessing.Pool
    except KeyError:
        pool = Pool() # type: multiprocessing.Pool
    # if all(map(lambda i: i > 1, (len(fastq_list), getattr(pool, '_processes')))):
    #     pass
    #   Close logfile
    args['logfile'].close()


if __name__ == '__main__':
    main()
