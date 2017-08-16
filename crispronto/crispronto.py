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
from collections import Counter, defaultdict, namedtuple

try:
    if PYTHON_VERSION is 3:
        from crispronto import plots
        from crispronto import toolkit
        from crispronto import nw_align
        from crispronto import analysis
        from crispronto import alignment
        from crispronto import quality_control
        from crispronto.arguments import make_argument_parser
        from multiprocessing import pool as Pool
    elif PYTHON_VERSION is 2:
        import plots
        import toolkit
        import nw_align
        import analysis
        import alignment
        import quality_control
        from arguments import make_argument_parser
        from multiprocessing import Pool
        from itertools import izip as zip
        from itertools import imap as map
        range = xrange
    else:
        sys.exit('Please use Python 2.7, 3.3, or higher')
except ImportError as error:
    sys.exit(error)


SNP = namedtuple('SNP', ('reference', 'target', 'position'))
Events = namedtuple('Events', ('num_reads', 'num_ins', 'num_del', 'num_mis'))

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
    fastq_file, reference, aligned_reference, args_dict, snp_info, output_directory = tup_args # type: str, toolkit.NamedSequence, toolkit.NamedSequence, Dict[str, Any], SNP, str
    #   Make a name for our FASTQ file
    fastq_name = os.path.basename(fastq_file) # type: str
    if fastq_name.count('.') == 2:
        fastq_name = fastq_name.split('.')[0]
    else:
        fastq_name = os.path.splitext(fastq_name)[0]
    output_prefix = os.path.join(output_directory, fastq_name)
    if not os.path.exists(output_prefix):
        os.makedirs(output_prefix)
    logging.info("FASTQ %s: Starting analysis...", fastq_name)
    analysis_start = time.time() # type: float
    #   Load the reads
    reads = toolkit.load_fastq(fastq_file=fastq_file) # type: Tuple[toolpack.Read]
    #   Determine the alignment direction for our reads
    unique_reads = Counter(map(lambda read: read.seq, reads)) # type: Counter
    do_reverse, fwd_score, rev_score, score_threshold = quality_control.determine_alignment_direction( # type: bool, float, float, float
        fastq_name=fastq_name,
        unique_reads=unique_reads.keys(),
        reference=aligned_reference.sequence,
        gap_open=args_dict['gap_opening'],
        gap_extension=args_dict['gap_extension'],
        pvalue_threshold=args_dict['pvalue_threshold']
    )
    if do_reverse:
        unique_reads = {toolkit.reverse_complement(sequence=read): count for read, count in unique_reads.items()} # type: Dict[str, count]
    #   Sort the reads by length
    reads_by_length = alignment.sort_reads_by_length(reads=unique_reads.keys(), fastq_name=fastq_name) # type: Dict[int, List[str]]
    alignments = alignment.align_recurse( # type: Dict[int, List[alignment.Alignment]]
        fastq_name=fastq_name,
        reads_by_length=reads_by_length,
        reference=aligned_reference.sequence,
        gap_open=args_dict['gap_opening'],
        gap_extension=args_dict['gap_extension']
    )
    logging.info("FASTQ %s: Starting read classifcation...", fastq_name)
    classifcation_start = time.time() # type: float
    counts = { # type: Dict[str, defaultdict]
        'deletions': defaultdict(list),
        'insertions': defaultdict(list),
        'mismatches': defaultdict(list),
        'matches': defaultdict(int)
    }
    classifications = (dict(), dict(), dict(), dict()) # type: Tuple[Dict[str, Events]]
    hdr, nhej, no_edit, discard = classifications # type: Dict[str, Events]
    for length, alignment_list in alignments.items(): # type: int, List[alignment.Alignment]
        for aligned in alignment_list: # type: alignment.Alignment
            num_ins, num_del, num_mis = 0, 0, 0 # type: int, int, int
            num_reads = unique_reads[str(aligned)] # type: int
            if aligned.score < score_threshold:
                discard[str(aligned)] = Events(num_reads=num_reads, num_ins=0, num_del=0, num_mis=0)
                continue
            ref_head, ref_tail = toolkit.trim_interval(seq=aligned.reference) # type: int, int
            al_ref_seq = aligned.reference[ref_head:ref_tail] # type: str
            al_read_seq = aligned.read[ref_head:ref_tail] # type: str
            read_head, read_tail = toolkit.trim_interval(seq=al_read_seq) # type: int, int
            insertion_list = toolkit.find_gaps(seq=al_ref_seq) # type: List
            if insertion_list:
                num_ins = len(insertion_list)
                temp = 0
                for gap in insertion_list:
                    position, gi_length = gap
                    position = position - temp
                    temp = temp + gi_length
                    counts['insertions'][position].extend([gi_length] * num_reads)
            if num_ins:
                ref_no_ins, read_no_ins = str(), str() # type: str, str
                for index in range(len(al_ref_seq)): # type: int
                    if al_ref_seq[index] == '-':
                        continue
                    else:
                        ref_no_ins = ref_no_ins + al_ref_seq[index] # type: str
                        read_no_ins = read_no_ins + al_read_seq[index] # type: str
                al_ref_seq = ref_no_ins # type: str
                al_read_seq = read_no_ins # type: str
            deletion_list = toolkit.find_gaps(seq=al_read_seq, head=read_head, tail=read_tail) # type: List
            if deletion_list:
                num_del = len(deletion_list)
                for gap in deletion_list:
                    position, gd_length = gap
                    counts['deletions'][position].extend([gd_length] * num_reads)
            mismatches, matches = toolkit.get_mismatch(
                seq_a=al_ref_seq,
                seq_b=al_read_seq,
                head=read_head,
                tail=read_tail,
                matches=True
            )
            for mis in mismatches:
                position, perm = mis
                source, snp = perm
                counts['mismatches'][position].extend([snp] * num_reads)
            for match_position in matches:
                counts['matches'][match_position] += num_reads
            results = Events(num_reads=num_reads, num_ins=num_ins, num_del=num_del, num_mis=num_mis) # type: Events
            if al_read_seq[snp_info.position] == snp_info.target:
                hdr[str(aligned)] = results
            elif all(map(lambda x: not x, (num_del, num_ins, num_mis))):
                no_edit[str(aligned)] = results
            else:
                nhej[str(aligned)] = results
    logging.debug("FASTQ %s: Read classifcation took %s seconds", fastq_name, round(time.time() - classifcation_start, 3))
    cummulative_deletions = analysis.events_report(
        fastq_name=fastq_name,
        events=counts,
        reference=reference.sequence,
        snp_index=snp_info.position,
        output_prefix=output_prefix
    )
    analysis.display_classification(
        fastq_name=fastq_name,
        classifications=classifications,
        unique_reads=unique_reads,
        snp_info=snp_info,
        fwd_score=fwd_score,
        rev_score=rev_score,
        score_threshold=score_threshold,
        output_prefix=output_prefix
    )
    plots.locus_plot(
        insertions=counts['insertions'],
        deletions=counts['deletions'],
        # deletions=cummulative_deletions,
        mismatches=counts['mismatches'],
        num_reads=len(reads),
        # num_reads=len(unique_reads),
        fastq_name=fastq_name,
        output_prefix=output_prefix
    )
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
    #   Make an output directory
    output_directory = os.path.join(args['outdirectory'], args['project']) # type: str
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)
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
    sigint_handler = signal.signal(signal.SIGINT, signal.SIG_IGN) # type: function
    signal.signal(signal.SIGINT, sigint_handler)
    #   Read in reference and template sequences
    logging.info("Quality control...")
    qc_start = time.time() # type: float
    reference = toolkit.load_seq(args['reference']) # type: toolkit.NamedSequence
    template = toolkit.load_seq(args['template']) # type: toolkit.NamedSequence
    #   Align template and reference sequences to determine alignment direction
    al_ref_seq, al_temp_seq = quality_control.align_reference( # type: str, str
        reference=reference.sequence,
        template=template.sequence,
        gap_penalty=args['gap_opening']
    )
    aligned_reference = toolkit.NamedSequence(name=reference.name, sequence=al_ref_seq) # type: toolkit.NamedSequence
    aligned_template = toolkit.NamedSequence(name=template.name, sequence=al_temp_seq) # type: toolkit.NamedSequence
    #   QC the alignments
    logging.info("Validating reference/template alignment...")
    alignment_validation = time.time() # type: float
    if '-' in set(aligned_reference.sequence):
        raise ValueError(logging.error("Cannot have insertions in the reference"))
    if '-' in set(toolkit.side_trimmer(seq=aligned_template.sequence)):
        raise ValueError(logging.error("Cannot have deletions in the template sequence"))
    template_reference_mismatch = toolkit.get_mismatch(seq_a=aligned_reference.sequence, seq_b=aligned_template.sequence) # type: List
    if not template_reference_mismatch:
        raise ValueError(logging.error("There must be at least one mismatch between the template and reference sequences"))
    if len(template_reference_mismatch) > 1 and args['analysis_mode'] == 'SNP':
        raise ValueError(logging.error("There can only be one mismatch between the template and reference sequences in SNP mode"))
    logging.debug("Reference/template aligmnent validation took %s seconds", round(time.time() - alignment_validation, 3))
    #   Get SNP information
    snp_index, reference_state, target_snp = quality_control.get_snp_states( # type: int, str, str
        reference=aligned_reference.sequence,
        template=aligned_template.sequence,
        mismatch=template_reference_mismatch,
        mode=args['analysis_mode']
    )
    snp = SNP(reference=reference_state, target=target_snp, position=snp_index) # type: SNP
    logging.debug("Quality control took %s seconds", round(time.time() - qc_start, 3))
    crispr_analysis(tup_args=(args['input_file'], reference, aligned_reference, args, snp, output_directory))
    #   Setup our multiprocessing pool
    #   Allow the user to specify the number of jobs to run at once
    #   If not specified, let multiprocessing figure it out
    # try:
    #     pool = Pool(processes=args['num_cores']) # type: multiprocessing.Pool
    # except KeyError:
    #     pool = Pool() # type: multiprocessing.Pool
    # if all(map(lambda i: i > 1, (len(fastq_list), getattr(pool, '_processes')))):
    #     pass
    #   Close logfile
    args['logfile'].close()


if __name__ == '__main__':
    main()
