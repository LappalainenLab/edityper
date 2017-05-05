#!/usr/bin/env python

"""Align and analyze CRISPR data"""

from __future__ import division
from __future__ import print_function

import sys
if sys.version_info.major is not 2 and sys.version_info.minor is not 7:
    sys.exit("Please use Python 2.7 for this program")


import os
import time
import logging
import itertools
from copy import deepcopy
from collections import namedtuple
from multiprocessing import freeze_support, Pool, Lock

try:
    import plots
    import arguments
    import configure
    import analysis as an
    import alignment as al
    import write_sam as ws
    import read_isolation as ri
    import quality_control as qc
    import genetic_toolpack as toolpack
except ImportError as error:
    sys.exit("Please keep this program in it's directory to load custom modules: " + error.message)


LOCK = Lock()
NamedSequence = namedtuple('NamedSequence', ('name', 'sequence'))

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


def load_data(conf_dict): # type: (Dict[Any]) -> (str, str, str, str, List[toolpack.FastQ]
    """Load our data"""
    logging.info("Loading data...")
    load_start = time.time()
    if 'sample_list' in conf_dict:
        with open(conf_dict['sample_list'], 'r') as listfile:
            fastq_files = [line.strip() for line in listfile] # type: List[str]
    elif 'input_file' in conf_dict:
        fastq_files = [conf_dict['input_file']] # type: List[str]
    else:
        raise ValueError("Cannot find the input FASTQ file(s)")
    logging.info("Loading %s as reference", conf_dict['reference'])
    ref_name, ref_seq = toolpack.load_seq(conf_dict['reference']) # type: str, str
    logging.info("Loading %s as template", conf_dict['template'])
    temp_name, temp_seq = toolpack.load_seq(conf_dict['template']) # type: str, str
    fastq_list = [toolpack.load_fastq(f) for f in fastq_files] # type: List[toolpack.FastQ]
    logging.debug("Data load took %s seconds", round(time.time() - load_start, 3))
    return ref_name, ref_seq, temp_name, temp_seq, fastq_list


def run_qc(
        conf_dict, # type: Dict[str, Any]
        reference, # type: str
        template, # type: str
        fastq, # type: toolpack.FastQ
):
    # type: (...) -> (str, Dict[str, List[toolpack.Read]], numpy.float64, int, str, bool)
    """Run the QC steps"""
    logging.info("Running quality control")
    qc_start = time.time()
    #   Validate our reference and template sequences by aligning them
    aligned_reference, aligned_template = qc.align_reference( # type: str, str
        reference=reference,
        template=template,
        gap_penalty=conf_dict['gap_opening']
    )
    reference_template_mismatch = qc.validate_reference_alignment( # type: List[List[int, Tuple[str]]]
        reference=aligned_reference,
        template=aligned_template,
        snp_mode=True if conf_dict['analysis_mode'] == 'SNP' else False
    )
    #   Get the reference SNP state and target SNP from our aligned reference and template
    snp_index, reference_state, target_snp = qc.get_snp_states( # type: int, str, str
        reference=aligned_reference,
        template=aligned_template,
        mismatch=reference_template_mismatch,
        mode=conf_dict['analysis_mode']
    )
    #   Determine alignment direction
    #   Need chain to flatten the list of lists
    raw_read_dict = ri.load_seqs( # type: Dict[str, List[toolpack.Read]]
        raw_reads=fastq.get_all_reads()
    ) # raw_read_dict is {unique sequences: [reads that have this sequence]}
    do_reverse, fwd_median, rev_median, score_threshold = qc.determine_alignment_direction( # type: bool, numpy.float64, numpy.float64, numpy.float64
        raw_sequences=raw_read_dict, # Effectively keys (unique sequences) only
        reference=aligned_reference,
        gap_open=conf_dict['gap_opening'],
        gap_extension=conf_dict['gap_extension'],
        pvalue_threshold=conf_dict['pvalue_threshold']
    )
    if do_reverse:
        reads_dict = ri.reverse_reads(reads_dict=raw_read_dict) # type: Dict[str, List[toolpack.Read]]
    else:
        reads_dict = deepcopy(raw_read_dict) # type: Dict[str, List[toolpack.Read]]
    logging.debug("Quality control took %s seconds", round(time.time() - qc_start, 3))
    return reads_dict, fwd_median, rev_median, score_threshold, snp_index, reference_state, target_snp, do_reverse


def run_alignment(
        reference, # type: str
        reads_dict, # type: Dict[str, int]
        conf_dict # type: Dict[str, Any]
):
    # type: (...) -> Dict[int, List[al.Alignment]]
    """Run the alignment steps"""
    logging.info("Aliging reads")
    alignment_start = time.time()
    sorted_reads = al.sort_reads_by_length(reads_dict=reads_dict) # type: Dict[int, List[al.ReadSummary]]
    alignments = al.align_recurse( # type: Dict[int, List[al.Alignment]]
        reads_by_length=sorted_reads,
        reference=reference,
        gap_open=conf_dict['gap_opening'],
        gap_ext=conf_dict['gap_extension']
    )
    logging.debug("Alignment took %s seconds", round(time.time() - alignment_start, 3))
    # return scores, aligned, sorted_reads
    return alignments


def make_sam_file(
        conf_dict, # type: Dict[str, Any]
        reads_dict, # type: Dict[str, List[toolpack.Read]]
        reference_name, # type: str
        reference_seq, # type: str
        alignments, # type: Iterable[al.Alignment]]
        is_reverse, # type: bool
        output_prefix # type: str
):
    # type: (...) -> None
    """Make a SAM file"""
    logging.info("Making SAM file")
    sam_start = time.time() # type: float
    bit_base = 0 # type: int
    if is_reverse:
        bit_base = 16 # type: int
    #   Get CIGAR strings
    logging.info("Creating CIGAR strings")
    cigar_start = time.time() # type: float
    cigars = map(ws.make_cigar, alignments) # type: List[str]
    logging.debug("Making CIGAR strings took %s seconds", round(time.time() - cigar_start, 3))
    #   Get read positions
    logging.info("Finding read positions on aligned reference")
    positions_start = time.time() # type: float
    positions = map(ws.calc_read_pos, alignments) # type: List[int]
    logging.debug("Finding read positions took %s seconds", round(time.time() - positions_start, 3))
    #   Convert the aligned reads to SAM-ready reads
    #   Get the heads and tails grouped by alignment
    #   Then, unzip the heads and tails into separate tuples
    logging.info("Creating SAM-ready sequences")
    ready_start = time.time()
    head_tail = map(toolpack.trim_interval, (aligned.get_aligned_read() for aligned in alignments)) # type: List[Tuple[Int]]
    heads, tails = zip(*head_tail) # type: Tuple[int], Tuple[int]
    sam_seq = map(ws.make_sam_sequence, alignments, heads, tails) # type: List[str]
    logging.debug("Creating SAM-ready sequences took %s seconds", round(time.time() - ready_start, 3))
    #   Get the number of times we're repeating our constants
    #   Apparently, Python 2 doesn't stop automatically...
    num_repeats = set(map(len, (alignments, sam_seq, cigars, positions))) # type: Set[int]
    if len(num_repeats) != 1:
        sys.exit(logging.critical("Unequal numbers of alignments, sequences, CIGAR strings, and positions; exiting..."))
    #   Get the number of times from the set
    num_repeats = num_repeats.pop() # type: int
    #   Make the alignment lines of the SAM file
    logging.info("Creating SAM lines")
    lines_start = time.time()
    sam_lines = map(
        ws.make_sam, # Func
        alignments, # alignment=
        sam_seq, # sam_seq=
        cigars, # cigar=
        itertools.repeat(bit_base, num_repeats), # bit_flag=
        itertools.repeat(reference_name, num_repeats), # ref_name=
        positions, # position=
        itertools.repeat(reads_dict, num_repeats) # reads_dict=
    )
    logging.debug("Making SAM lines took %s seconds", round(time.time() - lines_start, 3))
    #   Unpack and sort our SAM lines
    #   The sort works because SAM objects have __le__ and __lt__ methods defined
    #   These definitions allow for coordinate sort based on rname and pos fields
    #   as specified in the SAMv1 specification (@HD SO definition on page 3)
    sam_lines = itertools.chain.from_iterable(sam_lines) # type: itertools.chain
    sam_lines = sorted(sam_lines) # type: List[ws.SAM]
    #   Make header lines
    logging.info("Making headers")
    header_start = time.time() # type: float
    rg_header = ws.make_read_group(sam_lines=sam_lines, conf_dict=conf_dict) # type: Tuple[str]
    sq_header = ws.make_sequence_header( # type: Tuple[str]
        sam_lines=sam_lines,
        ref_seq_dict={reference_name: reference_seq}
    )
    logging.debug("Making headers took %s seconds", round(time.time() - header_start, 3))
    #   Write the SAM file
    sam_name = output_prefix + '.sam' # type: str
    logging.info("Writing SAM information to %s", sam_name)
    write_start = time.time()
    with open(sam_name, 'w') as samfile:
        #   Write the all the header lines
        #   Chain the headers together so only one forloop for all header lines is needed
        #   itertools.chain keeps the order
        headers = itertools.chain({ws.HD_HEADER}, sq_header, rg_header)
        for header_line in headers:
            samfile.write(header_line)
            samfile.write('\n')
            samfile.flush()
        for sam in sam_lines:
            samfile.write(str(sam))
            samfile.write('\n')
            samfile.flush()
    logging.debug("Writing SAM information to file took %s seconds", round(time.time() - write_start, 3))
    logging.debug("Making SAM file took %s seconds", round(time.time() - sam_start, 3))


def crispr(tup_args): # type: Tuple(...) -> None
    """Run the alignment and analysis on a single FASTQ file"""
    LOCK.acquire()
    fastq, suppression, conf_dict, seq_dict = tup_args # type: toolpack.FastQ, Dict[str, bool], Dict[str, Any], Dict[str, Sequence]
    logging.info("Starting alignment and analysis of %s", str(fastq))
    fastq_start = time.time()
    if str(fastq).count('.') <= 2:
        fastq_base = str(fastq).split('.')[0]
    else:
        fastq_base = os.path.splitext(str(fastq))[0]
    outdir = conf_dict['outdirectory'][:-1] if conf_dict['outdirectory'][-1] == '/' else conf_dict['outdirectory']
    outdir = '/'.join((outdir, fastq_base))
    configure.mkdir(outdir)
    output_prefix = configure.make_prefix_name( # type: str
        directory=outdir,
        base='_'.join((fastq_base, conf_dict['project']))
    )
    #   Run through quality control steps
    reads_dict, fwd_median, rev_median, score_threshold, snp_index, reference_state, target_snp, do_reverse = run_qc( # type: Dict[str, List[toolpack.Read]], numpy.float64, numpy.float64, numpy.float64, int, str, str, bool
        conf_dict=conf_dict,
        reference=seq_dict['reference'].sequence,
        template=seq_dict['template'].sequence,
        fastq=fastq
    )
    total_reads = sum((len(read_list) for read_list in reads_dict.values())) # type: int
    #   Align the reads
    alignments = run_alignment( # type: Dict[int, List[al.Alignment]]
        reference=seq_dict['reference'].sequence,
        reads_dict=reads_dict,
        conf_dict=conf_dict
    )
    #   Run analyses
    report, read_classifications = an.run_analysis( # type: an.Reporter, Tuple[defaultdict[int, List[al.Alignment]]]
        reads_dict=reads_dict,
        alignments=itertools.chain.from_iterable(alignments.values()),
        score_threshold=score_threshold,
        snp_index=snp_index,
        target_snp=target_snp
    )
    #   Start outputs
    if suppression['suppress_classification'] or suppression['suppress_tables']:
        logging.warning("Read classification suppressed, not writing classification table")
    else:
        an.display_classification(
            read_classification=read_classifications,
            total_reads=total_reads,
            snp_position=snp_index,
            ref_state=reference_state,
            target_snp=target_snp,
            fwd_score=fwd_median,
            rev_score=rev_median,
            score_threshold=score_threshold,
            output_prefix=output_prefix
        )
    if suppression['suppress_sam']:
        logging.warning("SAM output suppressed, not writing SAM file")
    else:
        make_sam_file(
            conf_dict=conf_dict,
            reads_dict=reads_dict,
            reference_name=seq_dict['reference'].name,
            reference_seq=seq_dict['reference'].sequence,
            alignments=tuple(al for al in itertools.chain.from_iterable(alignments.values())),
            is_reverse=do_reverse,
            output_prefix=output_prefix
        )
    if suppression['suppress_events'] or suppression['suppress_tables']:
        logging.warning("Events output suppressed, not writing events table")
    else:
        an.create_report(
            reporter=report,
            reference=seq_dict['reference'].sequence,
            snp_index=snp_index,
            output_prefix=output_prefix
        )
    #   Plotting
    if suppression['suppress_plots']:
        logging.warning("Plots suppressed, not creating plots")
    else:
        plots.locus_plot(
            insertions=report.insertions,
            deletions=report.deletions,
            mismatches=report.mismatches,
            num_reads=total_reads,
            fastq_name=str(fastq),
            output_prefix=output_prefix
        )
    logging.debug("Alignment and analysis of %s took %s seconds", str(fastq), round(time.time() - fastq_start, 3))
    LOCK.release()


def main(args): # type: (Dict) -> None
    """Run the program"""
    try:
        if args['subroutine'] == 'CONFIG':
            configure.write_config(args)
        elif args['subroutine'] == 'ALIGN':
            suppressions = {key: value for key, value in args.items() if 'suppress' in key} # type: Dict[str, bool]
            conf_dict = configure.read_config(args['config_file']) # type: Dict
            output_prefix = configure.make_prefix_name(directory=conf_dict['outdirectory'], base=conf_dict['project'])
            #   Load the data
            ref_name, ref_seq, temp_name, temp_seq, fastq_list = load_data( # type: str, str str, str, List[toolpack.FastQ]]
                conf_dict=conf_dict
            )
            sequences = { # type: Dict[str, NamedSequence]
                'reference': NamedSequence(name=ref_name, sequence=ref_seq),
                'template': NamedSequence(name=temp_name, sequence=temp_seq)
            }
            if args['xkcd']:
                plots._XKCD = True
            pool = Pool() # type: multiprocessing.Pool
            pool.map(
                crispr,
                itertools.izip(
                    fastq_list,
                    itertools.repeat(suppressions),
                    itertools.repeat(conf_dict),
                    itertools.repeat(sequences)
                )
            )
            #     plots.quality_plot(
            #         alignments=tuple(al for al in itertools.chain.from_iterable(alignments.values())),
            #         output_prefix=output_prefix
            #     )
    except IOError as error:
        sys.exit(logging.critical("Cannot find file %s, exiting", error.filename))
    except:
        raise


if __name__ == '__main__':
    PARSER = arguments.make_argument_parser() # type: argparse.ArgumentParser
    if not sys.argv[1:]:
        sys.exit(PARSER.print_help())
    #   Collect our arguments
    ARGS = {key: value for key, value in vars(PARSER.parse_args()).items() if value is not None} # type: Dict
    #   Set logging parameters
    logging.basicConfig(
        format='%(asctime)s %(levelname)s:\t%(message)s',
        stream=ARGS['logfile'],
        level=_set_verbosity(level=ARGS['verbosity']),
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    #   Run the program
    main(args=ARGS)
    #   Close our log file
    ARGS['logfile'].close()
else:
    logging.basicConfig(
        format='%(asctime)s %(levelname)s:\t%(message)s',
        level=logging.INFO,
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    logging.info("Logger level set to 'INFO'")
