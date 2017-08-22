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

LOCK = Lock()

try:
    if PYTHON_VERSION is 3:
        from crispronto import plots
        from crispronto import toolkit
        from crispronto import nw_align
        from crispronto import analysis
        from crispronto import alignment
        from crispronto import sam
        from crispronto import quality_control
        from crispronto.arguments import make_argument_parser
        from multiprocessing import pool as Pool
    elif PYTHON_VERSION is 2:
        import plots
        import toolkit
        import nw_align
        import analysis
        import alignment
        import sam
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


def crispr_analysis(
        tup_args # type: Tuple[str, toolkit.NamedSequence, toolkit.NamedSequence, Dict[str, Any], SNP, str]
):
    # type: (...) -> (Tuple[alignment.Alignment], Dict[str, Any])
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
    #   Classify the reads
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
    if not (args_dict['suppress_classification'] or args_dict['suppress_tables']):
        LOCK.acquire()
        hdr_indels, total_counts = analysis.display_classification(
            fastq_name=fastq_name,
            classifications=classifications,
            unique_reads=unique_reads,
            snp_info=snp_info,
            fwd_score=fwd_score,
            rev_score=rev_score,
            score_threshold=score_threshold,
            output_prefix=output_prefix
        )
        LOCK.release()
    #   Count events
    if not (args_dict['suppress_events'] or args_dict['suppress_tables']):
        analysis.events_report(
            fastq_name=fastq_name,
            events=counts,
            reference=reference.sequence,
            snp_index=snp_info.position,
            output_prefix=output_prefix
        )
    #   Make the locus plot
    if not args_dict['suppress_plots']:
        plots.locus_plot(
            insertions=counts['insertions'],
            deletions=counts['deletions'],
            mismatches=counts['mismatches'],
            num_reads=len(reads),
            fastq_name=fastq_name,
            output_prefix=output_prefix
        )
    #   Make a SAM file
    if not args_dict['suppress_sam']:
        logging.info("FASTQ %s: Making SAM file...", fastq_name)
        sam_start = time.time()
        sam_name = os.path.join(output_prefix, fastq_name + '.sam')
        count = 1 # type: int
        sam_lines = list() # type: List[write_sam.SAM]
        bit_base = 16 if do_reverse else 0 # type: int
        reads_dict = defaultdict(list) # type: Mapping[str, List[toolkit.Read]]
        for read in reads:
            reads_dict[read.seq].append(read)
        for aligned in itertools.chain.from_iterable(alignments.values()): # type: alignment.Alignment
            head, tail = toolkit.trim_interval(seq=aligned.read) # type: int, int
            unaligned = toolkit.reverse_complement(sequence=aligned.unaligned) if do_reverse else aligned.unaligned # type: str
            supporting_reads = tuple(reads_dict.pop(unaligned)) # type: Tuple[toolkit.Read]
            position = sam.calc_read_pos(alignment=aligned) # type: int
            cigar = sam.make_cigar(alignment=aligned) # type: str
            sam_seq = sam.make_sam_sequence(alignment=aligned, head=head, tail=tail) # type: str
            read_count = 1 # type: int
            sams = map( # type: Map[write_sam.SAM]
                sam.SAM,
                (read.name for read in supporting_reads), # qname=
                itertools.repeat(bit_base), # flag=
                itertools.repeat(reference.name), # rname=
                itertools.repeat(position), # pos=
                itertools.repeat(255), # mapq=
                itertools.repeat(cigar), # cigar=
                itertools.repeat('*'), # rnex=
                itertools.repeat(0), # pnex=
                itertools.repeat(0), # tlen=
                itertools.repeat(sam_seq), # seq=
                (read.qual for read in supporting_reads) # qual=
            )
            sam_lines.extend(sams)
            count += 1
        if reads_dict:
            unaligned_reads = tuple(itertools.chain.from_iterable(reads_dict.values())) # type: Tuple[toolkit.Read]
            logging.warning("FASTQ %s: %s reads unaligned...", fastq_name, len(unaligned_reads))
        sam_lines = tuple(sorted(sam_lines))
        #   Make the header
        rg_header = sam.make_read_group(sam_lines=sam_lines, conf_dict=args_dict) # type: Tuple[str]
        sq_header = sam.make_sequence_header(
            sam_lines=sam_lines,
            ref_seq_dict={reference.name: reference.sequence}
        )
        #   Write the SAM file
        LOCK.acquire()
        with open(sam_name, 'w') as sfile:
            logging.info("FASTQ %s: Writing %s SAM lines to %s", fastq_name, len(sam_lines), sam_name)
            headers = itertools.chain({sam.HD_HEADER}, sq_header, rg_header) # type: itertools.chain[str]
            for header_line in headers: # type: str
                sfile.write(header_line)
                sfile.write('\n')
                sfile.flush()
            for samline in sam_lines: # type: write_sam.SAM
                sfile.write(str(samline))
                sfile.write('\n')
                sfile.flush()
        LOCK.release()
        logging.debug("FASTQ %s: Making SAM file took %s seconds", fastq_name, round(time.time() - sam_start, 3))
        for samline in sam_lines:
            del samline
    logging.debug("FASTQ %s: Analysis took %s seconds", fastq_name, round(time.time() - analysis_start, 3))
    #   Assemble return values
    if not (args_dict['suppress_classification'] or args_dict['suppress_events'] or args_dict['suppress_tables']):
        mismatch_bases = Counter(itertools.chain.from_iterable(counts['mismatches'].values()))
        total_mismatch = sum(mismatch_bases.values())
        fastq_summary = { # type: Dict[str, Any]
            'filename'          :   fastq_name,
            'total_reads'       :   len(reads),
            'unique_reads'      :   len(unique_reads),
            'discarded'         :   total_counts['DISCARD'],
            'no_edit'           :   total_counts['NO_EDIT'],
            'no_edit_perc'      :   analysis.percent(num=total_counts['NO_EDIT'], total=len(reads)),
            'hdr_clean'         :   total_counts['HDR'] - hdr_indels,
            'hdr_clean_perc'    :   analysis.percent(num=total_counts['HDR'] - hdr_indels, total=total_counts['HDR']),
            'hdr_gap'           :   hdr_indels,
            'hdr_gap_perc'      :   analysis.percent(num=hdr_indels, total=total_counts['HDR']),
            'nhej'              :   total_counts['NHEJ'],
            'nhej_perc'         :   analysis.percent(num=total_counts['NHEJ'], total=len(reads)),
            'perc_a'            :   analysis.percent(num=mismatch_bases['A'], total=total_mismatch),
            'perc_t'            :   analysis.percent(num=mismatch_bases['T'], total=total_mismatch),
            'perc_c'            :   analysis.percent(num=mismatch_bases['C'], total=total_mismatch),
            'perc_g'            :   analysis.percent(num=mismatch_bases['G'], total=total_mismatch)
        }
    else:
        fastq_summary = dict() # type: Dict[str, Any]
    return tuple(toolkit.unpack(collection=alignments.values())), fastq_summary


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
    logging.info("Welcome to %s!", os.path.basename(sys.argv[0]))
    program_start = time.time()
    #   Make an output directory
    output_directory = os.path.join(args['outdirectory'], args['project']) # type: str
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)
    #   Check suppression values
    if _check_suppressions(suppressions=args): # All output suppressed? Error
        sys.exit(logging.critical("All output suppressed, not running"))
    if args['suppress_sam']: # Suppressed SAM output?
        logging.warning("SAM output suppressed, not writing SAM file")
    if args['suppress_events'] or args['suppress_tables']: # Suppressed events table?
        logging.warning("Events output suppressed, not writing events table")
    if args['suppress_classification'] or args['suppress_tables']: # Suppressed classification table?
        logging.warning("Read classification suppressed, not writing classification table")
    if args['suppress_plots']: # Suppressed plots?
        logging.warning("Plots suppressed, not creating plots")
    if args['xkcd']:
        plots._XKCD = True
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
    #   Collect FASTQ information
    if 'sample_list' in args:
        if not os.path.exists(args['sample_list']):
            raise ValueError(logging.critical("Cannot find sample list %s", args['sample_list']))
        with open(args['sample_list'], 'r') as listfile:
            fastq_list = tuple(line.strip() for line in listfile if not line.startswith('#')) # type: Tuple[str]
    elif 'input_file' in args:
        fastq_list = (args['input_file'],) # type: Tuple[str]
    else:
        sys.exit(logging.critical("No input file provided"))
    zipped_args = zip(
        fastq_list,
        itertools.repeat(reference),
        itertools.repeat(aligned_reference),
        itertools.repeat(args),
        itertools.repeat(snp),
        itertools.repeat(output_directory)
    )
    #   Setup our multiprocessing pool
    #   Allow the user to specify the number of jobs to run at once
    #   If not specified, let multiprocessing figure it out
    try:
        pool = Pool(processes=args['num_cores']) # type: multiprocessing.Pool
    except KeyError:
        pool = Pool() # type: multiprocessing.Pool
    #   If we have multiple FASTQ files AND multiple processes running
    #   use pool.map_async; else use generic map to avoid timeout issues
    if all(map(lambda i: i > 1, (len(fastq_list), getattr(pool, '_processes')))):
        try:
            #   Use map_async and get with a large timeout
            #   to allow for KeyboardInterrupts to be caught
            #   and handled with the try/except
            res = pool.map_async(crispr_analysis, zipped_args) # type: multiprocessing.pool.MapResult
            pool.close()
            results = res.get(9999)
        except KeyboardInterrupt: # Handle ctrl+c
            pool.terminate()
            sys.exit('\nkilled')
        except SystemExit: # Handle sys.exit() calls
            pool.terminate()
            raise
    #   Otherwise, don't bother with pool.map() make life easy
    else:
        #   Clean up the pool
        pool.close(); pool.terminate()
        #   Use standard map (or itertools.imap)
        results = map(crispr_analysis, zipped_args)
    alignments, summaries = zip(*results)
    alignments = toolkit.unpack(collection=alignments)
    output_prefix = os.path.join(output_directory, args['project'])
    if not args['suppress_plots']:
        plots.quality_plot(
            alignments=alignments,
            output_prefix=output_prefix
        )
    if not (args['suppress_classification'] or args['suppress_events'] or args['suppress_tables']):
        summary_name = output_prefix + '.summary'
        summary_header = (
            '#FASTQ',
            'TOTAL_READS',
            'UNIQ_READS',
            'DISCARDED',
            'SNP_POS',
            'REF_STATE',
            'TEMP_SNP',
            'NO_EDIT',
            'PERC_NO_EDIT',
            'HDR_CLEAN',
            'PERC_HDR_CLEAN',
            'HDR_GAP',
            'PERC_HDR_GAP',
            'NHEJ',
            'NHEJ_GAP',
            'PERC_MIS_A',
            'PERC_MIS_T',
            'PERC_MIS_C',
            'PERC_MIS_G'
        )
        logging.info("Writing summary to %s", summary_name)
        summary_start = time.time()
        with open(summary_name, 'w') as summfile:
            summfile.write('\t'.join(summary_header) + '\n')
            summfile.flush()
            for sum_dict in sorted(summaries, key=lambda d: d['filename']):
                out = (
                    sum_dict['filename'],
                    sum_dict['total_reads'],
                    sum_dict['unique_reads'],
                    sum_dict['discarded'],
                    snp.position + 1,
                    snp.reference,
                    snp.target,
                    sum_dict['no_edit'],
                    sum_dict['no_edit_perc'],
                    sum_dict['hdr_clean'],
                    sum_dict['hdr_clean_perc'],
                    sum_dict['hdr_gap'],
                    sum_dict['hdr_gap_perc'],
                    sum_dict['nhej'],
                    sum_dict['nhej_perc'],
                    sum_dict['perc_a'],
                    sum_dict['perc_t'],
                    sum_dict['perc_c'],
                    sum_dict['perc_g']
                )
                out = map(str, out)
                summfile.write('\t'.join(out))
                summfile.write('\n')
                summfile.flush()
        logging.debug("Writing summary took %s seconds", round(time.time() - summary_start, 3))
    #   Close logfile
    logging.debug("Entire program took %s seconds to run", round(time.time() - program_start, 3))
    args['logfile'].close()


if __name__ == '__main__':
    main()
