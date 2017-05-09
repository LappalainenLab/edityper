#!/usr/bin/env python2

"""Analysis of alignments for the CRISPR program"""

from __future__ import division
from __future__ import print_function

import sys
if sys.version_info.major is not 2 and sys.version_info.minor is not 7:
    sys.exit("Please use Python 2.7 for this module: " + __name__)


import re
import time
import logging
import warnings
import itertools
from collections import Counter, defaultdict, namedtuple

try:
    import genetic_toolpack as toolpack
except ImportError as error:
    sys.exit("Please keep this program in it's directory to load custom modules: " + error.message)

try:
    import numpy as np
except ImportError:
    sys.exit("Please install Numpy for this module: " + __name__)


_DISP_BREAK = '-----------------------------------------------------------------------------------------'
NA = 'NA'
Reporter = namedtuple('Reporter', ('deletions', 'insertions', 'mismatches', 'matches'))

def percent(num, total): # type: (int, int) -> float
    """Calculate a percent"""
    return round(num * 100 / total, 2)


def summarize(data, rounding=None): # type: (Iterable[int], Optional[int]) -> Tuple[int, float, float]
    '''Get the sum, mean, and standard deviation of a collection of data'''
    total = sum(data)
    avg = np.mean(data)
    std = np.std(data)
    if rounding:
        avg = round(avg, rounding)
        std = round(std, rounding)
    return total, float(avg), float(std)


def find_insertions(
        ref_seq, # type: str
        read_seq, # type: str
        num_reads # type: int
):
    # type: (...) -> Dict[int, List[int]], str, str, int
    """Find and log insertions"""
    insertions = defaultdict(list)
    insertion_gap_list = toolpack.find_gaps(seq=ref_seq) # type: List?
    if insertion_gap_list:
        nm_ins = len(insertion_gap_list) # type: int
        temp = 0 # type: int
        for gap in insertion_gap_list:
            position, gi_length = gap # type: int, int
            position = position - temp # type: int
            temp = temp + gi_length # type: int
            insertions[position].extend(itertools.repeat(gi_length, num_reads))
        to_remove = {m.start() for m in re.finditer('-', ref_seq)} # type: Set[int]
        ref_no_ins = ''.join((base for index, base in enumerate(ref_seq) if index not in to_remove)) # type: str
        read_no_ins = ''.join((base for index, base in enumerate(read_seq) if index not in to_remove)) # type: str
    else:
        nm_ins = 0 # type: int
        ref_no_ins = ref_seq # type: str
        read_no_ins = read_seq # type: str
    return dict(insertions), ref_no_ins, read_no_ins, nm_ins


def find_deletions(
        read_seq, # type: str
        num_reads, # type: int
        head, # type: int
        tail # type: int
):
    # type: (...) -> Dict[int, List[int]], int
    """Find and log deletions"""
    deletions = defaultdict(list)
    deletion_gap_list = toolpack.find_gaps(seq=read_seq, head=head, tail=tail)
    if deletion_gap_list:
        nm_del = len(deletion_gap_list)
        for gap in deletion_gap_list:
            position, gd_length = gap
            deletions[position].extend(itertools.repeat(gd_length, num_reads))
    else:
        nm_del = 0
    return dict(deletions), nm_del


def find_mismatches(
        ref_seq, # type: str
        read_seq, # type: str
        num_reads, # type: int
        read_head, # type: int
        read_tail # type: int
):
    # type: (...) -> Dict[int, List[str]], collections.Counter, int
    """Find and log mismatches
    'ref_seq' is the reference sequence without insertions
    'read_seq' is the aligned read sequence
    'num_reads' is the number of reads supporting this alignment
    'read_head' is where the aligned read starts (end of leading '-')
    'read_tail' is where the aligned read ends (start of trailing '-')"""
    mismatches = defaultdict(list)
    matches = Counter()
    mis_list, match_list = toolpack.get_mismatch(
        seq_a=ref_seq,
        seq_b=read_seq,
        head=read_head,
        tail=read_tail,
        matches=True
    )
    nm_mis = len(mis_list)
    for mismatch in mis_list:
        position, perm = mismatch
        source, snp = perm; del source
        mismatches[position].extend(itertools.repeat(snp, num_reads))
    for match in match_list:
        matches[match] += num_reads
    return dict(mismatches), matches, nm_mis



def run_analysis(
        reads_dict, # type: Dict[str, int]
        alignments, # type: Iterable[al.Alignment]
        score_threshold, # type: numpy.float64
        snp_index, # type: int
        target_snp # type: str
):
    # type: (...) -> (Reporter, Tuple[defaultdict[int, List[al.Alignment]]])
    """Analyze the alignment results
    'reads_dict' is a dictionary of unique sequences and the number of times they appear
    'alignments' is a diciontary where the values are lists of Alignments
    'score_threshold' is a floating point score threshold
    'snp_index', is an int representing the index where the SNP is
    'target_snp' is a signle character str reperesenting the target SNP for this experiment"""
    logging.info("Analyzing alignments")
    analyze_start = time.time() # type: float
    total_deletions = defaultdict(list)
    total_insertions = defaultdict(list)
    total_mismatches = defaultdict(list)
    total_matches = defaultdict(int)
    #   A tuple of results
    #   Unpack so I never have to refer to indecies
    class_reads = (defaultdict(list), defaultdict(list), defaultdict(list), defaultdict(list))
    hdr, nhej, no_edit, disc = class_reads
    out_idx = 0 # type: int
    #   Start doing things
    for alignment in alignments: # type: al.Alignment
        out_idx += 1
        num_reads = len(reads_dict[alignment.get_unaligned()]) # type: int
        #   Alignment is below score threshold, discard it
        if alignment.get_score() < score_threshold:
            alignment.set_stats(num_reads=num_reads, nm_del=0, nm_ins=0, nm_mis=0)
            disc[out_idx].append(alignment)
            continue
        #   Get get rid of head/tail gaps to better find indels and mismatches
        ref_head, ref_tail = toolpack.trim_interval(seq=alignment.get_aligned_reference()) # type: int, int
        aligned_ref = alignment.get_aligned_reference()[ref_head:ref_tail] # type: str
        aligned_read = alignment.get_aligned_read()[ref_head:ref_tail] # type: str
        read_head, read_tail = toolpack.trim_interval(aligned_read) # type: int, int
        #   Find insertions
        insertions, ref_no_ins, read_no_ins, nm_ins = find_insertions( # type: Dict[int, List[int]], str, str, int
            ref_seq=aligned_ref,
            read_seq=aligned_read,
            num_reads=num_reads
        )
        for position, ins_list in insertions.items():
            total_insertions[position].extend(ins_list)
        #   Find deletions
        deletions, nm_del = find_deletions( # type: Dict[int, List[int]], int
            read_seq=read_no_ins,
            num_reads=num_reads,
            head=read_head,
            tail=read_tail
        )
        for position, del_list in deletions.items():
            total_deletions[position].extend(del_list)
        #   Find matches and mismatches
        mismatches, matches, nm_mis = find_mismatches( # type: Dict[int, List[str]], collections.Counter, int
            ref_seq=ref_no_ins,
            read_seq=read_no_ins,
            num_reads=num_reads,
            read_head=read_head,
            read_tail=read_tail
        )
        for position, mis_list in mismatches.items():
            total_mismatches[position].extend(mis_list)
        for position, count in matches.items():
            total_matches[position] += count
        #   Classification
        #   Yay tuple unpacking!
        alignment.set_stats(num_reads=num_reads, nm_del=nm_del, nm_ins=nm_ins, nm_mis=nm_mis)
        if read_no_ins[snp_index] == target_snp:
            hdr[out_idx].append(alignment)
        elif not (nm_del and nm_ins and nm_mis):
            no_edit[out_idx].append(alignment)
        else:
            nhej[out_idx].append(alignment)
    report = Reporter( # type: Reporter
        deletions=total_deletions,
        insertions=total_insertions,
        mismatches=total_mismatches,
        matches=total_matches
    )
    logging.debug("Alignment analysis took %s seconds", round(time.time() - analyze_start, 3))
    return report, class_reads


def display_classification(
        read_classification, # type: Tuple[defaultdict[int, List[Alignment]]]
        total_reads, # type: int
        snp_position, # type: int
        ref_state, # type: str
        target_snp, # type: str
        fwd_score, # type: float
        rev_score, # type: float
        score_threshold, # type: float
        output_prefix # type: str
):
    # type: (...) -> None
    """Display the report of the read classifification
    'read_classification' is a tuple of four dictionaries where the values are
    lists of Alignment objects and the four dictionaries are (in order):
    'HDR', 'NHEJ', 'NO EDIT', and 'DISCARD' read classifications;
    'total_reads' is the total number of reads for the entire sample
    'snp_position' is the position of the SNP
    'ref_state' is the state of the reference sequence at 'snp_position'
    'target_snp' is the ideal SNP state
    'fwd_score' is the score of the forward alignment (from QC steps)
    'rev_score' is the score of the reverse alignment (from QC steps)
    'score_threshold' is the score threshold (from QC steps)
    'output_prefix' is the output directory + basename for the events report"""
    class_header = "################################################"
    warnings.simplefilter('error')
    full_class_name = output_prefix + '.classifications'
    logging.info("Writing full classification breakdown to %s", full_class_name)
    logging.warning(class_header)
    logging.warning("--------------Read Classifications--------------")
    num_unique = len( # type: int
        tuple(
            itertools.chain.from_iterable(
                (alignment_list for class_dict in read_classification for alignment_list in class_dict.values())
            )
        )
    )
    snp_header = ('##SNP', 'POS:%s' % (snp_position + 1), 'REF:%s' % ref_state, 'TEMPLATE:%s' % target_snp) # type: Tuple[str]
    read_header = ( # type: Tuple[str]
        '##READS',
        'TOTAL:%s' % total_reads,
        'UNIQUE:%s' % num_unique,
        'PERC_UNIQ:%s' % percent(num=num_unique, total=total_reads)
    )
    score_header = ('##SCORE', 'FWD:%s' % fwd_score, 'REV:%s' % rev_score, 'THESHOLD:%s' % score_threshold) # type: Tuple[str]
    category_header = ( # type: Tuple[str]
        '#TAG',
        'COUNT',
        'PERC_COUNT',
        'INS_EVENTS',
        'AVG_INS',
        'STD_DEV_INS',
        'DEL_EVENTS',
        'AVG_DEL',
        'STD_DEV_DEL',
        'MISMATCH_EVENTS',
        'AVG_MIS',
        'STD_DEV_MIS',
        'NO_INDELS',
        'PERC_NO_INDELS',
        'INS_ONLY',
        'PERC_INS_ONLY',
        'DEL_ONLY',
        'PERC_DEL_ONLY',
        'INDELS',
        'PERC_INDELS'
    )
    #   Four categories, and read_classification is a tuple, so need index, not names
    #   A dictionary of classifications, numbered for easy access
    iter_tag = { # type: Dict[int, str]
        0: 'HDR',
        1: 'NHEJ',
        2: 'NO_EDIT',
        3: 'DISCARD'
    }
    counted_total = 0 # type: int
    with open(full_class_name, 'w') as cfile:
        cfile.write('\t'.join(snp_header) + '\n')
        cfile.write('\t'.join(read_header) + '\n')
        cfile.write('\t'.join(score_header) + '\n')
        cfile.write('\t'.join(category_header) + '\n')
        cfile.flush()
        for index, tag in sorted(iter_tag.items(), key=lambda tup: tup[0]): # type: int, str
            #   Some holding values
            count = 0 # type: int
            event_lists = dict.fromkeys(('indels', 'insertions', 'deletions', 'mismatches'), list()) # type: Dict[str, List[int]]
            event_counts = dict.fromkeys(('none', 'deletions', 'insertions', 'indels'), 0) # type: Dict[str, int]
            for alignment in itertools.chain.from_iterable(read_classification[index].values()): # type: al.Alignment
                #   Create summaries
                num_reads, nm_del, nm_ins, nm_mis = alignment.get_stats()#; del nm_mis
                event_lists['indels'].append(nm_ins + nm_del)
                event_lists['insertions'].extend([nm_ins] * num_reads)
                event_lists['deletions'].extend([nm_del] * num_reads)
                event_lists['mismatches'].extend([nm_mis] * num_reads)
                if nm_ins > 0 and nm_del > 0:
                    event_counts['indels'] += num_reads
                elif nm_ins > 0 and nm_del <= 0:
                    event_counts['insertions'] += num_reads
                elif nm_del > 0 and nm_ins <= 0:
                    event_counts['deletions'] += num_reads
                else:
                    event_counts['none'] += num_reads
                count += num_reads
                counted_total += num_reads
            #   Do this to avoid Numpy warnings
            try:
                avg_indels = np.mean(event_lists['indels']) # type: numpy.float64
            except RuntimeWarning:
                avg_indels = 0
            perc_count = percent(num=count, total=total_reads)
            #   Display our summaries
            logging.warning("%s: count %s", iter_tag[index], count)
            logging.warning("%s: avg indels %s", iter_tag[index], round(avg_indels, 3))
            #   Reporting for HDR and NHEJ
            if tag in {iter_tag[0], iter_tag[1]}:
                total_ins, avg_ins, std_ins = summarize(data=event_lists['insertions'], rounding=2)
                total_del, avg_del, std_del = summarize(data=event_lists['deletions'], rounding=2)
                total_mis, avg_mis, std_mis = summarize(data=event_lists['mismatches'], rounding=2)
            else:
                total_ins, avg_ins, std_ins = NA, NA, NA
                total_del, avg_del, std_del = NA, NA, NA
                total_mis, avg_mis, std_mis = NA, NA, NA
            #   HDR-specific reporting
            if tag == iter_tag[0]:
                none_total, perc_none = event_counts['none'], percent(num=event_counts['none'], total=count)
                del_total, perc_del = event_counts['deletions'], percent(num=event_counts['deletions'], total=count)
                ins_total, perc_ins = event_counts['insertions'], percent(num=event_counts['insertions'], total=count)
                indel_total, perc_indel = event_counts['indels'], percent(num=event_counts['indels'], total=count)
            else:
                none_total, perc_none = NA, NA
                del_total, perc_del = NA, NA
                ins_total, perc_ins = NA, NA
                indel_total, perc_indel = NA, NA
            #   Assemble our output line
            out = ( # type: Tuple[Any]
                tag,
                count,
                perc_count,
                total_ins,
                avg_ins,
                std_ins,
                total_del,
                avg_del,
                std_del,
                total_mis,
                avg_mis,
                std_mis,
                none_total,
                perc_none,
                del_total,
                perc_del,
                ins_total,
                perc_ins,
                indel_total,
                perc_indel
            )
            out = map(str, out)
            cfile.write('\t'.join(out) + '\n')
            cfile.flush()
            #   Write full classifications
        if counted_total == total_reads:
            logging.warning("Classified all reads")
        else:
            logging.error("%s reads missing after classification", total_reads - counted_total)
        logging.warning(class_header)


def create_report(
        reporter, # type: Reporter
        reference, # type: str
        snp_index, # type: int
        output_prefix, # type: str
        interesting_only=True # type: bool
):
    # type: (...) -> None
    """Create the events report for the CRISPR program
    'reporter' is a namedtuple with four dictionaries of lists of alignments in the following order:
        'deletions', 'insertions', 'mismatches', 'matches'
    'reference' is the reference sequence
    'snp_index' is the position of the SNP
    'output_prefix' is the output directory + basename for the events report
    'interesting_only' tells us to report bases in the reference sequence where an event happens"""
    events_header = ( # type: Tuple[str]
        '#POS',
        'BASE',
        'COV',
        'DEL',
        'AVG_DEL',
        'INS',
        'AVG_INS',
        'A_MIS',
        'T_MIS',
        'C_MIS',
        'G_MIS'
    )
    #   Unpack our report
    deletions, insertions, mismatches, matches = reporter
    #   Start counting coverage
    cumulative_deletions = Counter() # type: collections.Counter
    for position, dist in deletions.items(): # type: int, List[int]
        for length in dist: # type: int
            for index in xrange(length): # type: int
                cumulative_deletions[position + index] += 1
    #   Prepare output file
    events_name = output_prefix + '.events'
    #   Start classifying events
    logging.info('Writing events log to %s', events_name)
    with open(events_name, 'w') as efile:
        efile.write('\t'.join(events_header) + '\n')
        for position, base in enumerate(reference): # type: int, str
            if position in insertions:
                count_ins = len(insertions[position]) # type: int
                mean_ins = round(np.mean(insertions[position]), 2) # type: float
            else:
                count_ins, mean_ins = 0, 0 # type: int, int
            if position in mismatches:
                nuc_counter = Counter(mismatches[position]) # type: collections.Counter
            else:
                nuc_counter = {base: 0 for base in 'ACGT'} # type: Dict[str, int]
            if position in deletions:
                count_del = len(deletions[position]) # type: int
                mean_del = round(np.mean(deletions[position]), 2) # type: float
            else:
                count_del, mean_del = 0, 0 # type: int, int
            covered = cumulative_deletions[position] + sum(nuc_counter.values()) # type: int
            if position in matches:
                covered += matches[position]
            if interesting_only:
                check = [ # type: List[int]
                    covered,
                    count_ins,
                    mean_ins,
                    count_del,
                    mean_del
                ]
                check.extend(nuc_counter.values())
                if sum(check) == 0:
                    continue
            msg = (
                position + 1,
                base,
                covered,
                count_del,
                mean_del,
                count_ins,
                mean_ins,
                nuc_counter['A'],
                nuc_counter['T'],
                nuc_counter['C'],
                nuc_counter['G']
            )
            if position == snp_index:
                efile.write(_DISP_BREAK + '\n')
            msg = map(str, msg) # type: Iterable[str]
            msg = '\t'.join(msg) # type: str
            efile.write(msg + '\n')
            if position == snp_index:
                efile.write(_DISP_BREAK + '\n')
