#!/usr/bin/env python

"""Analysis of alignments for EdiTyper"""

from __future__ import division
from __future__ import print_function

import sys
PYTHON_VERSION = sys.version_info.major

import os
import re
import time
import logging
import warnings
import itertools
from math import floor, ceil
from collections import Counter, defaultdict, namedtuple

if PYTHON_VERSION is 2:
    import toolkit
    from itertools import imap as map
    range = xrange
elif PYTHON_VERSION is 3:
    from edityper import toolkit
else:
    raise SystemExit("Please use Python 2 or 3 for this module: " + __name__)

try:
    import numpy
except ImportError as error:
    raise SystemExit(error)


_DISP_BREAK = '-----------------------------------------------------------------------------------------'
NA = 'NA'
Reporter = namedtuple('Reporter', ('deletions', 'insertions', 'mismatches', 'matches'))

def _filter_to_dict(filtered_dict): # type: (Iterable[Tuple[Any, Any]]) -> Dict[Any, Any]
    return {pair[0]: pair[1] for pair in filtered_dict}


def _fastq_header(fastq_name, fastq_path): # type: (str, str) -> str
    header = ( # type: Tuple[str]
        '##FASTQ',
        'Name:%s' % fastq_name,
        'Path:%s' % fastq_path
    )
    return '\t'.join(header)


def _snp_header(snp_info): # type: (SNP) -> str
    header = ( # type: Tuple[str]
        '##SNP',
        'POS:%s' % (snp_info.position + 1),
        'REF:%s' % snp_info.reference,
        'TEMPLATE:%s' % snp_info.target
    )
    return '\t'.join(header)


def cummulative_deletions(deletions): # type: Dict[int, List[int]] -> Dict[int, int]
    """Calculate cummulative deletions"""
    cummul_del = defaultdict(int) # type: defaultdict
    for position, dist in deletions.items(): # type: int, List[int]
        for length in dist: # type: int
            for i in range(length): # type: int
                cummul_del[position + i] += 1
    return dict(cummul_del)


def calc_coverage(cummul_del, mismatches, matches): # type: (Dict[int, int], Dict[int, List[str]], Dict[int, int]) -> Dict[int, int]
    """Calculate coverage"""
    coverage = defaultdict(int) # type: defaultdict
    for base, cummul_count in cummul_del.items(): # type: int, int
        coverage[base] += cummul_count
    for base, mismatch_list in mismatches.items(): # type: int, List[str]
        coverage[base] += len(mismatch_list)
    for base, match_count in matches.items(): # type: int, int
        coverage[base] += match_count
    return coverage


def percent(num, total): # type: (int, int) -> float
    """Calculate a percent"""
    perc = num * 100 / total if total is not 0 else 0
    return float(round(perc, 2))


def summarize(data, rounding=None): # type: (Iterable[Union[int, float]], Optional[int]) -> Union[int, float], float, float
    '''Get the sum, mean, and standard deviation of a collection of data'''
    total = sum(data)
    if len(data) > 0:
        avg = numpy.mean(data)
        std = numpy.std(data)
    else:
        avg, std = 0, 0 # type: int, int
    if rounding:
        avg = round(avg, rounding)
        std = round(std, rounding)
    return total, float(avg), float(std)


def events_report(
        fastq_name, # type: str
        fastq_path, # type: str
        events, # type: Dict[str, defaultdict]
        cummul_del, # type: Dict[int, int]
        coverage, # type: Dict[int, int]
        reference, # type: str
        snp_info, # type: SNP
        output_prefix # type: str
):
    # type: (...) -> None
    """Create the events table"""
    logging.info("FASTQ %s: Creating events table...", fastq_name)
    events_start = time.time() # type: float
    #   Header information
    header = ( # type: Tuple[str]
        '#POS',
        'REF',
        'COV',
        'DEL',
        'AVG_DEL',
        'DCOUNT',
        'INS',
        'AVG_INS',
        'A',
        'T',
        'C',
        'G'
    )
    #   Create output file
    output_name = os.path.join(output_prefix, fastq_name + '.events')
    with open(output_name, 'w') as efile:
        logging.info("FASTQ %s: Writing events table to %s", fastq_name, output_name)
        efile.write(_fastq_header(fastq_name=fastq_name, fastq_path=fastq_path) + '\n')
        efile.write(_snp_header(snp_info=snp_info) + '\n')
        efile.write('\t'.join(header) + '\n')
        efile.flush()
        for index, base in enumerate(reference):
            #   Get the mismatches
            try:
                nucleotides = dict(Counter(events['mismatches'][index])) # type: Dict[str, int]
            except KeyError:
                nucleotides = dict.fromkeys(('A', 'C', 'G', 'T'), 0) # type: Dict[str, int]
            #   Get deletions
            deletions = events['deletions'].get(index, []) # type: List[int]
            deletion_count = len(deletions) # type: int
            avg_deletion = numpy.mean(deletions) if deletion_count else 0
            # Get insertions
            insertions = events['insertions'].get(index, []) # type: List[int]
            insertion_count = len(insertions) # type: int
            avg_insertion = numpy.mean(insertions) if insertion_count else 0
            #   Matches
            nucleotides[base] = events['matches'].get(index, 0) # type: int
            #   Assemble and write
            results = ( # type: Tuple[Any]
                index + 1,
                base,
                coverage.get(index, 0),
                deletion_count,
                round(avg_deletion, 2),
                cummul_del.get(index, 0),
                insertion_count,
                round(avg_insertion, 2),
                nucleotides.get('A', 0),
                nucleotides.get('T', 0),
                nucleotides.get('C', 0),
                nucleotides.get('G', 0)
            )
            results = map(str, results) # type: Tuple[str]
            efile.write('\t'.join(results))
            efile.write('\n')
            efile.flush()
    logging.debug("FASTQ %s: Creating events table took %s seconds", fastq_name, round(time.time() - events_start, 3))


def display_classification(
        fastq_name, # type: str
        fastq_path, # type: str
        classifications, # type: Tuple[Dict[str, Events]]
        unique_reads, # type: Mapping[str, int]
        snp_info, # type: SNP
        fwd_score, # type: float
        rev_score, # type: float
        score_threshold, # type: float
        output_prefix # type: str
):
    # type: (...) -> (int, Dict[str, int])
    """Display the report of the read classifification"""
    #   Make some headers for the display
    class_header = "################################################"
    pre_repeat = int(floor((len(class_header) - len(fastq_name)) / 2)) # type: int
    post_repeat = int(ceil((len(class_header) - len(fastq_name)) / 2)) # type: int
    name_header = ''.join(itertools.repeat('-', pre_repeat)) + fastq_name + ''.join(itertools.repeat('-', post_repeat))
    #   Create an output name
    output_name = os.path.join(output_prefix, fastq_name + '.classifications')
    logging.info("FASTQ %s: Writing full classification breakdown to %s", fastq_name, output_name)
    #   Quick statistics
    num_unique = len(unique_reads)
    total_reads = sum(unique_reads.values())
    #   Display our classifications
    logging.warning(class_header)
    logging.warning("--------------Read Classifications--------------")
    logging.warning(name_header)
    read_header = ( # type: Tuple[str]
        '##READS',
        'TOTAL:%s' % total_reads,
        'UNIQUE:%s' % num_unique,
        'PERC_UNIQ:%s' % percent(num=num_unique, total=total_reads)
    )
    score_header = ('##SCORE', 'FWD:%s' % fwd_score, 'REV:%s' % rev_score, 'THRESHOLD:%s' % score_threshold) # type: Tuple[str]
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
    #   Five categories, and read_classification is a tuple, so need index, not names
    #   A dictionary of classifications, numbered for easy access
    iter_tag = { # type: Dict[int, str]
        0: 'HDR',
        1: 'MIX',
        2: 'NHEJ',
        3: 'NO_EDIT',
        4: 'DISCARD'
    }
    counted_total = 0 # type: int
    hdr_indels = 0 # type: int
    total_counts = dict.fromkeys(iter_tag.values(), 0) # type: Dict[str, int]
    with open(output_name, 'w') as cfile:
        cfile.write(_fastq_header(fastq_name=fastq_name, fastq_path=fastq_path) + '\n')
        cfile.write(_snp_header(snp_info=snp_info) + '\n')
        cfile.write('\t'.join(read_header) + '\n')
        cfile.write('\t'.join(score_header) + '\n')
        cfile.write('\t'.join(category_header) + '\n')
        cfile.flush()
        for index, tag in sorted(iter_tag.items(), key=lambda tup: tup[0]): # type: int, str
            #   Some holding values
            count = 0 # type: int
            event_lists = defaultdict(list) # type: defaultdict[List]
            event_counts = dict.fromkeys(('none', 'deletions', 'insertions', 'indels'), 0) # type: Dict[str, int]
            for event in classifications[index].values(): # type: Event
                #   Create summaries
                event_lists['indels'].append(event.num_ins + event.num_del)
                event_lists['insertions'].extend([event.num_ins] * event.num_reads)
                event_lists['deletions'].extend([event.num_del] * event.num_reads)
                event_lists['mismatches'].extend([event.num_mis] * event.num_reads)
                if event.num_ins > 0 and event.num_del > 0:
                    event_counts['indels'] += event.num_reads
                elif event.num_ins > 0 and event.num_del <= 0:
                    event_counts['insertions'] += event.num_reads
                elif event.num_del > 0 and event.num_ins <= 0:
                    event_counts['deletions'] += event.num_reads
                else:
                    event_counts['none'] += event.num_reads
                count += event.num_reads
                counted_total += event.num_reads
            avg_indels = numpy.mean(event_lists['indels']) if event_lists['indels'] else 0
            if tag == 'DISCARD':
                perc_count = NA
            else:
                perc_count = percent(num=count, total=total_reads)
            #   Display our summaries
            logging.warning("%s: count %s", tag, count)
            logging.warning("%s: avg indels %s", tag, round(avg_indels, 3))
            #   Reporting for HDR/MIX and NHEJ
            if tag in {iter_tag[0], iter_tag[1], iter_tag[2]}:
                total_ins, avg_ins, std_ins = summarize(data=event_lists['insertions'], rounding=2)
                total_del, avg_del, std_del = summarize(data=event_lists['deletions'], rounding=2)
                total_mis, avg_mis, std_mis = summarize(data=event_lists['mismatches'], rounding=2)
            else:
                total_ins, avg_ins, std_ins = NA, NA, NA
                total_del, avg_del, std_del = NA, NA, NA
                total_mis, avg_mis, std_mis = NA, NA, NA
            #   HDR/MIX-specific reporting
            if tag in {iter_tag[0], iter_tag[1]}:
                none_total, perc_none = event_counts['none'], percent(num=event_counts['none'], total=count)
                del_total, perc_del = event_counts['deletions'], percent(num=event_counts['deletions'], total=count)
                ins_total, perc_ins = event_counts['insertions'], percent(num=event_counts['insertions'], total=count)
                indel_total, perc_indel = event_counts['indels'], percent(num=event_counts['indels'], total=count)
                hdr_indels = indel_total
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
            out = map(str, out) # type: Iterable[str]
            cfile.write('\t'.join(out) + '\n')
            cfile.flush()
            total_counts[tag] += count
            #   Write full classifications
        if counted_total == total_reads:
            logging.warning("Classified all reads")
        else:
            logging.error("%s reads missing after classification", total_reads - counted_total)
        logging.warning(class_header)
    return hdr_indels, total_counts
