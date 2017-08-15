#!/usr/bin/env python

"""Tools for CRISPRonto"""

from __future__ import print_function
from __future__ import division

import sys

PYTHON_VERSION = sys.version_info.major

import os
import re
import time
import gzip
import logging
from collections import namedtuple

if PYTHON_VERSION is 2:
    from string import maketrans
    range = xrange
elif PYTHON_VERSION is 3:
    maketrans = str.maketrans
else:
    sys.exit("Please use Python 2 or 3 for this module: " + __name__)

try:
    from Bio.SeqIO.QualityIO import FastqGeneralIterator
except ImportError as error:
    sys.exit("Please install BioPython")


Read = namedtuple('Read', ('name', 'seq', 'qual'))
NamedSequence = namedtuple('NamedSequence', ('name', 'sequence'))

def load_fastq(fastq_file): # type: (str, Optional[str]) -> Tuple[Read]:
    """Load a FASTQ file"""
    logging.info("Reading in FASTQ file '%s'...", fastq_file)
    read_start = time.time() # type: float
    reads = [] # type: List[Read]
    if os.path.splitext(fastq_file)[-1] == '.gz':
        my_open = gzip.open
    else:
        my_open = open
    try:
        with my_open(fastq_file, 'rt') as ffile:
            for read in FastqGeneralIterator(ffile):
                name, seq, qual = read # type: str, str, str
                reads.append(Read(name=name, seq=seq, qual=qual))
    except:
        sys.exit(logging.critical("Cannot find or read FASTQ file '%s'", fastq_file))
    logging.debug("Reading in FASTQ file '%s' took %s seconds", fastq_file, round(time.time() - read_start, 3))
    return tuple(reads)


def load_seq(seq_file): # type: (str) -> NamedSequence
    """Load reference and template"""
    logging.info("Loading sequence file '%s'...", seq_file)
    load_start = time.time() # type: float
    output, name = str(), str() # type: str, str
    if os.path.splitext(seq_file)[-1] == '.gz':
        my_open = gzip.open
    else:
        my_open = open
    try:
        with my_open(seq_file, 'rt') as sfile:
            for line in sfile:
                if line.startswith('>'):
                    name += line.strip()
                    continue
                output += line.strip().replace(' ', '').upper()
    except:
        sys.exit(logging.critical("Cannot find or read sequence file '%s'", seq_file))
    if not name:
        name = os.path.basename(seq_file) # type:str
        if name.count('.') == 2:
            name = name.split('.')[0] # type: str
        else:
            name = os.path.splitext(name)[0] # type: str
    name = name.split(' ')[0].replace('>', '') # type: str
    logging.debug("Loading sequence took %s seconds", round(time.time() - load_start, 3))
    return NamedSequence(name=name, sequence=output)


def unpack(collection): # type: (Iterable[Any]) -> List[Any]
    """Unpack a series of nested lists, sets, or tuples"""
    result = [] # type: List
    for item in collection: # type: Any
        if hasattr(item, '__iter__') and not isinstance(item, str):
            result.extend(unpack(collection=item))
        else:
            result.append(item)
    return result

def reverse_complement(sequence): # type: (str) -> str
    """Reverse complement a nucleotide sequence"""
    table = maketrans('ACGT', 'TGCA') # type: Dict
    return sequence.translate(table)[::-1]


def get_mismatch(
        seq_a, # type: str
        seq_b, # type: str
        head=None, # type: Optional[int]
        tail=None, # type: Optional[int]
        matches=False # type: bool
):
    # type: (...) -> (List[int, Tuple[str]], Optional[List[int]])
    '''Get mismatches between two sequences after trimming head/tail gaps
    Optionally, get indecies where the sequences match as well'''
    if len(seq_a) != len(seq_b):
        raise ValueError("Sequences must be the same length")
    mis_list, match_list = list(), list() # type: Tuple[int, Tuple[str]], List[int]
    if head and tail:
        seq_range = range(head, tail + 1)
    else:
        seq_range = range(len(seq_a))
    for index in seq_range: # type: int
        if seq_a[index] == '-' or seq_b[index] == '-':
            #   Gap, not mismatch, continue
            continue
        if seq_a[index] == seq_b[index]:
            #   Match, not mismatch
            match_list.append(index)
            continue
        #   If neither a gap nor a match, it's a mismatch
        mis_list.append((index, (seq_a[index], seq_b[index])))
    if matches:
        return mis_list, match_list
    else:
        return mis_list


def find_gaps(seq, head=None, tail=None): # type: (str, Optional[int], Optional[int]) -> List[Tuple[int]]
    '''Return absolute position and length of all the gaps in a sequence
    Use known head and tail, or optionally find head and tail'''
    #Adjust the interval of study to discard head/tail gaps due to alignment
    if not (head and tail):
        head, tail = trim_interval(seq=seq) # type: int, int
    head = max(0, head) # type: int
    tail = min(tail, len(seq)) # type: int
    spans = (m.span() for m in re.finditer(r'(-+)', seq[head:tail]))
    gap_list = [(span[0] + head, span[1] - span[0]) for span in spans]
    return gap_list # type: List(Tuple[int])


def trim_interval(seq): # type: (str) -> (int, int)
    '''Define the interval discarding head/tail gaps caused by alignment'''
    # HEAD TRIMMING
    head = 0 # type: int
    while seq[head] == '-':
        head += 1
    #   Tail trimming
    tail = -1 # type: int
    while seq[tail] == '-':
        tail -= 1
    tail = (len(seq) + 1) + tail # type: int
    return head, tail


def side_trimmer(seq): # type: (str) -> str
    '''Trim only side gaps of an aligned sequence, return trimmed aligned sequence'''
    trimmed_seq = str() # type: str
    head, tail = trim_interval(seq=seq) # type: int, int
    if tail == -1:
        trimmed_seq = seq[head:] # type: str
    else:
        trimmed_seq = seq[head:tail] # type: str
    return trimmed_seq


def sim_seq(seq1, seq2): # type: (Iterable, Iterable) -> int
    '''Find up to which index seq2 is similar to seq1'''
    sim_index = 0 # type: int
    while seq1[sim_index] == seq2[sim_index]:
        sim_index += 1
    return sim_index
