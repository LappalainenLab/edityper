# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

#AUTHOR: ALEXANDRE YAHI
#AFFILATION: LAPPALAINEN LAB - New York Genome Center - Columbia University
#DATE: 2015-2016
#VERSION: 1.1

'''Python function for sequencing data - 1.1'''

from __future__ import print_function

import sys
if sys.version_info.major is not 2 and sys.version_info.minor is not 7:
    sys.exit("Please use Python 2.7 for this module: " + __name__)


import os
import gzip
import logging
from string import maketrans

def rvcomplement(seq): # type: (str) -> str
    '''Produces the reverse complement of a nucleotide sequence'''
    table = maketrans('ACGT', 'TGCA') # type: str
    return seq.translate(table)


class Read(object):
    """A FastQ read"""

    def __init__(self, readid, sequence, quality): # type: (str, str, str) -> None
        self._readid = str(readid)
        self._sequence = str(sequence)
        self._quality = str(quality)
        self._source = None
        if len(self._sequence) != len(self._quality):
            raise ValueError("'sequence' and 'quality' must be the same length!")

    def __repr__(self): # type: (None) -> str
        return self.get_readid()

    def __eq__(self, other): # type: (Union[Read, str]) -> bool
        if isinstance(other, Read):
            return self._readid == other._readid and self._sequence == other._sequence and self._quality == other._quality
        elif isinstance(other, str):
            return self._readid == other or self._sequence == other or self._quality == other
        else:
            return NotImplemented

    def __len__(self): # type: (None) -> int
        return len(self._sequence)

    def rvcomplement(self): # type: (None) -> None
        """Reverse complement the sequence for this read"""
        self._sequence = rvcomplement(seq=self._sequence)

    def set_source(self, source_file): # type: (str) -> None
        """Set the source file"""
        self._source = source_file

    def get_readid(self): # type: (None) -> str
        """Get the read ID"""
        return self._readid

    def get_sequence(self): # type: (None) -> str
        """Get the sequence of this read"""
        return self._sequence

    def get_quality(self): # type: (None) -> str
        """Get the quality scores for this read"""
        return self._quality

    def get_source(self): # type: (None) -> Optional[str]
        """Get the source file"""
        return self._source


class FastQ(object):
    """A FASTQ file (collection of reads)"""

    def __init__(self, fastq_file): # type: (str) -> None
        self._name = fastq_file
        self._reads = dict()

    def __repr__(self): # type: (None) -> str
        return self._name

    def __len__(self): # type: (None) -> int
        return len(self._reads)

    def add_read(self, read): # type: (Read) -> None
        """Add a Read to our FASTQ file"""
        if not isinstance(read, Read):
            raise TypeError("'read' must be of type Read")
        if read.get_readid() in self._reads:
            raise ValueError("Cannot have multiple Reads with the same read ID")
        read.set_source(source_file=self._name)
        self._reads[read.get_readid()] = read

    def get_read(self, readid): # type: (str) -> Read
        """Get a Read by read ID"""
        return self._reads[readid]

    def get_all_reads(self): # type: (None) -> List[Read]
        """Get all Reads from this FastQ"""
        return self._reads.values()

    def get_ids(self): # type: (None) -> Tuple(str)
        """Get all the read IDs from this FastQ"""
        return tuple(readid for readid in self._reads)

    def get_seqs(self): # type: (None) -> Tuple[str]
        """Get all the sequence data in this FASTQ"""
        return (read.get_sequence() for read in self._reads.values())


def load_fastq(fastq_file): # type: (str) -> FastQ
    '''Loading reads from a fastq file, list of lists including the +'''
    fastq, temp = FastQ(os.path.basename(fastq_file)), list()
    if 'gz' in fastq_file:
        fastq_handle = gzip.open(fastq_file, 'rb')
    else:
        fastq_handle = open(fastq_file, 'rb')
    for i, line in enumerate(fastq_handle):
        if i % 4 == 2:
            continue
        temp.append(line.decode().rstrip('\n'))
        if i % 4 == 3:
            this_read = Read(readid=temp[0], sequence=temp[1], quality=temp[2])
            fastq.add_read(this_read)
            temp = list()
    fastq_handle.close()
    return fastq


def load_seq(seq_file): # type: (str) -> (str, str)
    '''Load reference and template'''
    output, name = str(), str() # type: str, str
    if 'gz' in seq_file:
        my_open = gzip.open
    else:
        my_open = open
    with my_open(seq_file, 'rb') as sfile:
        for line in sfile:
            if line.startswith('>'):
                name += line.strip()
                continue
            output += line.strip().replace(' ', '').upper()
    if not name:
        name = os.path.basename(seq_file) # type:str
        if name.count('.') == 2:
            name = name.split('.')[0] # type: str
        else:
            name = os.path.splitext(name)[0] # type: str
    name = name.split(' ')[0].replace('>', '')
    return name, output


def get_mismatch(seq_a, seq_b, head=None, tail=None, matches=False):
    '''Get mismatches between two sequences after trimming head/tail gaps
    Optionally, get indecies where the sequences match as well'''
    if len(seq_a) != len(seq_b):
        raise ValueError("Sequences must be the same length")
    mis_list, match_list = list(), list() # type: Tuple[int, Tuple[str]], List[int]
    if head and tail:
        seq_range = xrange(head, tail + 1) # type: xrange
    else:
        seq_range = xrange(len(seq_a)) # type: xrange
    for index in seq_range: # type: int
        try:
            if seq_a[index] == '-' or seq_b[index] == '-':
                #   Gap, not mismatch, continue
                continue
            if seq_a[index] == seq_b[index]:
                #   Match, not mismatch
                match_list.append(index)
                continue
            #   If neither a gap nor a match, it's a mismatch
            mis_list.append((index, (seq_a[index], seq_b[index])))
        except IndexError:
            logging.error("No sequence found at base %i", index)
            break
    if matches:
        return mis_list, match_list
    else:
        return mis_list


def trim_interval(seq): # type: (str) -> (int, int)
    '''Define the interval discarding head/tail gaps caused by alignment'''
    # HEAD TRIMMING
    head = 0
    while seq[head] == '-':
        head += 1
    tail = -1
    while seq[tail] == '-':
        tail -= 1
    tail = len(seq) + tail
    return head, tail


def side_trimmer(seq): # type: (str) -> str
    '''Trim only side gaps of an aligned sequence, return trimmed aligned sequence'''
    trimmed_seq = str() # type: str
    head, tail = trim_interval(seq) # type: int, int
    if tail == -1:
        trimmed_seq = seq[head:]
    else:
        trimmed_seq = seq[head:tail+1]
    return trimmed_seq


def find_gaps(seq, head=None, tail=None): # type: (str, Optional[int], Optional[int]) -> List[List[int]]
    '''Return absolute position and length of all the gaps in a sequence
    Use known head and tail, or optionally find head and tail'''
    #Adjust the interval of study to discard head/tail gaps due to alignment
    if not (head and tail):
        head, tail = trim_interval(seq) # type: int, int
    gap_list = list() # type: List
    gap_open = False # type: bool
    # Temporary values
    index, length = 0, 0 # type: int, int
    for i in range(head, tail + 1):
        try:
            if seq[i] == '-': # GAP
                if not gap_open: # Open the gap
                    gap_open = True # type: bool
                    length = 1 # type: int
                    index = i # type: int
                else: #continue previous gap
                    length += 1
            else:
                if gap_open: # Close the gap
                    gap_list.append((index, length))
                    index, length = 0, 0 # type: int, int
                    gap_open = False # type: bool
                else: # Move along
                    continue
        except IndexError:
            logging.error("No sequence at %s", i)
            break
    return gap_list # type: List(Tuple[int])


def sim_seq(seq1, seq2): # type: (Iterable, Iterable) -> int
    '''Find up to which index seq2 is similar to seq1'''
    sim_index = 0 # type: int
    while seq1[sim_index] == seq2[sim_index]:
        sim_index += 1
    return sim_index
