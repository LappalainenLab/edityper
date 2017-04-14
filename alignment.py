#!/usr/bin/env python2

'''Alignment functions for the CRISPR program'''

from __future__ import print_function

import sys
if sys.version_info.major is not 2 and sys.version_info.minor is not 7:
    sys.exit("Please use Python 2.7 for this module: " + __name__)


import time
import logging
from collections import defaultdict, namedtuple

try:
    import NW_py as NW
    import genetic_toolpack as toolpack
except ImportError as error:
    sys.exit("Please keep this program in it's directory to load custom modules: " + error.message)


ReadSummary = namedtuple("ReadSummary", ('sequence', 'names'))

class Alignment(object):

    """An alignment just to make my life easier"""

    def __init__(
            self,
            ref_align, # type: str
            read_align, # type: str
            score, # type: str
            names=None # type: Optional[Iterable[str]]
    ):
        # type: (...) -> None
        self._ref = ref_align
        self._read = read_align
        self._score = score
        self._names = names
        self._unaligned = None
        self._num_reads = None
        self._nmdel = None
        self._nmins = None
        self._nmmis = None

    def __repr__(self): # type: (None) -> str
        return self._unaligned

    def set_unaligned(self, sequence): # type: (str) -> None
        """Provide an unaligned sequence"""
        self._unaligned = sequence

    def set_stats(self, num_reads, nm_del, nm_ins, nm_mis):  # type: (int, int, int, int) -> None
        """Set the stats"""
        self._num_reads = num_reads
        self._nmdel = nm_del
        self._nmins = nm_ins
        self._nmmis = nm_mis

    def get_aligned_reference(self): # type: (None) -> str
        """Get the aligned reference"""
        return self._ref

    def get_aligned_read(self): # type: (None) -> str
        """Get the aligned read"""
        return self._read

    def get_score(self): # type: (None) -> int
        """Get the alignment score"""
        return self._score

    def get_names(self): # type: (None) -> Tuple[str]
        """Get the name of this alignment"""
        return tuple(self._names)

    def get_stats(self): # type: (None) -> (int, int, int)
        """Get the stats for this alignment"""
        return self._num_reads, self._nmdel, self._nmins, self._nmmis

    def get_unaligned(self): # type: (None) -> str
        """Get the unaligned sequence"""
        return self._unaligned


def sort_reads_by_length(reads_dict): # type: (Dict[str, List[toolpack.Read]]) -> Dict[int, List[ReadSummary]]
    """Sort a list of reads by their length"""
    logging.info("Sorting reads by length")
    sort_start = time.time() # type: float
    reads_by_length = defaultdict(list)
    for seq, reads_list in reads_dict.items(): # type: str, List[toolpack.Read]
        length = len(seq) # type: int
        seq_info = ReadSummary(sequence=seq, names=tuple(read.get_readid() for read in reads_list)) # type: ReadSummary
        reads_by_length[length].append(seq_info)
        reads_by_length[length].sort(key=lambda summary: summary.sequence)
    logging.debug("Sorting the reads took %s seconds", round(time.time() - sort_start, 3))
    return reads_by_length


def align_recurse(
        reads_by_length, # type: Dict[int, List[ReadSummary]]
        reference, # type: str
        gap_open, # type: int
        gap_ext, # type: int
):
    # type: (...) -> Dict[int, List[Alignment]]
    """Align using the recurisve method"""
    logging.info("Aligning reads using the recursive method")
    align_start = time.time() # type: float
    #   Keep track of matrix re-used lines
    alignments = defaultdict(list) # type: defaultdict[List]
    reuse = 0 # type: int
    for length, reads_list in reads_by_length.items(): # type: int, List[ReadSummary]
        count, temp = 0, '' # type: int, str
        total = len(reads_list)
        for read in reads_list: # type: ReadSummary
            count += 1
            # seq = read.get_sequence() # type: str
            seq = read.sequence
            if not temp: # basically, first sequence to be aligned
                al_ref, al_read, score = NW.align_aff_mem(reference, seq, gap_open, gap_ext, 0, 0) # type: str, str, int
                aligned = Alignment(ref_align=al_ref, read_align=al_read, score=score, names=read.names) # type: Alignment
                aligned.set_unaligned(sequence=seq)
                alignments[length].append(aligned)
                temp = seq # type: str
                continue
            index = toolpack.sim_seq(seq1=al_ref, seq2=al_read) # type: int
            reuse += index
            if count is total:
                al_ref, al_read, score = NW.align_aff_mem(reference, seq, gap_open, gap_ext, index, 1) # type: str, str, int
            else:
                al_ref, al_read, score = NW.align_aff_mem(reference, seq, gap_open, gap_ext, index, 0) # type: str, str, int
            aligned = Alignment(ref_align=al_ref, read_align=al_read, score=score, names=read.names)
            aligned.set_unaligned(sequence=seq)
            alignments[length].append(aligned)
            temp = seq # type: str
    logging.debug("Aligning reads using the recursive method took %s seconds", round(time.time() - align_start, 3))
    return dict(alignments)
