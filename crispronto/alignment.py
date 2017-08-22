#!/usr/bin/env python

'''Alignment functions for CRISPRonto'''

from __future__ import division
from __future__ import print_function

import sys
PYTHON_VERSION = sys.version_info.major

import time
import logging
from collections import defaultdict

try:
    if PYTHON_VERSION is 3:
        from crispronto import toolkit
        from crispronto import nw_align
    elif PYTHON_VERSION is 2:
        import toolkit
        import nw_align
    else:
        sys.exit("Please use Python 2.7 or 3.5 or higher for this module: " + __name__)
except ImportError:
    raise

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
        self._ref = ref_align # type: str
        self._read = read_align # type: str
        self._score = score # type: str
        self._names = names # type: Iterable[str]
        self._source = None # type: Optional[str]
        self._unaligned = None # type: Optional[str]
        self._num_reads = None # type: Optional[int]
        self._nmdel = None # type: Optional[int]
        self._nmins = None # type: Optional[int]
        self._nmmis = None # type: Optional[int]

    def __repr__(self): # type: (None) -> str
        return self._unaligned if self._unaligned else self._read

    def _set_unaligned(self, sequence): # type: (str) -> None
        self._unaligned = sequence

    def _set_names(self, *args): # type: (Iterable[str]) -> None
        args = toolkit.unpack(collection=args)
        self._names = args

    def _set_source(self, source): # type: (str) -> None
        self._source = source

    def _get_aligned_reference(self): # type: (None) -> str
        return self._ref

    def _get_aligned_read(self): # type: (None) -> str
        return self._read

    def _get_score(self): # type: (None) -> int
        return self._score

    def _get_names(self): # type: (None) -> Tuple[str]
        return tuple(self._names)

    def _get_unaligned(self): # type: (None) -> str
        return self._unaligned

    def _get_source(self): # type: (None) -> str
        return self._source

    reference = property(fget=_get_aligned_reference, doc='Aligned reference sequence')
    read = property(fget=_get_aligned_read, doc='Aligned read sequence')
    score = property(fget=_get_score, doc='Alignment score')
    unaligned = property(fget=_get_unaligned, fset=_set_unaligned, doc='Set an unaligned sequence')
    names = property(fget=_get_names, fset=_set_names, doc='Names of supporting reads')
    source = property(fget=_get_source, fset=_set_source, doc='Name of source file')


def sort_reads_by_length(reads, fastq_name): # type: (Iterable[str], str) -> Dict[int, List[str]]
    """Sort a list of reads by their length"""
    logging.info("FASTQ %s: Sorting reads by length", fastq_name)
    sort_start = time.time()
    reads_by_length = defaultdict(list)
    for read in reads:
        reads_by_length[len(read)].append(read)
    for length in reads_by_length:
        reads_by_length[length].sort()
    logging.debug("FASTQ %s: Sorting reads took %s seconds", fastq_name, round(time.time() - sort_start, 3))
    return reads_by_length


def align_recurse(
        fastq_name, # type: str
        reads_by_length, # type: Dict[int, List[str]]
        reference, # type: str
        gap_open, # type: int
        gap_extension, # type: int
):
    # type: (...) -> Dict[int, List[Alignment]]
    """Align using the recurisve method"""
    #   Keep track of matrix re-used lines
    logging.info("FASTQ %s: Aligning reads", fastq_name)
    alignment_start = time.time()
    alignments = defaultdict(list) # type: defaultdict[List]
    reuse = 0 # type: int
    for length, reads_list in reads_by_length.items(): # type: int, List[str]
        count, temp = 0, '' # type: int, str
        total = len(reads_list)
        for read in reads_list: # type: read
            count += 1
            # seq = summary.sequence # type: str
            # names = tuple(read.get_readid() for read in summary.reads) # type: Tuple[str]
            # fastqs = {read.get_source() for read in summary.reads} # type: Set[str]
            if not temp: # basically, first sequence to be aligned
                al_ref, al_read, score = nw_align.align_aff_mem(
                    seq_1=reference,
                    seq_2=read,
                    gap_op=gap_open,
                    gap_ext=gap_extension,
                    sim=0,
                    terminate=0
                ) # type: str, str, int
                aligned = Alignment(ref_align=al_ref, read_align=al_read, score=score) # type: Alignment
                aligned.unaligned = read
                aligned.source = fastq_name
                alignments[length].append(aligned)
                temp = read # type: str
                continue
            index = toolkit.sim_seq(seq1=al_ref, seq2=al_read) # type: int
            reuse += index
            if count == total:
                al_ref, al_read, score = nw_align.align_aff_mem(
                    seq_1=reference,
                    seq_2=read,
                    gap_op=gap_open,
                    gap_ext=gap_extension,
                    sim=index,
                    terminate=1
                ) # type: str, str, int
            else:
                al_ref, al_read, score = nw_align.align_aff_mem(
                    seq_1=reference,
                    seq_2=read,
                    gap_op=gap_open,
                    gap_ext=gap_extension,
                    sim=index,
                    terminate=0
                ) # type: str, str, int
            aligned = Alignment(ref_align=al_ref, read_align=al_read, score=score)
            aligned.unaligned = read
            aligned.source = fastq_name
            alignments[length].append(aligned)
            temp = read # type: str
    logging.debug("FASTQ %s: Alignment took %s seconds", fastq_name, round(time.time() - alignment_start, 3))
    return dict(alignments)
