#!/usr/bin/env python2

"""Isolate reads"""

from __future__ import print_function

import sys
if sys.version_info.major is not 2 and sys.version_info.minor is not 7:
    sys.exit("Please use Python 2.7 for this module: " + __name__)


import time
import logging
from copy import deepcopy
from collections import defaultdict

try:
    # from genetic_toolpack import rvcomplement
    import genetic_toolpack as toolpack
except ImportError as error:
    sys.exit("Please keep this program in it's directory to load custom modules: " + error.message)


def load_seqs(raw_reads): # type: (Iterable[toolpack.Read]) -> Dict[str, List[toolpack.Read]]
    """Get a dictionary of unique sequences with their count"""
    logging.info('Loading reads')
    loading_start = time.time() # type: float
    reads_dict = defaultdict(list)
    for read in raw_reads: # type: toolpack.Read
        reads_dict[read.get_sequence()].append(read)
    logging.debug("Loading took %s seconds", round(time.time() - loading_start, 3))
    return dict(reads_dict)


def reverse_reads(reads_dict): # type: (Dict[str, List[toolpack.Read]]) -> Dict[str, List[toolpack.Read]]
    """Reverse a dictionary (counter) of reads"""
    logging.info("Reversing reads")
    reverse_start = time.time() # type: float
    rev_dict = defaultdict(list)
    for seq, reads_list in reads_dict.items(): # type: str, List[toolpack.Read]
        rev_reads = deepcopy(reads_list) # type: List[toolpack.Read]
        for read in rev_reads: # type: toolpack.Read
            read.rvcomplement()
        rev_dict[toolpack.rvcomplement(seq=seq)].extend(rev_reads)
    logging.debug("Reversing the reads took %s seconds", round(time.time() - reverse_start, 3))
    return dict(rev_dict)
