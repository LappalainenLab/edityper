#!/usr/bin/env python

"""Plotting utilities for the CRISPR program"""

from __future__ import division

import sys
if sys.version_info.major != 2 and sys.version_info.minor != 7:
    sys.exit("Please use Python 2.7 for this module: " + __name__)


import os
from collections import defaultdict, Counter

# from analysis import Reporter

try:
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.font_manager as fm
except ImportError as error:
    sys.exit("Please install Numpy and Matplotlib for this module: " + error.message)


_LOCUS_WIDTH = 0.5
_INDEL_COLOR = 'r'
_MISMATCH_COLOR = 'b'
_XKCD = False

_check_fonts = lambda: 'Humor-Sans.ttf' in map(os.path.basename, fm.findSystemFonts()) if _XKCD else _XKCD


def locus_plot(
        insertions, # type: Mapping[int, List[int]]
        deletions, # type: Mapping[int, List[int]]
        mismatches, # type: Mapping[int, List[str]]
        reference, # type: str
        output_prefix # type: str
):
    # type: (...) -> None
    """Make a locus plot"""
    if _check_fonts():
        plt.xkcd()
    indels = Counter() # type: Counter
    fig, ax = plt.subplots()
    # indels = dict.fromkeys(xrange(len(reference)), 0)
    # mis_dict = {position: len(mismatches[position]) if position in mismatches else 0 for position in xrange(len(reference))}
    mis_dict = {position: len(bases) for position, bases in mismatches.items()}
    # xlocations = np.arange(len(reference)) # type: numpy.ndarray
    for this_dict in (insertions, deletions): # type: Mapping[int, List[int]]
        for position, counts in this_dict.items(): # type: int, List[int]
            indels[position] += sum(counts)
    indel_counts = map(lambda tup: tup[1], sorted(indels.items())) # type: List[int]
    mis_counts = map(lambda tup: tup[1], sorted(mis_dict.items())) # type: List[int]
    indel_bars = ax.bar(
        left=sorted(indels),
        height=indel_counts,
        width=_LOCUS_WIDTH,
        color=_INDEL_COLOR
    )
    mis_bars = ax.bar(
        left=map(lambda x: x + _LOCUS_WIDTH, sorted(mismatches)),
        height=mis_counts,
        width=_LOCUS_WIDTH,
        color=_MISMATCH_COLOR
    )
    max_count = max(indel_counts + mis_counts)
    plt.ylim(0, max_count + (max_count / 8))
    plt.ylabel('Number of Reads')
    plt.show()


def quality_plot(
        alignments, # type: Iterable[alignment.Alignment]
        output_prefix # type: str
):
    # type: (...) -> None
    """Make a violin plot of the alignment scores"""
    if _check_fonts():
        plt.xkcd()
    scores = tuple(al.get_score() for al in alignments)
    vlnplt = plt.violinplot(scores)
    plt.ylabel('Alignment Score')
    plt.xlabel('')
    plt.show()
