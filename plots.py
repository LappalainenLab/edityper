#!/usr/bin/env python

"""Plotting utilities for the CRISPR program"""

from __future__ import division

import sys
if sys.version_info.major != 2 and sys.version_info.minor != 7:
    sys.exit("Please use Python 2.7 for this module: " + __name__)


import os
from collections import defaultdict

try:
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
        num_reads, # type: int
        output_prefix # type: str
):
    # type: (...) -> None
    """Make a locus plot"""
    if _check_fonts():
        plt.xkcd()
    #   Create our subplots
    fig, ax = plt.subplots()
    #   Get our data ready for plotting
    mis_dict = {position: len(bases) for position, bases in mismatches.items()} # type: Dict[int, int]
    indels = defaultdict(int) # type: Mapping[int, int]
    for this_dict in (insertions, deletions): # type: Mapping[int, List[int]]
        for position, lengths in this_dict.items(): # type: int, List[int]
            indels[position] += len(lengths)
    #   Sort our data
    indel_counts = map(lambda tup: tup[1], sorted(indels.items())) # type: List[int]
    mis_counts = map(lambda tup: tup[1], sorted(mis_dict.items())) # type: List[int]
    #   Create the bar graphs
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
    #   Set the y limits and labels
    plt.ylim(0, num_reads)
    plt.ylabel('Number of Reads')
    #   Add a percent y axis
    ax2 = ax.twinx()
    ax2.set_ylabel('Percent')
    ax2.set_yticks(map(lambda x: round(x * 100), ax2.get_yticks()))
    #   Yield the plot
    plt.show()


def quality_plot(
        alignments, # type: Iterable[alignment.Alignment]
        output_prefix # type: str
):
    # type: (...) -> None
    """Make a violin plot of the alignment scores"""
    if _check_fonts():
        plt.xkcd()
    #   Assemble our scores
    scores = tuple(al.get_score() for al in alignments)
    #   Plot the scores
    vlnplt = plt.violinplot(scores)
    #   Set labels
    plt.ylabel('Alignment Score')
    plt.xlabel('')
    #   Yield the plot
    plt.show()
