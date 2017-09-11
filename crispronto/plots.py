#!/usr/bin/env python

"""Plotting utilities for the CRISPR program"""

from __future__ import division
from __future__ import print_function

import sys
PYTHON_VERSION = sys.version_info.major

import os
import time
import logging
from collections import defaultdict

try:
    if PYTHON_VERSION is 3:
        from crispronto import toolkit
    elif PYTHON_VERSION is 2:
        import toolkit
    else:
        raise SystemExit("Please use Python 2.7 or 3.5 or higher for this module: " + __name__)
except ImportError as error:
    raise SystemExit(error)

try:
    import matplotlib.pyplot as plt
    import matplotlib.patches as ptch
    import matplotlib.font_manager as fm
except ImportError as error:
    raise SystemExit(error)


_LOCUS_WIDTH = 0.5
_INDEL_COLOR = 'r'
_INS_COLOR = 'r'
_DEL_COLOR = 'g'
_MISMATCH_COLOR = 'b'
_XKCD = False

def _check_fonts(): # type: (None) -> bool
    humor = 'Humor-Sans.ttf' in map(os.path.basename, fm.findSystemFonts()) # type: bool
    if _XKCD and not humor:
        logging.error("Cannot find 'Humor-Sans.ttf' font. Please install and clear your matplotlib cache")
    return _XKCD and humor


def locus_plot(
        insertions, # type: Mapping[int, List[int]]
        deletions, # type: Mapping[int, List[int]]
        mismatches, # type: Mapping[int, List[str]]
        num_reads, # type: int
        fastq_name, # type: str
        output_prefix # type: str
):
    # type: (...) -> None
    """Make a locus plot"""
    logging.info("Making locus plot")
    locus_start = time.time() # type: float
    if _check_fonts():
        plt.xkcd()
    #   Make a name for our locus plot
    plot_name = os.path.join(output_prefix, fastq_name + '_locus.pdf') # type: str
    #   Create our subplots
    fig, ax = plt.subplots()
    #   Get our data ready for plotting
    mis_dict = {position: len(bases) for position, bases in mismatches.items()} # type: Dict[int, int]
    ins_dict = {position: len(counts) for position, counts in insertions.items()} # type: Dict[int, int]
    del_dict = {position: len(counts) for position, counts in deletions.items()} # type: Dict[int, int]
    #   Sort our data
    mis_counts = map(lambda tup: tup[1], sorted(mis_dict.items())) # type: List[int]
    ins_counts = map(lambda tup: tup[1], sorted(ins_dict.items())) # type: List[int]
    del_counts = map(lambda tup: tup[1], sorted(del_dict.items())) # type: List[int]
    #   Create the bar graphs
    ins_bars = ax.bar(
        left=sorted(insertions),
        height=ins_counts,
        width=_LOCUS_WIDTH,
        color=_INS_COLOR
    )
    del_bars = ax.bar(
        left=map(lambda x: x + _LOCUS_WIDTH, sorted(deletions)),
        height=del_counts,
        width=_LOCUS_WIDTH,
        color=_DEL_COLOR
    )
    mis_bars = ax.bar(
        left=map(lambda x: x + (_LOCUS_WIDTH * 2), sorted(mismatches)),
        height=mis_counts,
        width=_LOCUS_WIDTH,
        color=_MISMATCH_COLOR
    )
    #   Add title and legend
    plt.title(fastq_name)
    #   Add patches for colors
    ins_patch = ptch.Patch(color=_INS_COLOR, label='Insertions')
    del_patch = ptch.Patch(color=_DEL_COLOR, label='Deletions')
    mismatch_patch = ptch.Patch(color=_MISMATCH_COLOR, label='Mismatches')
    plt.legend(handles=(ins_patch, del_patch, mismatch_patch))
    #   Set the y limits and labels
    plt.ylim(0, num_reads)
    plt.ylabel('Number of Reads')
    #   Add a percent y axis
    ax2 = ax.twinx()
    ax2.set_ylabel('Percent')
    ax2.set_yticks(map(lambda x: round(x * 100), ax2.get_yticks()))
    #   Adjust the plot area to ensure everything is shown
    plt.tight_layout()
    #   Yield the plot
    logging.info("Saving plot to %s", plot_name)
    plt.savefig(plot_name, format='pdf')
    logging.debug("Making locus plot took %s seconds", round(time.time() - locus_start, 3))


def quality_plot(
        alignments, # type: Iterable[alignment.Alignment]
        output_prefix # type: str
):
    # type: (...) -> None
    """Make a violin plot of the alignment scores"""
    logging.info("Making quality scores plot")
    quality_start = time.time() # type: float
    if _check_fonts():
        plt.xkcd()
    plot_name = output_prefix + '_quality.pdf' # type: str
    alignments_by_fastq = defaultdict(list) # type: Mapping[str, List[alignment.Alignment]]
    for aligned in alignments: # type: alignment.Alignment
        alignments_by_fastq[aligned.source].append(aligned)
    #   Assemble our scores
    scores = {fastq: tuple(al.score for al in al_list) for fastq, al_list in alignments_by_fastq.items()} # type: Dict[str, Tuple[int]]
    #   Plot the scores
    fig, ax = plt.subplots(nrows=1, ncols=1)
    ax.violinplot(scores.values())
    #   Determine rotation of xtick text
    if len(scores) == 1:
        rotation = 'horizontal'
    else:
        rotation = 45
    #   Set labels
    plt.title("Alignment Score Distribution by FASTQ File")
    plt.ylabel('Alignment Score')
    plt.xlabel('FASTQ')
    ax.set_xticks(range(1, len(scores) + 1))
    ax.set_xticklabels(scores.keys(), rotation=rotation, fontsize='small')
    #   Adjust the plot area to ensure everything is shown
    plt.tight_layout()
    #   Yield the plot
    logging.info("Saving plot to %s", plot_name)
    plt.savefig(plot_name, format='pdf')
    logging.debug("Making quality scores plot took %s seconds", round(time.time() - quality_start, 3))
