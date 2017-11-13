#!/usr/bin/env python

"""Plotting utilities for EdiTyper"""

from __future__ import division
from __future__ import print_function

import sys
PYTHON_VERSION = sys.version_info.major

import os
import re
import time
import logging
from collections import defaultdict

try:
    if PYTHON_VERSION is 3:
        from edityper import toolkit
        from edityper.analysis import percent
    elif PYTHON_VERSION is 2:
        import toolkit
        from analysis import percent
        range = xrange
    else:
        raise SystemExit("Please use Python 2.7 or 3.5 or higher for this module: " + __name__)
except ImportError as error:
    raise SystemExit(error)

try:
    import matplotlib.pyplot as plt
    import matplotlib.patches as ptch
    import matplotlib.font_manager as fm
    from matplotlib.backends.backend_pdf import PdfPages
except ImportError as error:
    raise SystemExit(error)


_LOCUS_WIDTH = 0.5
_INDEL_COLOR = 'r'
_INS_COLOR = 'r'
_DEL_COLOR = 'g'
_MISMATCH_COLOR = 'b'
# _COVERAGE_COLOR = '#00000022'
# _COVERAGE_COLOR = '%08.x' % 0x00000022
_COVERAGE_COLOR = '#e5e5e8'
_THRESHOLD_LENGTH = 0.05
_CHUNK_DEFAULT = 5
_XKCD = False

def _check_fonts(): # type: (None) -> bool
    humor = 'Humor-Sans.ttf' in map(os.path.basename, fm.findSystemFonts()) # type: bool
    if _XKCD and not humor:
        logging.error("Cannot find 'Humor-Sans.ttf' font. Please install and clear your matplotlib cache")
    return _XKCD and humor


def _splitstr(string, cutoff=20): # type: (str, int) -> str
    if len(string) <= cutoff:
        return string
    pts = tuple(i.start() for i in re.finditer(r'([-_\.])', string)) # type: Tuple[int]
    midpoint = len(string) / 2 # type: float
    closest = min(pts, key=lambda x: abs(midpoint - x)) + 1 # type: int
    return string[:closest] + '\n' + string[closest:]


def ichunk(x, chunksize): # type: (Iterable[Any], int) -> Iterator[Iterable[Any]]
    """Iterate by chunks"""
    chunk_range = range(0, len(x), chunksize) # type: range
    for i in chunk_range: # type: int
        yield x[i:(i + chunksize)]


def locus_plot(
        insertions, # type: Mapping[int, List[int]]
        deletions, # type: Mapping[int, List[int]]
        mismatches, # type: Mapping[int, List[str]]
        coverage, # type: Mapping[int, int]
        num_reads, # type: int
        fastq_name, # type: str
        output_prefix # type: str
):
    # type: (...) -> None
    """Make a locus plot"""
    logging.info("FASTQ %s: Making locus plot", fastq_name)
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
    mis_counts = tuple(map(lambda tup: tup[1], sorted(mis_dict.items()))) # type: Tuple[int]
    ins_counts = tuple(map(lambda tup: tup[1], sorted(ins_dict.items()))) # type: Tuple[int]
    del_counts = tuple(map(lambda tup: tup[1], sorted(del_dict.items()))) # type: Tuple[int]
    cov_counts = tuple(map(lambda tup: tup[1], sorted(coverage.items()))) # type: Tuple[int]
    #   Create the bar graphs
    ax.bar( # Coverage
        left=sorted(coverage),
        height=cov_counts,
        width=1.0,
        color=_COVERAGE_COLOR,
        alpha=0.5
    )
    ax.bar( # Insertions
        left=sorted(insertions),
        height=ins_counts,
        width=_LOCUS_WIDTH,
        color=_INS_COLOR
    )
    ax.bar( # Deletions
        left=tuple(map(lambda x: x + _LOCUS_WIDTH, sorted(deletions))),
        height=del_counts,
        width=_LOCUS_WIDTH,
        color=_DEL_COLOR
    )
    ax.bar( # Mismatches
        left=tuple(map(lambda x: x + (_LOCUS_WIDTH * 2), sorted(mismatches))),
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
    coverage_patch = ptch.Patch(color=_COVERAGE_COLOR, label='Coverage')
    plt.legend(handles=(ins_patch, del_patch, mismatch_patch, coverage_patch))
    #   Set the y limits and labels
    ax.set_ylim(0, num_reads)
    ax.set_ylabel('Number of Reads')
    #   Add a percent y axis
    ax2 = ax.twinx()
    ax2.set_ylabel('Percent')
    ax2.set_yticks(tuple(map(lambda x: round(x * 100), ax2.get_yticks())))
    #   Zoomed plot
    fig_z, ax_z = plt.subplots()
    ax_z.bar( # Coverage
        left=sorted(coverage),
        height=cov_counts,
        width=1.0,
        color=_COVERAGE_COLOR,
        alpha=0.5
    )
    ax_z.bar( # Insertions
        left=sorted(insertions),
        height=ins_counts,
        width=_LOCUS_WIDTH,
        color=_INS_COLOR
    )
    ax_z.bar( # Deletions
        left=tuple(map(lambda x: x + _LOCUS_WIDTH, sorted(deletions))),
        height=del_counts,
        width=_LOCUS_WIDTH,
        color=_DEL_COLOR
    )
    ax_z.bar( # Mismatches
        left=tuple(map(lambda x: x + (_LOCUS_WIDTH * 2), sorted(mismatches))),
        height=mis_counts,
        width=_LOCUS_WIDTH,
        color=_MISMATCH_COLOR
    )
    #   Add title and legend
    plt.title(fastq_name)
    plt.legend(handles=(ins_patch, del_patch, mismatch_patch, coverage_patch))
    #   Set y label and add percent?
    ax_z.set_ybound(lower=0)
    ax_z.set_ylabel('Number of Reads')
    ax_z2 = ax_z.twinx()
    ax_z2.set_ylabel('Percent')
    ax_z2_percent = percent(num=max(ax_z.get_ylim()), total=num_reads) # type: float
    ax_z2.set_yticks(tuple(map(lambda x: round(float(x) * ax_z2_percent, 2), ax_z2.get_yticks())))
    #   Adjust the plot area to ensure everything is shown
    plt.tight_layout()
    #   Yield the plots
    with PdfPages(plot_name) as pdf:
        logging.info("FASTQ %s: Saving plot to %s", fastq_name, plot_name)
        for figure in (fig, fig_z):
            pdf.savefig(figure)
    plt.close('all')
    logging.debug("FASTQ %s: Making locus plot took %s seconds", fastq_name, round(time.time() - locus_start, 3))


def quality_plot(
        alignments, # type: Iterable[alignment.Alignment]
        thresholds, # type: Dict[str, float]
        output_prefix # type: str
):
    # type: (...) -> None
    """Make a violin plot of the alignment scores"""
    logging.info("Making quality scores plot")
    quality_start = time.time() # type: float
    if _check_fonts():
        plt.xkcd()
    plot_name = output_prefix + '_quality.pdf' # type: str
    scores_by_fastq = defaultdict(list) # type: Mapping[str, List[int]]
    for aligned in alignments: # type: alignment.Alignment
        scores_by_fastq[aligned.source].append(aligned.score)
    #   Get FASTQ names
    fastqs = sorted(set(scores_by_fastq.keys())) # type: List[str]
    #   Assemble our scores
    scores = {f: tuple(scores_by_fastq[f]) for f in fastqs} # type: Dict[str, Tuple[int]]
    with PdfPages(plot_name) as pdf:
        logging.info("Saving plot to %s", plot_name)
        #   Set chunksize
        chunksize = _CHUNK_DEFAULT if len(fastqs) > _CHUNK_DEFAULT else len(fastqs)
        logging.debug("Plotting %s FASTQ files per page", chunksize)
        #   Put at most 5 plots per page
        for fastq_chunk in ichunk(x=fastqs, chunksize=chunksize):
            logging.debug("Plotting quality scores for %s", ', '.join(fastq_chunk))
            fig, ax = plt.subplots(nrows=1, ncols=1)
            #   Plot the scores for this chunk
            scores_chunk = tuple(scores[f] for f in fastq_chunk) # type: Tuple[Tuple[int]]
            ax.violinplot(scores_chunk)
            #   Set labels for this page
            if len(scores) == 1:
                rotation = 'horizontal' # type: st
            else:
                rotation = 45 # type: int
            plt.title("Alignment Score Distribution by FASTQ File")
            plt.ylabel('Alignment Score')
            plt.xlabel('FASTQ')
            ax.set_xticks(range(1, len(scores) + 1))
            mod_names = tuple(_splitstr(string=f) for f in fastq_chunk) # type: Tuple[str]
            ax.set_xticklabels(mod_names, rotation=rotation, fontsize='small')
            #   Add threshold cutoff lines
            ax.hlines(
                y=tuple(thresholds[f] for f in fastq_chunk),
                xmin=tuple(map(lambda x: x - _THRESHOLD_LENGTH, ax.get_xticks())),
                xmax=tuple(map(lambda x: x + _THRESHOLD_LENGTH, ax.get_xticks()))
            )
            #   Ensure everything is shown
            plt.tight_layout()
            #   Save this figure
            pdf.savefig(fig)
    logging.debug("Making quality scores plot took %s seconds", round(time.time() - quality_start, 3))
