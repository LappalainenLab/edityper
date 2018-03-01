#!/usr/bin/env python

"""Plotting utilities for EdiTyper"""

from __future__ import division
from __future__ import print_function

import sys
PYTHON_VERSION = sys.version_info.major

# import gc
import os
import re
import time
import logging
import itertools
import subprocess
import pkg_resources
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

# try:
#     import matplotlib.pyplot as plt
#     import matplotlib.patches as ptch
#     import matplotlib.font_manager as fm
#     from matplotlib.backends.backend_pdf import PdfPages
# except ImportError as error:
#     raise SystemExit(error)


NAN = 'NaN'

_LOCUS_HEADER = (
    # 'Position',
    'Insertions',
    'Deletions',
    'Mismatches',
    'Coverage'
) # type: Tuple[str]

_LOCUS_WIDTH = 0.5 # type: float
_INDEL_COLOR = 'r' # type: str
_INS_COLOR = 'r' # type: str
_DEL_COLOR = 'g' # type: str
_MISMATCH_COLOR = 'b' # type: str
_COVERAGE_COLOR = '#e5e5e8' # type: str
_THRESHOLD_LENGTH = 0.05 # type: float
_CHUNK_DEFAULT = 5 # type: int
_XKCD = False # type: bool

# def _check_fonts(): # type: (None) -> bool
#     humor = 'Humor-Sans.ttf' in map(os.path.basename, fm.findSystemFonts()) # type: bool
#     if _XKCD and not humor:
#         logging.error("Cannot find 'Humor-Sans.ttf' font. Please install and clear your matplotlib cache")
#     return _XKCD and humor


def _splitstr(string, cutoff=20): # type: (str, int) -> str
    if len(string) <= cutoff:
        return string
    pts = tuple(i.start() for i in re.finditer(r'([-_\.])', string)) # type: Tuple[int]
    midpoint = len(string) / 2 # type: float
    closest = min(pts, key=lambda x: abs(midpoint - x)) + 1 # type: int
    return string[:closest] + '\n' + string[closest:]


# def ichunk(x, chunksize): # type: (Iterable[Any], int) -> Iterator[Iterable[Any]]
#     """Iterate by chunks"""
#     chunk_range = range(0, len(x), chunksize) # type: range
#     for i in chunk_range: # type: int
#         yield x[i:(i + chunksize)]


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
    #   Get our data ready for plotting
    mis_dict = {position: len(bases) for position, bases in mismatches.items()} # type: Dict[int, int]
    ins_dict = {position: len(counts) for position, counts in insertions.items()} # type: Dict[int, int]
    del_dict = {position: len(counts) for position, counts in deletions.items()} # type: Dict[int, int]
    all_counts = list()
    table_name = os.path.join(output_prefix, fastq_name + '_locus.txt')
    for position in range(max(itertools.chain(mis_dict, ins_dict, del_dict, coverage)) + 1): # type: int
        pos_counts = ( # type: Tuple[Union[int, str]]
            ins_dict.get(position, 0),
            del_dict.get(position, 0),
            mis_dict.get(position, 0),
            coverage.get(position, 0)
        )
        pos_counts = tuple(map(str, pos_counts)) # type: Tuple[str]
        all_counts.append(pos_counts)
    logging.debug("Making table of locus events")
    with open(table_name, 'w') as tfile:
        tfile.write('\t'.join(_LOCUS_HEADER))
        tfile.write('\n')
        tfile.flush()
        for line in all_counts:
            tfile.write('\t'.join(line))
            tfile.write('\n')
            tfile.flush()
    logging.debug("Assembling locus plot command")
    rscript_exec = toolkit.which('Rscript')
    lp_script = pkg_resources.resource_filename('edityper', 'locusPlot.R')
    cmd = [rscript_exec, lp_script, table_name, num_reads]
    cmd = map(str, cmd)
    logging.debug("Starting subprocess")
    logging.debug(cmd)
    proc = subprocess.Popen(cmd, stderr=subprocess.STDOUT, stdout=subprocess.PIPE)
    proc.wait()
    proc.communicate()
    logging.info("Plot located at ...")
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
    # if _check_fonts():
    #     plt.xkcd()
    # import code; code.interact(local=locals()); sys.exit()
    # plot_name = output_prefix + '_quality.pdf' # type: str
    scores_by_fastq = defaultdict(list) # type: Mapping[str, List[int]]
    for aligned in alignments: # type: alignment.Alignment
        scores_by_fastq[aligned.source].append(aligned.score)
    #   Get FASTQ names
    fastqs = sorted(set(scores_by_fastq.keys())) # type: List[str]
    #   Assemble our scores
    scores = {f: tuple(scores_by_fastq[f]) for f in fastqs} # type: Dict[str, Tuple[int]]
    table_name = output_prefix + '_quality.txt' # type: str
    max_num_scores = max(map(len, scores_by_fastq.values()))
    with open(table_name, 'w') as tfile:
        for fastq, scores in scores_by_fastq.items(): # type: str, List[int]
            scores = toolkit.unpack(collection=(scores, itertools.repeat(NAN, max_num_scores - len(scores)))) # type: Tuple[Union[int, str]]]
            scores = toolkit.unpack(collection=(thresholds[fastq], scores)) # type: Tuple[Untion[int, str, float]]
            scores = tuple(map(str, scores)) # type: Tuple[str]
            tfile.write(fastq + '\t')
            tfile.write('\t'.join(scores))
            tfile.write('\n')
            tfile.flush()
    rscript_exec = toolkit.which('Rscript')
    qp_script = pkg_resources.resource_filename('edityper', 'qualityPlot.R')
    cmd = [rscript_exec, qp_script, table_name]
    rc = subprocess.check_call(cmd)
    if rc != 0:
        logging.error("Something failed when making locus plot")
    # with PdfPages(plot_name) as pdf:
    #     logging.info("Saving plot to %s", plot_name)
    #     #   Set chunksize
    #     chunksize = _CHUNK_DEFAULT if len(fastqs) > _CHUNK_DEFAULT else len(fastqs)
    #     logging.debug("Plotting %s FASTQ files per page", chunksize)
    #     #   Put at most 5 plots per page
    #     for fastq_chunk in ichunk(x=fastqs, chunksize=chunksize):
    #         logging.debug("Plotting quality scores for %s", ', '.join(fastq_chunk))
    #         fig, ax = plt.subplots(nrows=1, ncols=1)
    #         #   Plot the scores for this chunk
    #         scores_chunk = tuple(scores[f] for f in fastq_chunk) # type: Tuple[Tuple[int]]
    #         ax.violinplot(scores_chunk)
    #         #   Set labels for this page
    #         if len(scores) == 1:
    #             rotation = 'horizontal' # type: st
    #         else:
    #             rotation = 45 # type: int
    #         plt.title("Alignment Score Distribution by FASTQ File")
    #         plt.ylabel('Alignment Score')
    #         plt.xlabel('FASTQ')
    #         ax.set_xticks(range(1, len(scores) + 1))
    #         mod_names = tuple(_splitstr(string=f) for f in fastq_chunk) # type: Tuple[str]
    #         ax.set_xticklabels(mod_names, rotation=rotation, fontsize='small')
    #         #   Add threshold cutoff lines
    #         ax.hlines(
    #             y=tuple(thresholds[f] for f in fastq_chunk),
    #             xmin=tuple(map(lambda x: x - _THRESHOLD_LENGTH, ax.get_xticks())),
    #             xmax=tuple(map(lambda x: x + _THRESHOLD_LENGTH, ax.get_xticks()))
    #         )
    #         #   Ensure everything is shown
    #         plt.tight_layout()
    #         #   Save this figure
    #         pdf.savefig(fig)
    #         gc.collect()
    # plt.close('all')
    logging.debug("Making quality scores plot took %s seconds", round(time.time() - quality_start, 3))
