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
import itertools
import subprocess
import pkg_resources
from collections import defaultdict

try:
    if PYTHON_VERSION == 3:
        from edityper import toolkit
    elif PYTHON_VERSION == 2:
        import toolkit
        range = xrange
    else:
        raise SystemExit("Please use Python 2.7 or 3.5 or higher for this module: " + __name__)
except ImportError as error:
    raise SystemExit(error)


NAN = 'NaN'

_LOCUS_HEADER = ( # type: Tuple[str]
    # 'Position',
    'Insertions',
    'Deletions',
    'Mismatches',
    'Coverage'
)

_LOCUS_WIDTH = 0.5 # type: float
_INDEL_COLOR = 'r' # type: str
_INS_COLOR = 'r' # type: str
_DEL_COLOR = 'g' # type: str
_MISMATCH_COLOR = 'b' # type: str
_COVERAGE_COLOR = '#e5e5e8' # type: str
_THRESHOLD_LENGTH = 0.05 # type: float
_CHUNK_DEFAULT = 5 # type: int
_XKCD = False # type: bool


def _splitstr(string, cutoff=20): # type: (str, int) -> str
    if len(string) <= cutoff:
        return string
    pts = tuple(i.start() for i in re.finditer(r'([-_\.])', string)) # type: Tuple[int]
    midpoint = len(string) / 2 # type: float
    closest = min(pts, key=lambda x: abs(midpoint - x)) + 1 # type: int
    return string[:closest] + '\\n' + string[closest:]


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
    all_counts = list() # type: List[Tuple[str]]
    table_name = os.path.join(output_prefix, fastq_name + '_locus.txt') # type: str
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
    rscript_exec = toolkit.which('Rscript') # type: str
    lp_script = pkg_resources.resource_filename('edityper', 'locusPlot.R') # type: str
    cmd = [rscript_exec, lp_script, table_name, num_reads] # type: List[Union[int, str]]
    cmd = list(map(str, cmd)) # type: List[str]
    proc = subprocess.Popen(cmd, stderr=subprocess.PIPE, stdout=subprocess.PIPE) # type: subprocess.Popen
    proc.wait()
    _, err = proc.communicate() # type: str, str
    if err:
        logging.error(err)
    logging.info("Plot located at %s", os.path.splitext(table_name)[0] + '.pdf')
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
    scores_by_fastq = defaultdict(list) # type: Mapping[str, List[int]]
    for aligned in alignments: # type: alignment.Alignment
        scores_by_fastq[aligned.source].append(aligned.score)
    #   Get FASTQ names
    fastqs = sorted(set(scores_by_fastq.keys())) # type: List[str]
    #   Assemble our scores
    scores = {f: tuple(scores_by_fastq[f]) for f in fastqs} # type: Dict[str, Tuple[int]]
    table_name = output_prefix + '_quality.txt' # type: str
    max_num_scores = max(map(len, scores_by_fastq.values()))
    logging.debug("Making a table of quality scores")
    with open(table_name, 'w') as tfile:
        for fastq, scores in scores_by_fastq.items(): # type: str, List[int]
            scores = toolkit.unpack(collection=(scores, itertools.repeat(NAN, max_num_scores - len(scores)))) # type: Tuple[Union[int, str]]]
            scores = toolkit.unpack(collection=(thresholds[fastq], scores)) # type: Tuple[Untion[int, str, float]]
            scores = tuple(map(str, scores)) # type: Tuple[str]
            tfile.write(_splitstr(fastq) + '\t')
            tfile.write('\t'.join(scores))
            tfile.write('\n')
            tfile.flush()
    r_exec = toolkit.which('R') # type: str
    rscript_exec = toolkit.which('Rscript') # type: str
    qp_script = pkg_resources.resource_filename('edityper', 'qualityPlot.R') # type: str
    #   Check for writeable R library location
    lib_cmd = '%s --no-restore --no-save --slave -e "cat(.libPaths())"' % r_exec # type: str
    proc = subprocess.Popen(lib_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE) # type: subprocess.Popen
    out, _ = proc.communicate() # type: str, str
    out = out.split() # type: list[str]
    for libdir in out: # type: str
        if os.access(libdir, os.W_OK):
            break
    else:
        libdir = os.path.dirname(qp_script) # type: str
        if not os.access(libdir, os.W_OK):
            libdir = None # type: NoneType
    if libdir: # type: Optional[str]
        os.environ['R_LIBS'] = libdir
        #   Assemble plotting command
        plt_cmd = [rscript_exec, qp_script, table_name] # type: List[str]
        proc = subprocess.Popen(plt_cmd, stderr=subprocess.PIPE, stdout=subprocess.PIPE) # type: subprocess.Popen
        proc.wait()
        _, err = proc.communicate()
        if err:
            logging.error(err)
        logging.info("Quality scores plot can be found at %s", os.path.splitext(table_name)[0] + '.pdf')
    else:
        logging.error("Cannot make quality plots as there is no writeable R library location")
    logging.debug("Making quality scores plot took %s seconds", round(time.time() - quality_start, 3))
