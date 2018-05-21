#!/usr/bin/env python

"""Tools for EdiTyper"""

from __future__ import division
from __future__ import print_function

import sys

PYTHON_VERSION = sys.version_info.major

import os
import re
import time
import copy
import gzip
import math
import logging
import itertools
from collections import namedtuple

if PYTHON_VERSION == 2:
    from string import maketrans
    from itertools import imap as map
    from itertools import ifilter as filter
    range = xrange
    FileNotFoundError = IOError
elif PYTHON_VERSION == 3:
    maketrans = str.maketrans
else:
    raise SystemExit("Please use Python 2 or 3 for this module: " + __name__)


try:
    from Bio.SeqIO.QualityIO import FastqGeneralIterator
except ImportError as error:
    raise SystemExit(error)


try:
    from pympler import muppy, summary, tracker
    TRACKER = tracker.SummaryTracker()
except ImportError:
    _MEM_PROFILE = False
else:
    _MEM_PROFILE = True


_DO_PROFILE = False
Read = namedtuple('Read', ('name', 'seq', 'qual'))
NamedSequence = namedtuple('NamedSequence', ('name', 'sequence'))

class ExitPool(Exception):
    """Something happend and the pool needs to exit"""
    def __init__(self, msg):
        super(ExitPool, self).__init__(msg)
        self.msg = msg


class StrippedFormatter(logging.Formatter):
    """A formatter where all ANSI formatting is removed"""

    def __init__(self, *args, **kwargs):
        logging.Formatter.__init__(self, *args, **kwargs)

    def format(self, record): # type: (logging.LogRecord) -> str
        """Strip ANSI formatting from log messages"""
        message = logging.Formatter.format(self, record) # type: str
        while True:
            #   In Python, '\x1b' == '\033', so both codes for ANSI formatting are covered
            start = message.find('\x1b') # type: int
            #   If no ASI formatting is found break
            if start == -1:
                break
            #   Find the first 'm' after the ANSI code start
            #   and remove everything between and including
            #   the ANSI code start and the 'm'
            m_pos = message.find('m', start) # type: int
            message = message[:start] + message[m_pos + 1:]
        return message


class ColoredFormatter(logging.Formatter):
    """A colorized formatter for logging"""

    _colors = { # type: Dict[int, str]
        50: '\x1b[1m\x1b[31m', # CRITICAL: bold red
        40: '\x1b[31m', # ERROR: red
        30: '\x1b[33m', # WARNING: yellow
        20: '\x1b[32m', # INFO: green
        10: '\x1b[36m' # DEBUG: cyan
    }

    _default = '\x1b[0m' # Anything else: reset

    def __init__(self, *args, **kwargs):
        logging.Formatter.__init__(self, *args, **kwargs)

    def format(self, record): # type: (logging.LogRecord) -> str
        """Colorize log messages"""
        message = logging.Formatter.format(self, record) # type: str
        if sys.platform not in ('win32', 'cygwin'):
            color_level = min(self._colors.keys(), key=lambda level: abs(level - record.levelno)) # type: int
            color_level = min((color_level, record.levelno)) # type: int
            color = self._colors.get(color_level, self._default) # type: str
            message = color + message + self._default # type: str
        return message


class ColoredStreamHandler(logging.StreamHandler):
    """A colorized stream handler for logging"""

    _colors = {
        50: '\x1b[1m\x1b[31m', # CRITICAL: bold red
        40: '\x1b[31m', # ERROR: red
        30: '\x1b[33m', # WARNING: yellow
        20: '\x1b[32m', # INFO: green
        10: '\x1b[36m' # DEBUG: cyan
    }

    _default = '\x1b[0m' # Anything else: reset

    def emit(self, record): # type: (logging.LogRecord) -> None
        """Colorize to console"""
        color_record = copy.copy(record) # type: logging.LogRecord
        if sys.platform not in ('win32', 'cygwin'):
            color_level = min(self._colors.keys(), key=lambda level: abs(level - color_record.levelno)) # type: int
            color_level = min((color_level, color_record.levelno)) # type: int
            color = self._colors.get(color_level, self._default) # type: str
            color_record.msg = color + str(color_record.msg) + self._default # type: str
        logging.StreamHandler.emit(self, color_record)


def profile(diff=False): # type: (bool) -> None
    """A simple profiler using stuff from pympler"""
    if _DO_PROFILE and _MEM_PROFILE:
        if diff:
            TRACKER.print_diff()
        else:
            summary.print_(summary.summarize(muppy.get_objects()))
    elif _DO_PROFILE and not _MEM_PROFILE:
        logging.error("Could not find 'pympler' module, please install with pip and run again")
    else:
        pass


def which(program): # type: (str) -> str
    """Like UNIX which, returns the first location of a program given using your system PATH
    If full location to program is given, uses that first"""
    dirname, progname = os.path.split(program) # type: str, str
    syspath = tuple([dirname] + os.environ['PATH'].split(os.pathsep)) # type: Tuple[str]
    syspath = tuple(filter(None, syspath)) # type: tuple[str]
    progpath = map(os.path.join, syspath, itertools.repeat(progname, len(syspath))) # type: map[str]
    try:
        extensions = tuple([''] + os.environ.get('PATHEXT').split(os.pathsep)) # type: Tuple[str]
        progpath = map(lambda t: ''.join(t), itertools.product(progpath, extensions)) # type: map[str]
    except AttributeError:
        pass
    progpath = tuple(filter(lambda e: os.path.isfile(e) and os.access(e, os.X_OK), progpath)) # type: Tuple[str]
    if not progpath:
        raise ValueError("Cannot find program '%s' in your PATH" % program)
    return progpath[0]


def median(x): # type: (Iterable[Union[int, float]]) -> float
    """Get the median of a dataset"""
    if isinstance(x, (int, float)):
        return float(x)
    if not all(map(lambda i: isinstance(i, (int, float)), x)):
        raise ValueError("'x' must be an iterable of integers or floats")
    x_sorted = sorted(x)
    midpoint = int(len(x_sorted) / 2)
    if len(x_sorted) % 2:
        return float(x_sorted[midpoint])
    elif len(x_sorted) == 2:
        lower, upper = x_sorted
    else:
        lower, upper = x_sorted[midpoint:(midpoint + 2)]
    return upper - ((upper - lower) / 2)


def mean(x): # type: (Iterable[Union[int, float]]) -> float
    """Get the mean of a dataset"""
    if not x:
        return 0.0
    if isinstance(x, (int, float)):
        x = (x,) # type: Tuple[float]
    if not all(map(lambda i: isinstance(i, (int, float)), x)):
        raise ValueError("'x' must be an iterable of integers or floats")
    return sum(x) / len(x)


def variance(x): # type: (Iterable[Union[int, float]]) -> float
    """Get the variance of a dataset"""
    if not x:
        return 0.0
    if isinstance(x, (int, float)):
        x = (x,) # type: Tuple[float]
    if not all(map(lambda i: isinstance(i, (int, float)), x)):
        raise ValueError("'x' must be an iterable of integers or floats")
    avg = mean(x=x)
    return sum(map(lambda i: (i - avg) ** 2, x)) / len(x)


def stdev(x): # type: (Iterable[Union[int, float]]) -> float
    """Get the standard deviation of a dataset"""
    return math.sqrt(variance(x=x))


def norm_pdf(x, mu=0, sigma=1): # type: (float, float, float) -> float
    """Calculate a the probability density of a normal distribution"""
    return math.exp(-(((x - mu) ** 2) / (2 * (sigma ** 2)))) / math.sqrt(2 * math.pi * (sigma ** 2))


def find_fastq(directory): # type: (str) -> Tuple[str]
    """Find FASTQ files in a directory"""
    directory = os.path.abspath(directory) # type: str
    if not os.path.isdir(directory):
        raise SystemExit(logging.critical("Cannot find FASTQ directory %s", directory))
    fastqs = re.findall(r'(.*\.(fq|fastq)(\.gz)?)', '\n'.join(os.listdir(directory))) # type: List[Tuple[str, str, str]]
    fastqs = unpack(collection=fastqs) # type: Tuple[str]
    fastqs = itertools.islice(fastqs, 0, None, 3) # type: itertools.islice
    fastqs = map(lambda name: os.path.join(directory, name), fastqs) # type: map
    fastqs = tuple(fastqs) # type: Tuple[str]
    if not fastqs:
        raise SystemExit(logging.critical("No FASTQs found in %s", directory))
    return fastqs


def load_fastq(fastq_file): # type: (str, Optional[str]) -> Tuple[Read]:
    """Load a FASTQ file"""
    logging.info("Reading in FASTQ file '%s'...", fastq_file)
    read_start = time.time() # type: float
    reads = [] # type: List[Read]
    if os.path.splitext(fastq_file)[-1] == '.gz':
        my_open = gzip.open # type: function
    else:
        my_open = open # type: function
    try:
        with my_open(fastq_file, 'rt') as ffile:
            for read in FastqGeneralIterator(ffile): # type: Tuple[str, str, str]
                name, seq, qual = read # type: str, str, str
                reads.append(Read(name=name, seq=seq.upper(), qual=qual))
    except (IOError, FileNotFoundError):
        raise ExitPool(logging.critical("Cannot find or read FASTQ file '%s'", fastq_file))
    logging.debug("Reading in FASTQ file '%s' took %s seconds", fastq_file, round(time.time() - read_start, 3))
    return tuple(reads)


def load_seq(seq_file, chrom=None): # type: (str, Optional[str]) -> NamedSequence
    """Load reference and template"""
    logging.info("Loading sequence file '%s'...", seq_file)
    load_start = time.time() # type: float
    output, name = str(), str() # type: str, str
    if os.path.splitext(seq_file)[-1] == '.gz':
        my_open = gzip.open
    else:
        my_open = open
    try:
        with my_open(seq_file, 'rt') as sfile:
            for line in sfile:
                if line.startswith('>'):
                    name += line.strip()
                    continue
                output += line.strip().replace(' ', '').upper()
    except (IOError, FileNotFoundError):
        raise SystemExit(logging.critical("Cannot find or read sequence file '%s'", seq_file))
    if chrom:
        name = chrom # type: str
    if not name:
        name = os.path.basename(seq_file) # type:str
        if name.count('.') == 2:
            name = name.split('.')[0] # type: str
        else:
            name = os.path.splitext(name)[0] # type: str
    name = name.split(' ')[0].replace('>', '') # type: str
    logging.debug("Loading sequence took %s seconds", round(time.time() - load_start, 3))
    return NamedSequence(name=name, sequence=output)


def unpack(collection): # type: (Iterable[Any]) -> Tuple[Any]
    """Unpack a series of nested lists, sets, or tuples"""
    result = [] # type: List
    for item in collection: # type: Any
        if hasattr(item, '__iter__') and not isinstance(item, str):
            result.extend(unpack(collection=item))
        else:
            result.append(item)
    return tuple(result)


def reverse_complement(sequence): # type: (str) -> str
    """Reverse complement a nucleotide sequence"""
    table = maketrans('ACGT', 'TGCA') # type: Dict
    return sequence.translate(table)[::-1]


def get_mismatch(
        seq_a, # type: str
        seq_b, # type: str
        head=None, # type: Optional[int]
        tail=None, # type: Optional[int]
        matches=False # type: bool
):
    # type: (...) -> (List[Tuple[int, Tuple[str]]], Optional[List[int]])
    '''Get mismatches between two sequences after trimming head/tail gaps
    Optionally, get indecies where the sequences match as well'''
    if len(seq_a) != len(seq_b):
        raise ValueError("Sequences must be the same length")
    mis_list, match_list = list(), list() # type: List[Tuple[int, Tuple[str, str]]], List[int]
    if head and tail:
        seq_range = range(head, tail + 1) # type: Iterable[int]
    else:
        seq_range = range(len(seq_a)) # type: Iterable[int]
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
            break
    if matches:
        return mis_list, match_list
    else:
        return mis_list


def find_gaps(seq, head=None, tail=None): # type: (str, Optional[int], Optional[int]) -> List[Tuple[int, int]]
    '''Return absolute position and length of all the gaps in a sequence
    Use known head and tail, or optionally find head and tail'''
    # Adjust the interval of study to discard head/tail gaps due to alignment
    if not (head and tail):
        head, tail = trim_interval(seq=seq) # type: int, int
    head = max(0, head) # type: int
    tail = min(tail, len(seq)) # type: int
    spans = (m.span() for m in re.finditer(r'(-+)', seq[head:tail])) # type: generator[Tuple[int, int]]
    gap_list = [(span[0] + head, span[1] - span[0]) for span in spans] # type: List[Tuple[int, int]]
    return gap_list


def trim_interval(seq): # type: (str) -> (int, int)
    '''Define the interval discarding head/tail gaps caused by alignment'''
    # HEAD TRIMMING
    head = 0 # type: int
    while seq[head] == '-':
        head += 1
    #   Tail trimming
    tail = -1 # type: int
    while seq[tail] == '-':
        tail -= 1
    tail = (len(seq) + 1) + tail # type: int
    return head, tail


def side_trimmer(seq): # type: (str) -> str
    '''Trim only side gaps of an aligned sequence, return trimmed aligned sequence'''
    head, tail = trim_interval(seq=seq) # type: int, int
    return seq[head:tail] # type: str


def sim_seq(seq1, seq2): # type: (Iterable, Iterable) -> int
    '''Find up to which index seq2 is similar to seq1'''
    sim_index = 0 # type: int
    while seq1[sim_index] == seq2[sim_index]:
        sim_index += 1
    return sim_index
