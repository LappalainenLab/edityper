#!/usr/bin/env python

"""CRISPRonto"""

from __future__ import print_function
from __future__ import division

import sys
PYTHON_VERSION = sys.version_info.major

import time
import logging
import signal
import itertools
from multiprocessing import Lock

try:
    if PYTHON_VERSION is 3:
        from crispronto import toolkit
        from crispronto import nw_align
        from crispronto.arguments import make_argument_parser
        from multiprocessing import pool as Pool
    elif PYTHON_VERSION is 2:
        import toolkit
        import nw_align
        from arguments import make_argument_parser
        from multiprocessing import Pool
        from itertools import izip as zip
        from itertools import imap as map
        range = xrange
    else:
        sys.exit('WTF MATE')
except ImportError as error:
    sys.exit("FAIL")


def _set_verbosity(level): # type: (str) -> int
    level = level.upper()
    if level == 'DEBUG':
        log_level = logging.DEBUG # type: int
    elif level == 'INFO':
        log_level = logging.INFO # type: int
    elif level == 'WARNING':
        log_level = logging.WARNING # type: int
    elif level == 'ERROR':
        log_level = logging.ERROR # type: int
    elif level == 'CRITICAL':
        log_level = logging.CRITICAL # type: int
    else:
        raise ValueError("'level' must be one of 'DEBUG', 'INFO', 'WARNING', 'ERROR', or 'CRITICAL'")
    return log_level


def _check_suppressions(suppressions): # type: (Dict[str, bool]) -> bool
    suppress_tables = suppressions['suppress_events'] and suppressions['suppress_classification']
    return suppressions['suppress_plots'] and suppressions['suppress_sam'] and (suppressions['suppress_tables'] or suppress_tables)


def main():
    """CRISPronto"""
    #   Setup CRISPRonto
    #   Parse arguments
    parser = make_argument_parser() # type: argparse.ArgumentParser
    if not sys.argv[1:]:
        sys.exit(parser.print_help())
    args = {key: value for key, value in vars(parser.parse_args()).items() if value is not None} # type: Dict[str, Any]
    #   Setup logger
    logging.basicConfig(
        format='%(asctime)s %(levelname)s:\t%(message)s',
        stream=args['logfile'],
        level=_set_verbosity(level=args['verbosity']),
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    #   Check suppression values
    suppressions = {key: value for key, value in args.items() if 'suppress' in key} # type: Dict[str, bool]
    if _check_suppressions(suppressions=suppressions): # All output suppressed? Error
        sys.exit(logging.critical("All output suppressed, not running"))
    if suppressions['suppress_sam']: # Suppressed SAM output?
        logging.warning("SAM output suppressed, not writing SAM file")
    if suppressions['suppress_events'] or suppressions['suppress_tables']: # Suppressed events table?
        logging.warning("Events output suppressed, not writing events table")
    if suppressions['suppress_classification'] or suppressions['suppress_tables']: # Suppressed classification table?
        logging.warning("Read classification suppressed, not writing classification table")
    if suppressions['suppress_plots']: # Suppressed plots?
        logging.warning("Plots suppressed, not creating plots")
    if args['xkcd']:
        # plots._XKCD = True
        pass
    #   Tell worker processes to ignore SIGINT (^C)
    #   by turning INTERUPT signals into IGNORED signals
    sigint_handler = signal.signal(signal.SIGINT, signal.SIG_IGN)
    signal.signal(signal.SIGINT, sigint_handler)
    #   Read in reference and template sequences
    reference = toolkit.load_seq(args['reference'])
    template = toolkit.load_seq(args['template'])
    sys.exit((reference, template))
    #   Setup our multiprocessing pool
    #   Allow the user to specify the number of jobs to run at once
    #   If not specified, let multiprocessing figure it out
    try:
        pool = Pool(processes=args['num_cores']) # type: multiprocessing.Pool
    except KeyError:
        pool = Pool() # type: multiprocessing.Pool
    if all(map(lambda i: i > 1, (len(fastq_list), getattr(pool, '_processes')))):
        pass
    #   Close logfile
    args['logfile'].close()


if __name__ == '__main__':
    main()
