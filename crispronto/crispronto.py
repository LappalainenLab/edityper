#!/usr/bin/env python

from __future__ import print_function
from __future__ import division

import sys
PYTHON_VERSION = sys.version_info.major

try:
    if PYTHON_VERSION is 3:
        from crispronto.arguments import make_argument_parser
    elif PYTHON_VERSION is 2:
        from arguments import make_argument_parser
    else:
        sys.exit('WTF MATE')
except ImportError:
    sys.exit("FAIL")


def main():
    """CRISPronto"""
    parser = make_argument_parser()
    if not sys.argv[1:]:
        sys.exit(parser.print_help())
    args = {key: value for key, value in vars(parser.parse_args()).items() if value}
    sys.exit(args)


if __name__ == '__main__':
    main()
