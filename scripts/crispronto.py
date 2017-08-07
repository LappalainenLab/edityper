#!/usr/bin/env python

from __future__ import print_function
from __future__ import division

import sys

try:
    from arguments import make_argument_parser
except ImportError:
    sys.exit("FAIL")

# try:
#     import argparse
# except ImportError:
#     sys.exit("Please use Python 2.7 or 3.X")


# try:
#     range = xrange
# except NameError:
#     pass


# try:
#     import ConfigParser as configparser
# except ImportError:
#     import configparser

def main():
    """CRISPronto"""
    parser = make_argument_parser()
    if not sys.argv[1:]:
        sys.exit(parser.print_help())
    args = {key: value for key, value in vars(parser.parse_args()).items() if value}
    sys.exit(args)


if __name__ == '__main__':
    main()
