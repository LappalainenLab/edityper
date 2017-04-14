#!/usr/bin/env python2

from __future__ import print_function

import sys
if sys.version_info.major is not 2 and sys.version_info.minor is not 7:
    sys.exit("Please use Python 2.7 for this module: " + __name__)


import itertools
from copy import copy
from write_sam import SAM, make_cigar

try:
    import regex
except ImportError as error:
    sys.exit("Please install %s for this module" % error.name)


def calc_read_pos(alignment): # type: alignment.Alignment -> int
    """Calculate alignment read position"""
    cigar = make_cigar(alignment=alignment, exact=True) # type: str
    reference = alignment.get_aligned_reference() # type: str
    read = alignment.get_aligned_read() # type: str
    if len(reference) != len(read):
        raise ValueError("The read and reference must be the same length")
    reference = reference.replace('-', '') # type: str
    read = read.replace('-', '') # type: str
    SAM.validate_cigar(cigar=cigar, sequence=read)
    operations = regex.findall(r'(\d+[%s])' % ''.join(SAM.VALID_OPS), cigar) # type Tuple[str]
    digits = regex.compile(r'(\d+)[%s]' % ''.join(SAM.VALID_OPS)).findall
    pattern = '(' # type: str
    possible_references = itertools.repeat(reference, 2) # type: itertools.repeat
    num_ref_adjust = 0 # type: int
    for oper in operations: # type: str
        this_op = ''.join((char for char in oper if not char.isdigit())) # type: str
        this_digit = int(digits(oper)[0]) # type: int
        if this_op in {'='}:
            pattern += read[:this_digit]
            read = read[this_digit:] # type: str
        elif this_op in {'X'}:
            i = 0
            while i < this_digit:
                pattern += '[^%s]' % read[0]
                read = read[1:] # type: str
                i += 1
        elif this_op in {'I', 'S', 'P'}: # Insertion in read relative to reference, soft clipping of read, or silent deletion of reference
            read = read[this_digit:] # type: str
            continue
        elif this_op in {'D'}: # Deletion in read relative to reference
            adjustments = list(regex.finditer(r'%s.{%s})' % (pattern, this_digit), ref) for ref in possible_references) # type: List
            adjustments = (match.groups() for adj in adjustments for match in adj) # type: generator
            adjustments = itertools.chain.from_iterable(adjustments) # type: itertools.chain
            adj_ref = copy(possible_references) # type: List
            for adj in adjustments:
                for ref in adj_ref: # type: str
                    temp = list() # type: list
                    try:
                        ref.index(adj)
                        ref.replace(adj, adj[:-this_digit])
                        temp.append(ref)
                    except ValueError:
                        continue
                adj_ref = copy(temp) # type: list
            possible_references = copy(adj_ref) # type: list
            num_ref_adjust += this_digit
        else:
            pass
        possible_references = (regex.findall(r'%s.*)' % pattern, ref) for ref in possible_references) # type: generator
        possible_references = tuple(itertools.chain.from_iterable(possible_references)) # type: Tuple[str]
        if len(possible_references) == 1:
            break
        elif len(possible_references) == 0:
            raise ValueError("Something happened")
    adjusted_ref = possible_references[0] # type: str
    return regex.search(r'(%s){e<=%d}' % (adjusted_ref, num_ref_adjust), reference, regex.BESTMATCH).start()


test_ref = 'TGAACAGGCCTGGAAGCCCAACCCGGTGCGGACCGCCCTCAGTCTCCTAGTCCGGGCGCGCTGAGCGCAGGGGGGCGCTGTCCTCACCCACGCG--TAAGTCGCCGCAGG--GCGTTTCCGGGCG-GAGAAAACCTAC---ACGTGATGGGCG-----CCCACCGAGTGCCAGCACCGCCTCCTCCAGCTCCGCCAAGTAGGTGGGATCCACTACTTTGCAGAGGAGGAAGCCGTTCAGGCCTCCCTCTGTGCACTTTGCAAAAGC-CTCTGCC--CTTCCATCTCGAA--TCCCTTGGCGCCGATCACACTTCCTCTGCCTGGAAACCTGGAGCCGTCTCTCGCGAGACGTCCTCCGCCCTGTAGAAGGCCGTTTCGGTTCTTCGTGCGCGGTAGCCGCCCCACTTGCGGGATTCCAAGGCCTCATCGAGTGCGGGTATCTGGCTGTGGATTCGCCGCCGTCCTGCTGGACGCCTGGAGGCTCGAACCCCGCCGCCCCCCTACCCCAGGCTTTACTCCCACCCCGGCTTCCGCCCACTGTGCTGCCCTTCCTCGGACCTGGGCTGTCGGGAGAGCTGGAGGTGAGCGCTCTTGGGAGAGCCCAGCCAATTTAAAGAG'

test_read = '----------------------------------------------------------------------------------------------ATTAGGT--------GTCG-GTCTATGGGCGTGAG-CTA-CTCCGGAACCTTA-GGGCGTTCACCCCGCCGA-TG--G----CGCGTGCT--------------TC--TTGG--------CTTTGC-----CGGAAGATG-TCCCGCCT--C-C--TGCA---G---AGCGCTCTCTGCCGAGGTCCA-----AAGGT-CC--GTCGCC-TTCACACTA-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------'

from alignment import Alignment
a = Alignment(ref_align=test_ref, read_align=test_read, score=3)
