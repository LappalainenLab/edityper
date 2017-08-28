#!/usr/bin/env python

"""Quality control for CRIPRonto"""

from __future__ import division
from __future__ import print_function

import sys
PYTHON_VERSION = sys.version_info.major

import time
import random
import logging

try:
    import numpy
    from scipy.stats import norm
except ImportError as error:
    sys.exit(error)


try:
    if PYTHON_VERSION is 3:
        from crispronto import toolkit
        from crispronto import nw_align
    elif PYTHON_VERSION is 2:
        import toolkit
        import nw_align
    else:
        sys.exit("Please use Python 2 or 3 for this module: " + __name__)
except ImportError as error:
    sys.exit(error)


random.seed(a=time.time())

def align_reference(reference, template, gap_penalty): # type: (str, str, int) -> (str, str)
    """Align our template to our reference"""
    logging.info("Aligning reference and template sequences")
    align_start = time.time()
    fwd_ref, fwd_template, fwd_qual_score = nw_align.align_glocal(
        seq_1=reference,
        seq_2=template,
        gap_penalty=gap_penalty
    )
    rev_ref, rev_template, rev_qual_score = nw_align.align_glocal(
        seq_1=reference,
        seq_2=toolkit.reverse_complement(sequence=template),
        gap_penalty=gap_penalty
    )
    logging.debug("Alignment took %s seconds", round(time.time() - align_start, 3))
    if rev_qual_score > fwd_qual_score:
        logging.warning("Using reverse alignment as it had higher quality")
        return rev_ref, rev_template
    else:
        logging.warning("Using forward alignment as it had higher quality")
        return fwd_ref, fwd_template


def get_snp_states(reference, template, mismatch, mode): # type: (str, str, List, str) -> (int, str, str)
    """Get the SNP states from our alignment"""
    logging.info("Finding reference and template SNP states")
    snp_start = time.time()
    if mode.split('+')[0] == 'SNP':
        snp_index = mismatch[0][0]
    else:
        snp_index = mismatch[1][0]
    reference_state = reference[snp_index]
    target_snp = template[snp_index]
    logging.debug("Finding reference and template SNP states took %s seconds", round(time.time() - snp_start, 3))
    return snp_index, reference_state, target_snp


def determine_alignment_direction(
        fastq_name, # type: str
        unique_reads, # type: Iterable[str]
        reference, # type: str,
        gap_open, # type: int
        gap_extension, # type: int
        pvalue_threshold # type: float
):
    # type: (...) -> None
    """Deterime if we're aligning our reads in the forward or reverse direction"""
    logging.info("FASTQ %s: determining alignment direction", fastq_name)
    determine_start = time.time()
    ten_percent = int(round(0.1 * len(unique_reads)) + 1)
    sampled_reads = random.sample(unique_reads, k=min((500, ten_percent))) # Sample at most 500 reads
    permutate = lambda read: ''.join(random.sample(read, k=len(read)))
    ref_rc = toolkit.reverse_complement(sequence=reference)
    norm_scores, perm_scores, rev_scores, perm_rev = list(), list(), list(), list()
    for read in sampled_reads:
        perm_read = permutate(read)
        _, _, score = nw_align.align_aff(seq_1=reference, seq_2=read, gap_op=gap_open, gap_ext=gap_extension)
        _, _, rev_score = nw_align.align_aff(seq_1=ref_rc, seq_2=read, gap_op=gap_open, gap_ext=gap_extension)
        _, _, perm_score = nw_align.align_aff(seq_1=reference, seq_2=perm_read, gap_op=gap_open, gap_ext=gap_extension)
        _, _, rev_perm = nw_align.align_aff(seq_1=ref_rc, seq_2=perm_read, gap_op=gap_open, gap_ext=gap_extension)
        norm_scores.append(score)
        perm_scores.append(perm_score)
        rev_scores.append(rev_score)
        perm_rev.append(rev_perm)
    norm_median = numpy.median(norm_scores)
    rev_median = numpy.median(rev_scores)
    do_reverse = rev_median > norm_median
    if do_reverse:
        threshold = numpy.std(perm_rev) * norm.pdf(1 - pvalue_threshold) + numpy.median(perm_rev)
    else:
        threshold = numpy.std(perm_scores) * norm.pdf(1 - pvalue_threshold) + numpy.median(perm_scores)
    msg = 'Aligning in the %s direction' % ('reverse' if do_reverse else 'forward')
    logging.warning("FASTQ %s: %s", fastq_name, msg)
    logging.warning("FASTQ %s: %s vs %s (norm vs reverse) - threshold: %s", fastq_name, norm_median, rev_median, threshold)
    logging.debug("FASTQ %s: Determinging alignment direction took %s seconds", fastq_name, round(time.time() - determine_start, 3))
    return do_reverse, norm_median, rev_median, threshold
