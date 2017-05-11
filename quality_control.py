#!/usr/bin/env python2

'''Quality control for the CRISPR program'''

from __future__ import print_function

import sys
if sys.version_info.major is not 2 and sys.version_info.minor is not 7:
    sys.exit("Please use Python 2.7 for this module: " + __name__)


import time
import random
import logging

from itertools import repeat

try:
    import numpy as np
    from scipy.stats import norm
except ImportError as error:
    sys.exit("Please install Numpy and SciPy for this module")

try:
    import genetic_toolpack as toolpack
    import NW_py as NW
except ImportError as error:
    sys.exit("Please keep this program in it's directory to load custom modules: " + error.message)


random.seed(a=time.time())


def align_reference(reference, template, gap_penalty): # type: (str, str, int) -> (str, str)
    """Align our template to our reference"""
    logging.info("Aligning reference and template sequences")
    align_start = time.time()
    fwd_ref, fwd_template, fwd_qual_score = NW.align_glocal(reference, template, gap_penalty)
    rev_ref, rev_template, rev_qual_score = NW.align_glocal(reference, toolpack.rvcomplement(seq=template), gap_penalty)
    logging.debug("Alignment took %s seconds", round(time.time() - align_start, 3))
    if rev_qual_score > fwd_qual_score:
        logging.warning("Using reverse alignment as it had higher quality")
        return rev_ref, rev_template
    else:
        logging.warning("Using forward alignment as it had higher quality")
        return fwd_ref, fwd_template


def validate_reference_alignment(reference, template, snp_mode=False): # type: (str, str, bool) -> List[List[int, Tuple[str]]]
    """Validate our reference/template alignment and find mismatches"""
    logging.info("Validating reference/template alignment")
    validate_start = time.time()
    ref_template_mismatch = toolpack.get_mismatch(seq_a=reference, seq_b=template) # type: List[List[int, Tuple[str]]]
    if not ref_template_mismatch:
        logging.error("There must be at least one mismatch between the template and reference sequences")
        sys.exit(1)
    if len(ref_template_mismatch) > 1 and snp_mode:
        sys.exit(logging.error("There can only be one mismatch between the template and reference sequences in SNP mode"))
    logging.debug("Validation took %s seconds", round(time.time() - validate_start, 3))
    return ref_template_mismatch


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
        fastq, # type: toolpack.FastQ
        raw_sequences, # type: List[str]
        reference, # type: str
        gap_open, # type: int
        gap_extension, # type: int
        pvalue_threshold # type: float
):
    # type: (...) -> (bool, float)
    """Determine if we're aligning our reads in the forward or reverse direction"""
    ten_percent = int(round(0.1 * len(raw_sequences)) + 1)
    sampled_reads = random.sample(raw_sequences, k=min([500, ten_percent])) # Sample at most 500 reads
    #   Set up our zip objects
    perm_reads = map(lambda read: ''.join(random.sample(read, k=len(read))), sampled_reads) # Permutate our reads
    ref_rc = toolpack.rvcomplement(reference) # Reverse complement our reference sequence
    norm_f = zip(repeat(reference), sampled_reads, repeat(gap_open), repeat(gap_extension)) # type: List
    perm_f = zip(repeat(reference), perm_reads, repeat(gap_open), repeat(gap_extension)) # type: List
    norm_r = zip(repeat(ref_rc), sampled_reads, repeat(gap_open), repeat(gap_extension)) # type: List
    perm_r = zip(repeat(ref_rc), perm_reads, repeat(gap_open), repeat(gap_extension)) # type: List
    #   Get our alignment scores
    align = lambda tup: NW.align_aff(tup[0], tup[1], tup[2], tup[3]) # type: function
    norm_scores = map(align, norm_f) # type: List[Tuple]
    perm_scores = map(align, perm_f) # type: List
    norm_scores_r = map(align, norm_r) # type: List
    perm_scores_r = map(align, perm_r) # type: List
    #   Analyze our alignment scores
    get_third = lambda iterable: iterable[2] # type: function
    fwd_median = np.median(map(get_third, norm_scores)) # type: numpy.float64
    rev_median = np.median(map(get_third, norm_scores_r)) # type: numpy.float64
    if rev_median > fwd_median:
        reverse = True # type: bool
        use_scores = map(get_third, perm_scores_r) # type: List[int]
    else:
        reverse = False # type: bool
        use_scores = map(get_third, perm_scores) # type: List[int]
    msg = "Aligning in the %s direction" % ('reverse' if reverse else 'forward') # type: str
    score_threshold = np.std(use_scores) * norm.pdf(1 - pvalue_threshold) + np.median(use_scores) # type: numpy.float64
    logging.warning("FASTQ %s: %s", str(fastq), msg)
    logging.warning("FASTQ %s: %s vs %s (norm vs reverse) - threshold: %s", str(fastq), fwd_median, rev_median, score_threshold)
    return reverse, fwd_median, rev_median, score_threshold
