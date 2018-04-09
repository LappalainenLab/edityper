#!/usr/bin/env python

"""Quality control for EdiTyper"""

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
    raise SystemExit(error)


try:
    if PYTHON_VERSION == 3:
        from edityper import toolkit
        from edityper import recnw
        from edityper.analysis import NA
    elif PYTHON_VERSION == 2:
        import toolkit
        import recnw
        from analysis import NA
    else:
        raise SystemExit("Please use Python 2 or 3 for this module: " + __name__)
except ImportError as error:
    raise SystemExit(error)


random.seed(a=time.time())

def align_reference(reference, template, gap_penalty): # type: (str, str, int) -> (str, str)
    """Align our template to our reference"""
    logging.info("Aligning reference and template sequences")
    align_start = time.time() # type: float
    fwd_ref, fwd_template, fwd_qual_score = recnw.nw_lin( # type: str, str, int
        reference,
        template,
        gap_penalty=gap_penalty
    )
    rev_ref, rev_template, rev_qual_score = recnw.nw_lin( # type: str, str, int
        reference,
        toolkit.reverse_complement(sequence=template),
        gap_penalty=gap_penalty
    )
    logging.debug("Alignment took %s seconds", round(time.time() - align_start, 3))
    if rev_qual_score > fwd_qual_score:
        logging.warning("Using reverse alignment as it had higher quality")
        return rev_ref, rev_template
    else:
        logging.warning("Using forward alignment as it had higher quality")
        return fwd_ref, fwd_template


def get_snp_states(reference, template, mismatch): # type: (str, str, List[int, Tuple[str]]) -> (int, str, str)
    """Get the SNP states from our alignment"""
    logging.info("Finding reference and template SNP states")
    snp_start = time.time() # type: float
    snp_index = mismatch[0] # type: int
    try:
        reference_state = reference[snp_index] # type: str
    except (TypeError, IndexError):
        reference = NA
    try:
        target_snp = template[snp_index] # type: str
    except (TypeError, IndexError):
        target_snp = NA
    logging.debug("Finding reference and template SNP states took %s seconds", round(time.time() - snp_start, 3))
    return snp_index, reference_state, target_snp


def determine_alignment_direction(
        fastq_name, # type: str
        unique_reads, # type: Iterable[str]
        reference, # type: str,
        gap_open, # type: int
        gap_extension, # type: int
        pvalue_threshold, # type: float
        args_dict # type: Dict[str, Any]
):
    # type: (...) -> (bool, float, float, float)
    """Deterime if we're aligning our reads in the forward or reverse direction"""
    logging.info("FASTQ %s: determining alignment direction", fastq_name)
    determine_start = time.time() # type: float
    ten_percent = int(round(0.1 * len(unique_reads)) + 1) # type: int
    sampled_reads = random.sample(unique_reads, k=min((500, ten_percent))) # Sample at most 500 reads
    permutate = lambda read: ''.join(random.sample(read, k=len(read))) # type: function
    ref_rc = toolkit.reverse_complement(sequence=reference) # type: str
    norm_scores, perm_scores, rev_scores, perm_rev = list(), list(), list(), list() # type: List[float], List[float], List[float], List[float]
    for read in sampled_reads: # type: str
        perm_read = permutate(read) # type: str
        _, _, score = recnw.nw_aff(reference, read, gap_op=gap_open, gap_ext=gap_extension) # type: _, _, float
        _, _, rev_score = recnw.nw_aff(ref_rc, read, gap_op=gap_open, gap_ext=gap_extension) # type: _, _, float
        _, _, perm_score = recnw.nw_aff(reference, perm_read, gap_op=gap_open, gap_ext=gap_extension) # type: _, _, float
        _, _, rev_perm = recnw.nw_aff(ref_rc, perm_read, gap_op=gap_open, gap_ext=gap_extension) # type: _, _, float
        norm_scores.append(score)
        perm_scores.append(perm_score)
        rev_scores.append(rev_score)
        perm_rev.append(rev_perm)
    norm_median = numpy.median(norm_scores) # type: float
    rev_median = numpy.median(rev_scores) # type: float
    do_reverse = rev_median > norm_median # type: bool
    try:
        threshold = args_dict['threshold']
    except KeyError:
        if do_reverse:
            threshold = numpy.std(perm_rev) * norm.pdf(1 - pvalue_threshold) + numpy.median(perm_rev) # type: float
        else:
            threshold = numpy.std(perm_scores) * norm.pdf(1 - pvalue_threshold) + numpy.median(perm_scores) # type: float
    msg = 'Aligning in the %s direction' % ('reverse' if do_reverse else 'forward') # type: str
    logging.warning("FASTQ %s: %s", fastq_name, msg)
    logging.warning("FASTQ %s: %s vs %s (norm vs reverse) - threshold: %s", fastq_name, norm_median, rev_median, threshold)
    logging.debug("FASTQ %s: Determinging alignment direction took %s seconds", fastq_name, round(time.time() - determine_start, 3))
    return do_reverse, norm_median, rev_median, threshold
