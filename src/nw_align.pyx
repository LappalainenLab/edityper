from edityper import toolkit
from libcpp.string cimport string

cdef extern from "py_align.h":
	string nw_align(string seq_1, string seq_2, int gap_penalty)
	string nw_free_trail(string seq_1, string seq_2)
	string nw_align_aff(string seq_1, string seq_2, int gap_op, int gap_ext)
	string nw_align_aff_param(string seq_1, string seq_2, int gap_op, int gap_ext, int match, int mismatch)
	string nw_align_mem(string seq_1, string seq_2, int gap_penalty, int sim, int terminate)
	string nw_align_aff_mem(string seq_1, string seq_2, int gap_op, int gap_ext, int sim, int terminate)

def align_glocal(seq_1, seq_2, gap_penalty):
	result = nw_align(seq_1, seq_2, gap_penalty)
	al1, al2, score = result.split('|')
	score = int(score)
	return al1, al2, score


def align_freetrail(seq_1, seq_2):
	result = nw_free_trail(seq_1, seq_2)
	al1, al2, score = result.split('|')
	score = int(score)
	return al1, al2, score


def align_mem(seq_1, seq_2, gap_penalty, sim, terminate):
	result = nw_align_mem(seq_1, seq_2, gap_penalty, sim, terminate)
	al1, al2, score = result.split('|')
	score = int(score)
	return al1, al2, score


def align_aff(seq_1, seq_2, gap_op, gap_ext):
	result = nw_align_aff(seq_1, seq_2, gap_op, gap_ext)
	al1, al2, score = result.split('|')
	score = int(score)
	return al1, al2, score


def align_aff_param(seq_1, seq_2, gap_op, gap_ext, match, mismatch):
	result = nw_align_aff_param(seq_1, seq_2, gap_op, gap_ext, match, mismatch)
	al1, al2, score = result.split('|')
	score = int(score)
	return al1, al2, score


def align_aff_mem(seq_1, seq_2, gap_op, gap_ext, sim, terminate):
	result = nw_align_aff_mem(seq_1, seq_2, gap_op, gap_ext, sim, terminate)
	al1, al2, score = result.split('|')
	score = int(score)
	return al1, al2, score
