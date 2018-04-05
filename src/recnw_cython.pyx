from libcpp.string cimport string
from libcpp cimport bool

cdef extern from "recnw.h":
	string recnw_affine(string seq_1, string seq_2, float gap_op, float gap_ext, float match, float mismatch, bool free_hgap_1, bool free_hgap_2, bool free_tgap_1, bool free_tgap_2, int sim, int terminate)
	string recnw_reg(string seq_1, string seq_2, float gap_penalty, float match, float mismatch, bool free_hgap_1, bool free_hgap_2, bool free_tgap_1, bool free_tgap_2, int sim, int terminate)


def nw_aff(ref, seq, gap_op=8, gap_ext=1, match=5, mismatch=-4, head_free=(False, True), tail_free=(False, True), sim=-1, terminate=1):
	'''RecNW with affine gap penalty'''
	result = recnw_affine(ref, seq, gap_op, gap_ext, match, mismatch, head_free[0], head_free[1], tail_free[0], tail_free[1], sim, terminate)
	al1, al2, score = result.split('|')
	score = float(score)
	return al1, al2, score


def nw_lin(ref, seq, gap_penalty=8, match=5, mismatch=-4, head_free=(False, True), tail_free=(False, True), sim=-1, terminate=1):
	'''RecNW with linear gap penalty'''
	result = recnw_reg(ref, seq, gap_penalty, match, mismatch, head_free[0], head_free[1], tail_free[0], tail_free[1], sim, terminate)
	al1, al2, score = result.split('|')
	score = float(score)
	return al1, al2, score