# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

#AUTHOR: ALEXANDRE YAHI
#AFFILATION: LAPPALAINEN LAB
#DATE: 2015-2017
#VERSION: 0.6


'''CRISPR gene editing results analysis'''

import os
from os import listdir
from os.path import isfile, join, isdir
import re
import sys
import csv
import gzip
# Data structure packages
from collections import defaultdict
from collections import Counter
from operator import itemgetter
# TIME MANAGEMENT packages
import time
import random
# MATH packages
import numpy as np
from scipy.stats import norm
# ALIGNMENT package
from genetic_toolpack import load_fastq, load_seq, get_mismatch, side_trimmer, rvcomplement, trim_interval, find_gaps, find_gaps_kn, get_mismatch_match, sim_seq
#Multiprocessing
import multiprocessing as mp
#C++ ALIGNMENT FUNCTION
from NW_py import align_glocal, align_freetrail, align_mem, align_aff, align_aff_mem
# PLOTING packages
from genetic_plotpack import violin_quality, locus_map, violin_matplot


# RELATIVE PATH
script_name =  os.path.basename(__file__)
abs_path = os.path.realpath(__file__).split(script_name)[0]

# INPUT ARGUMENTS
inst = "USAGE : settings=setting_file.txt"
try:
	args = dict([argv.split('=') for argv in sys.argv[1:]])
	if len(args.keys()) == 1:
		raise ValueError
except:
	print inst


def display(text):
	sys.stdout.write("\r %s" % text)
	sys.stdout.flush()


def get_settings(input_path, setting_file):
	"""Read setting file specified"""
	f_input = open(setting_file,'rb')

	input_name = ''			# Default processes all fastq files in input/ directory
	reference = ''			# Reference - Default will infer amplicon from file name with "amplicon-X_R1_blah.fastq" template
	template = ''           # Template - Default will copy the name of the amplicon
	gap_opening = 8			# Default gap opening is 8
	gap_extension = 1 		# Default gap extension is 1
	analysisMode = 'SNP'    # Default unique SNP mode (supported: 'SNP+PAM')
	mode = 'uniq'			# Determined by the input_name, (default, "uniq")

	dir_list = [ f for f in listdir(input_path) if isdir(join(input_path,f))]
	file_list = [ f for f in listdir(input_path) if isfile(join(input_path,f)) and '.fastq' in f]


	settings_dict = {'analysisMode':analysisMode,'mode':mode, 'input_name':input_name, 'reference':reference, 'template':template, 'gap_opening':gap_opening, 'gap_extension':gap_extension}

	for line in f_input:
		if '#' == line[0]:
			continue
		else:
			line = line.replace(' ', '')
			line = line.replace('\n', '')
			parameter, value = line.split('=')
			print parameter, ' - ', value
			if value == '':
				continue
			else:
				settings_dict[parameter] = value

	# Mode: uniq/pair/list/batch - depend on the input_name
	if '.fastq' in settings_dict['input_name'] and settings_dict['input_name'] in file_list:
		settings_dict['mode'] = 'uniq'
	elif ';' in settings_dict['input_name']:
		settings_dict['mode'] = 'list'
	elif settings_dict['input_name'] in dir_list:
		settings_dict['mode'] = 'batch'
	else:
		settings_dict['mode'] = 'pair'

	return settings_dict


def CRISPR_analysis(analysisMode, inputFile, inputReference, inputTemplate, fileRadical, outputPath, gap_open, gap_ext):
	'''Main function to align and analyze the data'''
	pvalTHR = 1E-3
	gap_open = int(gap_open)
	gap_ext = int(gap_ext)
	gap_penalty = gap_open

########	 LOADING DATA  ########
	display('LOADING DATA...')
	load_start = time.time()

	loaded_fastq = load_fastq(inputFile)			#LOAD THE FASTQ FILE - Header / Read / Quality
	ref = load_seq(inputReference)					#LOAD THE REFERENCE SEQUENCE
	template = load_seq(inputTemplate)				#LOAD THE TEMPLATE SEQUENCE

	display('LOADING DATA...OK\n')
	load_time = round(time.time() - load_start,2)

########################  QUALITY CONTROL  ########################
	display('QUALITY CONTROL...')
	#Find how to match the reference and the template
	al_ref_n, al_temp_n, q_score_n = align_glocal(ref, template, gap_penalty)
	al_ref_rv, al_temp_rv, q_score_rv = align_glocal(ref, rvcomplement(template), gap_penalty)
	if q_score_rv > q_score_n:
		al_ref = al_ref_rv
		al_temp = al_temp_rv
	else:
		al_ref = al_ref_n
		al_temp = al_temp_n

	if '-' in al_ref:
		raise AttributeError('Insertions in the template - wrong reference')
	if '-' in side_trimmer(al_temp):
		raise AttributeError('Deletions in the template - wrong ')

	temp2ref_mismatch = get_mismatch(al_ref,al_temp) #RETURN A LIST [INDEX, MISMATCH OF THE TEMPLATE]
	if len(temp2ref_mismatch) > 1 and analysisMode == 'SNP':
		raise AttributeError('More than one mismatch - Not supported!')
	if temp2ref_mismatch == []:
		raise AttributeError('No differences between reference and template - Not supported!')


	if analysisMode == 'SNP':					# analysisMode=SNP
		snp_index = temp2ref_mismatch[0][0]
	else:
		if analysisMode.split('+')[0] == 'SNP': # analysisMode=SNP+PAM
			snp_index = temp2ref_mismatch[0][0]
		else: 									# analysisMode=PAM+SNP
			snp_index = temp2ref_mismatch[1][0]

	snp_ref = al_ref[snp_index]						# Nucleotide of reference
	snp_target = al_temp[snp_index]					# Wanted mutation

	display('QUALITY CONTROL...OK\n')

########################  ISOLATE THE READS  ########################
	display('READS PRE-PROCESSING...')

	reads_raw = [elem[1] for elem in loaded_fastq]	# Isolate the raw reads in one list
	num_reads = len(reads_raw) 						# Total number of raw reads

	if num_reads == 0:
		return ''
	# LOOKING FOR READS DUPLICATES

	reads_dict = defaultdict(int)					# Dictionary of reads -> [read] = count
	for seq in reads_raw:
		reads_dict[seq] += 1

	reads_uniq = reads_dict.keys() 					# List of uniques read
	uniq_perc = round(float(len(reads_uniq)/float(num_reads))*100,2)

	display('READS PRE-PROCESSING... %s reads - %s unique reads (%s %%)\n' % (num_reads, len(reads_uniq),uniq_perc))
########################  FIND HOW TO ALIGN THE READS  ########################
	display('FIND HOW TO ALIGN...')

	score_normal = list()
	perm_score_normal = list()
	score_reverse = list()
	perm_score_reverse = list()

	random.shuffle(reads_raw)

	to_sample = int(round(0.1*num_reads,0)+1)

	if to_sample < 500:
		sampled_reads = reads_raw[:to_sample]
	else:
		sampled_reads = reads_raw[:500]					# Randomly sample 500 reads

	for s_read in sampled_reads:
		_, _, score = align_aff(ref, s_read, gap_open, gap_ext)
		_, _, perm_score = align_aff(ref, ''.join(random.sample(s_read,len(s_read))), gap_open, gap_ext)
		score_normal.append(score)
		perm_score_normal.append(perm_score)

	reverse_ref = rvcomplement(ref)

	for s_read in sampled_reads:
		_, _, score = align_aff(reverse_ref, s_read, gap_open, gap_ext)
		_, _, perm_score = align_aff(reverse_ref,''.join(random.sample(s_read,len(s_read))), gap_open, gap_ext)
		score_reverse.append(score)
		perm_score_reverse.append(perm_score)

	norm_median = np.median(score_normal)
	rev_median = np.median(score_reverse)

	reverse_status = 'Normal'
	if norm_median <= rev_median:
		reverse = True
		reverse_status = 'Reverse'
		ScoreTHR = np.std(perm_score_reverse)*norm.ppf(1-pvalTHR) + np.median(perm_score_reverse)
	else:
		reverse = False
		ScoreTHR = np.std(perm_score_normal)*norm.ppf(1-pvalTHR) + np.median(perm_score_normal)

	ScoreTHR = round(ScoreTHR,3)


	if reverse:
		# Reverse complement safe update of data structures
		proc_check = dict()
		for i, old_read in enumerate(reads_uniq):
			# print reads_dict[old_read] ## CHECK
			new_read = rvcomplement(old_read)
			proc_check[new_read] = reads_dict[old_read]
			# reads_uniq[i] = new_read
			# reads_dict[new_read] = reads_dict.pop(old_read)
			# print reads_dict[new_read] ## CHECK
		print 'After reverse...'
		reads_dict = dict(proc_check)
		reads_uniq = reads_dict.keys()


	display('FIND HOW TO ALIGN... %s - %s vs %s - thres. %s \n' % (reverse_status, norm_median, rev_median, ScoreTHR))
########################  ALIGNMENT  ########################
	display('ALIGNING READS...')

	start_time = time.time()

	Scores_li = defaultdict(list)								# Score dict - by length
	Align_li = defaultdict(list)								# dict of aligned outputs list - by length

	#Check if reads have unique length:
	len2uniq_temp = defaultdict(list)
	for u_read in reads_uniq:
		length = len(u_read)
		len2uniq_temp[length].append(u_read)

	# Sort the unique reads - doesn't matter if several lengths
	len2uniq = dict()
	for length in len2uniq_temp.keys():
		read_batch = len2uniq_temp[length]
		read_batch.sort()							# Sorting the unique reads
		len2uniq[length] = read_batch

	memory = True #THE RECURSIVE ALIGNMENT

	# Keep track of matrix re-used lines
	reuse_lines = int()

	if memory:
		for length in len2uniq.keys():
			count = int()
			temp = str()
			total = len(len2uniq[length])

			for sample in len2uniq[length]:
				count += 1
				if temp == str():
					al1, al2, score = align_aff_mem(ref, sample, gap_open, gap_ext, 0, 0)

					al1_temp = al1
					al2_temp = al2
					temp = sample
					Align_li[length].append([al1,al2])
					Scores_li[length].append(score)
					continue

				index = sim_seq(temp,sample)
				reuse_lines += index

				if count == total:
					al1, al2, score = align_aff_mem(ref, sample, gap_open, gap_ext, index, 1) # de-allocate memory - delete the alignment matrix
				else:
					al1, al2, score = align_aff_mem(ref, sample, gap_open, gap_ext, index, 0)

				# Store the aligned sequences and score
				Align_li[length].append([al1,al2])
				Scores_li[length].append(score)

				temp = sample

	else:	#Not keeping memory
		for length in len2uniq.keys():
			for sample in len2uniq[length]:
				al1, al2, score = align_aff(ref, sample, gap_open, gap_ext)
				Align_li[length].append([al1,al2])
				Scores_li[length].append(score)

	align_time = round(time.time() - start_time,2)

	display('ALIGNING READS... OK - %s seconds\n' % (align_time))
	######################################################
	display('ANALYSING RESULTS...')
	analysis_start = time.time()

	## OUTPUTS ##
	# Counting by index position
	deletion_index = [ list() for _ in xrange(len(ref))]
	insertion_index = [ list() for _ in xrange(len(ref))]
	mismatch_index = [ list() for _ in xrange(len(ref))]
	match_index = defaultdict(int) #Counter, key is ref position, count is coverage

	# Classifying reads
	class_reads = dict()
	#############

	for category in ['HDR', 'NHEJ', 'NO_EDIT', 'Disc']:
		class_reads[category] = defaultdict(list)

	out_idx = int() # Index to reference final results - different from index per length category
	super_check = int()
	for length_cat in len2uniq.keys():
		for idx,read_seq in enumerate(len2uniq[length_cat]):							# Number of total unique reads
			out_idx += 1
			tag = str()
			nm_del = int()
			nm_ins = int()
			nm_mis = int()
			multi = reads_dict[read_seq] # Get multiplicity of the read
			super_check += multi
			if Scores_li[length_cat][idx] < ScoreTHR:
				class_reads['Disc'][out_idx] = [multi, 0, 0, 0]
				continue

			#LOAD THE ALIGNMENTS
			al_ref = Align_li[length_cat][idx][0] #Aligned reference
			al_read = Align_li[length_cat][idx][1] #Aligned read

			#TRIM INTERVAL OF THE READ
			# use head-tail of trimmed aligned ref to trim al_read
			ref_head, ref_tail = trim_interval(al_ref)
			# Trimming aligned sequences
			al_ref = al_ref[ref_head:ref_tail]
			al_read = al_read[ref_head:ref_tail]
			# Get head and tail of the newly trimmed al_read
			head, tail = trim_interval(al_read)

			# #############
			# ## TEST ###
			# ##########

			# head_ref, tail_ref = trim_interval(al_ref)
			# al_ref = al_ref[head_ref:]
			# al_read = al_read[head_ref:]
			# #############
			# ## TEST END###
			# ##########

			#COUNT INSERTIONS and POSITIONS
			ins_gap_list = find_gaps(al_ref)

			if ins_gap_list == []:
				nm_ins = 0
			else:
				nm_ins = len(ins_gap_list)
				temp = int(0)	#Len of the previous insertion
				for gap in ins_gap_list:
					pos, gi_length = gap
					pos = pos - temp
					temp = temp + gi_length
					new_entry = [gi_length]*multi
					try:
						insertion_index[pos].extend(new_entry)
					except IndexError:
						print pos, temp

			#REMOVING INSERTIONS SEGMENTS
			if nm_ins != 0:
				new_ref = str()
				new_read = str()

				for i in xrange(len(al_ref)):
					if al_ref[i] == '-':
						continue
					else:
						new_ref = new_ref + al_ref[i]
						new_read = new_read + al_read[i]
				al_ref = new_ref
				al_read = new_read

			#COUNT DELETIONS
			del_gap_list = find_gaps_kn(al_read, head, tail) #With known head and tail before the insertions removal

			if del_gap_list == []:
				nm_del = 0
			else:
				nm_del = len(del_gap_list)
				for gap in del_gap_list:
					pos, gd_length = gap
					new_entry = [gd_length]*multi
					deletion_index[pos].extend(new_entry)

			#COUNT and INDEX MISMATCHES
			mis_list, match_list = get_mismatch_match(al_ref, al_read, head, tail)

			nm_mis = len(mis_list)

			for mismatch in mis_list:

				pos, perm = mismatch
				source, SNP = perm
				new_entry = [SNP]*multi
				try:#MODIF
					mismatch_index[pos].extend(new_entry)
				except IndexError:
					print 'length cat: ', length_cat
					print pos, perm
					print al_ref
					print al_read
					print head, tail
					print 'POSITION'

			# COVERAGE: count when the reference is matched

			for match_pos in match_list:
				match_index[match_pos] += multi

			#CLASSIFICATION

			if al_read[snp_index] == snp_target:
				tag = 'HDR'
			else:
				if nm_del == 0 and nm_ins == 0 and nm_mis == 0:
				#if al_read == al_ref:
					tag = 'NO_EDIT'
				else:
					tag = 'NHEJ'
			class_reads[tag][out_idx] = [multi, nm_ins, nm_del, nm_mis]	#Accounting for the multiplicity of the read

	print '\n'
	print 'SUPER CHECK: ', super_check
	counted_total = int()
	for tag in class_reads.keys():
		count = int()
		indels = list()
		for read_id in class_reads[tag].keys():
			multi, ins_n, del_n, mis_n = class_reads[tag][read_id]
			ind = ins_n + del_n
			indels.append(ind)
			count += multi
			counted_total += multi
		print tag, ' : ', count, ' avg indels: ', round(np.mean(indels),2)
	if counted_total == num_reads:
		print 'CHECKSUM: OK!'
	else:
		print 'CHECKSUM: ERROR! %s MISSING READS AFTER CLASSIFICATION' % (num_reads - counted_total)


	display('ANALYSING RESULTS... DONE!\n')
	analysis_time = round(time.time() - analysis_start,2)
	print 'TOTAL TIME FOR ANALYSIS... ', analysis_time, ' sec'


###########  REPORTING  #############
	display('REPORTING......')
	reporting_start = time.time()

	cumulDel = [ int() for _ in xrange(len(ref))]

	for i, dist in enumerate(deletion_index):
		for length in deletion_index[i]:
			for elem in range(length):
				cumulDel[i+elem] += 1

	## Saving unique reads: read_ID, multiplicity, reads
	#f = open(abs_path + 'output_analysis/%s_.txt' % fileRadical, 'wb')


	## Saving Ins, Del, Mism: Position, Length (or position, and #of nucleotides)
	f = open(outputPath + '%s_events.txt' % fileRadical, 'wb')
	f.write('Position\tReference\tCoverage\tDelStart\tDelLengthAvg\tDelNucCount\tInsStart\tInsLengthAvg\tA\tT\tC\tG\n')

	# NOTE: INSERT LINE IF '-' BEFORE AND AFTER THE SNP POSITION TO MAKE IT MORE LEGIBLE
	locus_ins_start = [ int() for _ in xrange(len(ref))]
	locus_mis_pos = [ int() for _ in xrange(len(ref))]


	for i in range(len(ref)):
		nuc_counter = Counter(mismatch_index[i])

		# For ploting later
		locus_mis_pos[i] = sum(nuc_counter.values())


		count_del = len(deletion_index[i])
		if count_del == 0:
			mean_del = 0
		else:
			mean_del = round(np.mean(deletion_index[i]),2)

		count_ins = len(insertion_index[i])
		if count_ins == 0:
			mean_ins = 0
		else:
			mean_ins = round(np.mean(insertion_index[i]),2)

		# For ploting later
		locus_ins_start[i] = count_ins

		try:
			nuc_counter[ref[i]] = match_index[i]
		except KeyError:
			nuc_counter[ref[i]] = 0

		# Compute coverage: match, mismatch, deletion
		covered = cumulDel[i]
		for cov_nuc in nuc_counter.keys():
			covered += nuc_counter[cov_nuc]


		if i == snp_index:
			f.write("--------------------------------------------------------------------\n")
			f.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (i+1,ref[i],covered,count_del,mean_del,cumulDel[i],count_ins,mean_ins,nuc_counter['A'],nuc_counter['T'],nuc_counter['C'],nuc_counter['G']))
			f.write("--------------------------------------------------------------------\n")
		else:
			f.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (i+1,ref[i],covered,count_del,mean_del,cumulDel[i],count_ins,mean_ins,nuc_counter['A'],nuc_counter['T'],nuc_counter['C'],nuc_counter['G']))



	f.close()


	# ## Saving read analysis: read_id, tag, #ins, #del, #mism.
	# detailed = open(outputPath + '%s_detailedCounts.txt' % fileRadical, 'wb')

	# detailed.close()
	## General report
	f = open(outputPath + '%s_general.txt' % fileRadical, 'wb')
	f.write("# SUMMARY OF THE ANALYSIS - CRISPR analysis v0.4 - Lappalainen Lab\n")
	f.write("# Template is different from Reference on position %s (%s -> %s) \n" % ((snp_index+1),snp_ref,snp_target))
	f.write("# Total number of reads: %s - %s unique reads (%s %%)\n" % (num_reads, len(reads_uniq),uniq_perc))
	f.write("# Normal of reverse complementary: %s - %s vs %s - threshold score to discard: %s \n" % (reverse_status, norm_median, rev_median, ScoreTHR))

	for tag in class_reads.keys():
		f.write("%s:\n" % tag)
		count = int()
		insertion = list()
		deletion = list()
		mismatch = list()

		no_indels = int()
		with_del = int()
		with_ins = int()
		with_indels = int()

		for read_id in class_reads[tag].keys():
			multi, ins_n, del_n, mis_n = class_reads[tag][read_id]
			insertion.extend([ins_n]*multi)
			deletion.extend([del_n]*multi)
			mismatch.extend([mis_n]*multi)
			count += multi
			if ins_n == 0 and del_n > 0:
				with_del += multi
			elif del_n == 0 and ins_n > 0:
				with_ins += multi
			elif ins_n == 0 and del_n == 0:
				no_indels += multi
			elif ins_n != 0 and del_n != 0:
				with_indels += multi
			else:
				print 'Not classified: ',tag, multi, ins_n, del_n, mis_n

		if tag == 'NO_EDIT':
			f.write("\tCount: %s - %s %% \n" % (count,round(count*100/float(num_reads),2)))
			NotEdited_count = count
		elif tag == 'Disc':	#Short version
			f.write("\tCount: %s - %s %% \n" % (count,round(count*100/float(num_reads),2)))
			DiscardedReads_count = count
		elif tag == 'HDR':
			f.write("\tCount: %s - %s %% \n" % (count,round(count*100/float(num_reads),2)))

			HDRnoIndels_count = no_indels
			HDRgap_count = count - no_indels

			if count != 0:
				f.write("\t HDR without indels: %s - %s %% \n" % (no_indels,round(no_indels*100/float(count),2)))
				f.write("\t HDR with deletions only: %s - %s %% \n" % (with_del,round(with_del*100/float(count),2)))
				f.write("\t HDR with insertions only: %s - %s %% \n" % (with_ins,round(with_ins*100/float(count),2)))
				f.write("\t HDR with indels: %s - %s %% \n" % (with_indels,round(with_indels*100/float(count),2)))
			else:
				f.write("\t HDR without indels: 0 - N/A %% \n")
				f.write("\t HDR without insertions only: 0 - N/A %% \n")
				f.write("\t HDR without deletions only: 0 - N/A %% \n")
				f.write("\t HDR with indels: 0 - N/A %% \n")

			f.write("\t Insertions events: %s (avg: %s | std: %s )\n" % (np.sum(insertion),round(np.mean(insertion),2),round(np.std(insertion),2)))
			f.write("\t Deletions events: %s (avg: %s | std: %s )\n" % (np.sum(deletion),round(np.mean(deletion),2),round(np.std(deletion),2)))
			f.write("\t Mismatches events: %s (avg: %s | std: %s )\n" % (np.sum(mismatch),round(np.mean(mismatch),2),round(np.std(mismatch),2)))
		else:
			NHEJ_count = count
			if count !=0:
				f.write("\tCount: %s - %s %% \n" % (count,round(count*100/float(num_reads),2)))
			else:
				f.write("\tCount: 0 - N/A %% \n")

			f.write("\t Insertions events: %s (avg: %s | std: %s )\n" % (np.sum(insertion),round(np.mean(insertion),2),round(np.std(insertion),2)))
			f.write("\t Deletions events: %s (avg: %s | std: %s )\n" % (np.sum(deletion),round(np.mean(deletion),2),round(np.std(deletion),2)))
			f.write("\t Mismatches events: %s (avg: %s | std: %s )\n" % (np.sum(mismatch),round(np.mean(mismatch),2),round(np.std(mismatch),2)))

	f.close()

	# PLOT LOCUS
	locus_map(cumulDel,locus_ins_start,locus_mis_pos,range(len(ref)), fileRadical, outputPath)
	# PLOT the score quality control violin plot
	#violin_matplot([score_normal, perm_score_normal, score_reverse, perm_score_reverse],ScoreTHR,fileRadical,outputPath)
	display('REPORTING...... DONE!\n')

	#Return summary data for batch summary
		#Get count mismatch total
	mismatch_list = list()
	for mismatch_pos in mismatch_index:
		mismatch_list.extend(mismatch_pos)

	mismatch_dict = Counter(mismatch_list)

	total_mismatch = sum(mismatch_dict.values())
	if total_mismatch != 0:
		perc_A = round(float(mismatch_dict['A'])*100/float(total_mismatch),2)
		perc_T = round(float(mismatch_dict['T'])*100/float(total_mismatch),2)
		perc_C = round(float(mismatch_dict['C'])*100/float(total_mismatch),2)
		perc_G = round(float(mismatch_dict['G'])*100/float(total_mismatch),2)
	else:
		perc_A = 0.0
		perc_T = 0.0
		perc_C = 0.0
		perc_G = 0.0

	HDRnoIndels_perc = round(HDRnoIndels_count*100/float(num_reads),2)
	NHEJ_perc = round(NHEJ_count*100/float(num_reads),2)
	NotEdited_perc = round(NotEdited_count*100/float(num_reads),2)
	HDRgap_perc = round(HDRgap_count*100/float(num_reads),2)

	stat_line = [num_reads,len(reads_uniq),DiscardedReads_count,(snp_index+1),snp_ref,snp_target,NotEdited_count, NotEdited_perc,HDRnoIndels_count,HDRnoIndels_perc,HDRgap_count,HDRgap_perc,NHEJ_count,NHEJ_perc,perc_A, perc_T, perc_C, perc_G]

	reporting_time = round(time.time() - reporting_start, 2)
	#return class_reads, insertion_index, deletion_index, match_index
	return [stat_line, load_time, align_time, analysis_time, reporting_time, reuse_lines, len(ref)]


def main():

	run_date = time.strftime("%m-%d-%Y_%H.%M.%S")

	try:
		input_directory = args['input']
		if input_directory[-1] != '/':
			input_directory += '/'
	except KeyError:
		input_directory = abs_path + 'input/'

	amplicon_directory = abs_path + 'amplicon_lib/'

	# New output directory management 'analyses'
	output_directory = abs_path + 'analyses/' + run_date + '/'

	if not os.path.exists(output_directory):
		os.makedirs(output_directory)


	settings_dict = get_settings(input_directory, args['settings'])

	# LOAD ARGUMENTS GENE LIST
	print settings_dict
	# settings_dict = {'mode':mode, 'input_name':input_name, 'reference':reference, 'template':template, 'gap_opening':gap_opening, 'gap_extension':gap_extension}

	if settings_dict['mode'] == 'uniq':
		fastq_list = [settings_dict['input_name']]
		output_directory = input_directory + 'analyses/' + settings_dict['reference'] + '_' + run_date + '/'
	elif settings_dict['mode'] == 'pair':
		fastq_list = [f for f in listdir(input_directory) if isfile(join(input_directory, f)) and settings_dict['input_name'] in f]
		output_directory = input_directory + 'analyses/' + settings_dict['reference'] + '_' + run_date + '/'
	elif settings_dict['mode'] == 'list':
		fastq_list = settings_dict['input_name'].split(';')
		output_directory = input_directory + 'analyses/' + settings_dict['reference'] + '_' + run_date + '/'
	else:   # Batch mode
		batch_directory = input_directory + settings_dict['input_name'] + '/'
		fastq_list = [f for f in listdir(batch_directory) if isfile(join(batch_directory, f)) and 'fastq' in f]
		output_directory = batch_directory + 'analyses/' + settings_dict['reference'] + '_' + run_date + '/'

	all_stats = list()
	log_data = list()

	for filename in fastq_list:
		if settings_dict['reference'] == '':
			reference = filename.split('-')[0]
		else:
			reference = settings_dict['reference']

		if settings_dict['template'] == '':
			template = reference
		else:
			template = settings_dict['template']

		inputReference = amplicon_directory + [ f for f in listdir(amplicon_directory) if isfile(join(amplicon_directory,f)) and '%s.ref' % reference in f][0]
		inputTemplate = amplicon_directory + [ f for f in listdir(amplicon_directory) if isfile(join(amplicon_directory,f)) and '%s.template' % template in f][0]
		if settings_dict['mode'] == 'batch':
			inputFile = batch_directory + filename
		else:
			inputFile = input_directory + filename
		fileRadical = filename.split('.fas')[0]
		outputPath = output_directory + fileRadical + '/'

		if not os.path.exists(outputPath):
			os.makedirs(outputPath)

		start_processing = time.time()

		print 'Analyzing file %s ...' % (filename)

		analysis_results = CRISPR_analysis(settings_dict['analysisMode'], inputFile, inputReference, inputTemplate, fileRadical, outputPath, settings_dict['gap_opening'], settings_dict['gap_extension'])

		if analysis_results == '':
			print 'EMPTY SAMPLE'
			all_stats.append([filename,0,0,0,'N/A','N/A','N/A',0,0,0,0,0,0,0,0,0,0,0,0])
			continue

		stat_line,load_time, align_time, analysis_time, reporting_time, reuse, ref_len = analysis_results

		total_time = round(time.time() - start_processing, 2)

		print '%s processed in: ' % filename, total_time, 's'
		print '\n'

		stat_line = [filename] + stat_line
		all_stats.append(stat_line)
		log_temp = [stat_line[0],stat_line[1],stat_line[2],reuse, ref_len,load_time, align_time, analysis_time, reporting_time,total_time]
		log_data.append(log_temp)

		# FileName, TotalReads, UniqueReads, DiscardedReads and %, SNPlocation, Reference, Template,
		# NHEJreads and %, NoIndelReads, HDRreads and %, mismatches %


	if settings_dict['mode'] == 'batch':
		f = open(output_directory + settings_dict['reference'] + '_' + 'BATCH_summary.txt', 'wb')
		header = ['FileName','totalReads','uniqueReads','discarded','SNPposition','Reference',\
		'Template','NotEdited','NotEdited%%','HDRclean','HDRclean%%','HDRgap','HDRgap%%','NHEJ','NHEJ%%','misA%%','misT%%','misC%%','misG%%']

		header = '\t'.join(header)
		f.write('%s\n' % header)
		for line in all_stats:
			line = [str(elem) for elem in line]
			to_write = '\t'.join(line)
			f.write('%s\n' % to_write)

		f.close()

	# Saving data to logs
	log_name = run_date + '_log.csv'
	g = open(abs_path + 'logs/%s' % log_name, 'wb')
	log_out = csv.writer(g)

	for line in log_data:
		log_out.writerow(line)
	g.close()

if __name__ == "__main__":
	main()
