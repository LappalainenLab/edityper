# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

#AUTHOR: ALEXANDRE YAHI
#AFFILATION: LAPPALAINEN LAB - New York Genome Center - Columbia University
#DATE: 2015
#VERSION: 0.2

'''Genetic ploting of CRISPR analytics tool 0.2'''

import os
from os import listdir
from os.path import isfile, join
import re
import sys
import csv
import gzip
# Data structure packages
from collections import defaultdict
from operator import itemgetter
# TIME MANAGEMENT packages
from datetime import date, datetime
from dateutil.relativedelta import relativedelta
import time
import random
# MATH packages
import numpy as np
from scipy.stats import norm
# PLOTTING
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
# Force matplotlib to not use any Xwindows backend.

import pandas as pd

# RELATIVE PATH
script_name =  os.path.basename(__file__)
abs_path = os.path.realpath(__file__).split(script_name)[0] 

def violin_quality(input_list,scoreTHR,fileRadical, outputPath):
	sns.set(style="whitegrid")

	labels = ['-','-*','RC','RC*']

	ar = np.array(input_list)
	ar = np.transpose(ar)

	df = pd.DataFrame(ar, columns=labels)

	f, ax = plt.subplots(figsize=(4, 6))
	sns.violinplot(df, palette="Set3",  bw=.5, cut=0.5, linewidth=1)
	ax.set(ylim=(-.7, max(df.max(axis=0))+100))
	sns.despine(left=True, bottom=True)

	plt.axhline(y=scoreTHR, linewidth=0.5, color = 'r')

	# Finalize the figure
	plt.savefig(outputPath + '%s_scores.pdf' % fileRadical)


def locus_map(cumul_del, ins_start, mis_pos, length, fileRadical, outputPath):

	ar = list()

	for i, elem in enumerate(cumul_del):
		ar.append([i, 'Del', float(elem)])
	for i, elem in enumerate(ins_start):
		ar.append([i, 'Ins start', float(elem)])
	for i, elem in enumerate(mis_pos):
		ar.append([i, 'Mismatch', float(elem)])


	ar = np.array([cumul_del,ins_start,mis_pos])
	#ar = np.array(ar)
	ar = np.transpose(ar)


	labels = ['Position','Event','counts']

	df = pd.DataFrame(ar, columns=['Deletion','Insertion','Mismatch'])

	f, ax = plt.subplots(figsize=(25, 5))
	df.plot(lw=1)
	#ax = sns.tsplot(data=df, time="Position", value="counts", condition="Event",err_style=None, n_boot=0)
	#ax = sns.tsplot(data=ar,err_style=None)

	plt.savefig(outputPath + '%s_locus.pdf' % fileRadical)
