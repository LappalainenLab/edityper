## CRISPR_TL v0.42
## Author: Alexandre Yahi
## Columbia University - New York Genome Center

This version has been developed with Python 2.7.

## 1. PYTHON PACKAGES INSTALLATION

Make sure the following Python packages are installed:
- matplotlib
- seaborn
- scipy
- numpy
- pandas

## 2. PREPARATING THE SETTINGS FILE

The settings must be specified in a text file at the root CRISPR_TL/, ex: my_settings.txt
Each variable/values pair is separated by ‘=‘ and there’s only one variable per line. Lines starting by ‘#’ are considered commented and won’t be read.

# 2.1. Analysis Modes

This tool currently support 3 analysis modes:
- SNP: detection of a SNP editing
- SNP+PAM: two edits, the SNP at the 5’, the PAM at the 3’
- PAM+SNP: same as previous mode but PAM at 5’, SNP at 3 ‘

# 2.2. Input name - processing mode

This tool infers the way it processes inputs based on the input name:
- Single file processing: input_name=sample.fastq
- Paired-end files processing: input_name=sample_R
	Will process all files starting by sample_R (e.g.,sample_R1.fastq, sample_R2.fastq)
- List of files: input_name=file_1.fastq;file_2.fastq;file_3.fastq
- Batch of file in a directory: input_name=directory_name
	Wille process all files under /CRISPR_TL/input/directory_name/

# 2.3. Amplicon library

Either specify the name of the amplicon directly, and then the tool will look for « AMPLICON.template.txt » and « AMPLICON.ref.txt » in the resources/ directory, or leave blank and the amplicon name will be inferred from the file name: « AMPLICON-X_Ri_blah.fastq ».

# 2.4. Gap penalties

We use affine gap penalty in our alignment algorithm. You can specify both gap opening and gap extension penalty: gap_opening/gap_extension. This value must be POSITIVE.

## 3. RUNNING AN ANALYSIS

Once the setting file prepared, simply enter in your command line window from the root directory CRISPR_TL/:

python CRISPR_TL.py settings=my_settings.txt

## 4. RESULTS

Each run creates a directory named with the timestamp of the run in CRISPR_TL/output/.
Each .fastq file has its own directory in the timestamp directory.
The current outputs are:
- fileName_events.txt
	Accounting of all the events with the reference’s index as relative reference.
Description of the columns:
	* Position: index of the position of the reference 
	* Reference: nucleotide of the reference
	* Coverage: count of the reads covering this section (exact nuc, mismatches, or deletions)
	* DelStart: count of deletion starting at this position
	* DelLengthAvg: average length of deletion starting at this position
	* DelNucCount: number of time this position is deleted
	* InsStart: count of insertion starting at this position
	* InsLengthAvg: average length of insertions starting at this position
	* A,T,C,G: count of mismatches and matches


- fileName_general.txt
	Summary of HDR, NHEJ, no edits and discarded reads

NO_EDIT: reads with no SNP edit, no insertions, no deletions. May contain mismatches.
Disc: discarded reads after quality filtering
HDR: reads with the wanted SNP edition
	* HDR without indels: HDR reads that may contain mismatches but no insertion nor deletion events.
	* HDR with deletions only: HDR reads that may contain mismatches, and at least one deletion (no inFisertions)
	* HDR with insertions only: HDR reads that may contain mismatches, and at least one insertion (no deletions)
	* HDR with indels: HDR reads with at least one insertion and one deletion, may contain mismatches
Insertion/Deletions/Mismatches events: total count of these events for each category
	- Avg: average number of given event per read, 
	- Std: standard deviation of this distribution of event.
NHEJ: reads without the wanted SNP edition, and at least one insertion or one deletion.

- fileName_locus.pdf
	A simple plot showing the various events relatively to the reference
- fileName_scores.pdf
	The violin plot showing the score distribution of the reads and their reverse complementary.

## 5. BATCH MODE SPECIAL

When in batch mode, a new tab-separated file is generated summarizing each file on one line.
The columns are:
FileName: name of the file
totalReads: total number of reads in this fastq file
uniqueReads: number of unique reads found in this file
discarded: number of reads discarded
SNPposition: position of the SNP relative to the reference
Reference: nucleotide at this position in the reference
Template: nucleotide at this position in the template (target SNP)
NotEdited: reads without indels and SNP edition, potential mismatches
HDRclean: number of HDR reads with no indels, may contain mismatches
NHEJ: number of reads without wanted SNP edition and at least one deletion, or one insertion
misA,T,C,G: percentage of mismatches found per nucleotide relative to the total number of mismatches 


