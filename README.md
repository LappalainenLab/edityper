# EdiTyper

A python program for aligning CRISPR-edited RNA-seq data. Written by Alexandre Yahi.

## Installation

To install, run the following command

```bash
python setup.py install
```

Please note, to install system-wide, you must have sudo access. If you are not installing system-wide, please install within your `$PYTHONPATH`.

## Dependencies

EdiTyper is compatible with Python 2.7+ and 3.3+; it also depends on the following Python modules:
 - [Cython](http://cython.org/)
 - [NumPy](http://www.numpy.org/)
 - [SciPy](https://www.scipy.org/)
 - [matplotlib](http://matplotlib.org/)
 - [Biopython](http://biopython.org/)
 - [regex](https://pypi.python.org/pypi/regex)

Each of these modules is available on [PyPi](https://pypi.python.org/) and installable using [pip](https://pip.pypa.io/en/stable/)

## Basic Usage

To get a basic help message, simply run

```bash
EdiTyper
```

The minimal arguments needed to run are `-r | --reference-sequence`, `-t | --template-sequence`, and one of the following: `-i | --input file`, `-l | --sample-list`, or `-d | --fastq-directory`

```bash
EdiTyper -r reference.fasta -t template.fasta -i sample.fastq
```

## Classification scheme

We classify the unique reads based on their alignment with the reference. HDR are reads with the right SNP edited to the target SNP as identified in the template sequence, and no deletion and no insertions but tolerate mismatches. If there are insertions or deletions then the read is classified as MIX, which is a HDR with indels. If the SNP of interest is not edited and there are no indels, then the read is unchanged (NO_EDIT). Otherwise, the read is a NHEJ.

## Arguments

### Alignment Arguments

| Parameter | Definition | Required? | Default |
| --------- | ---------- | --------- | ------- |
| `-p | --pvalue-threshold` | Set the p-value threshold for the alignment | No | 1 * 10<sup>-3</sup> |
| `-g | --gap-opening` | Set the gap opening penalty | No | 8 |
| `-e | --gap-extension` | Set the gap extension penalty | No | 1 |
| `--parallel` | Run EdiTyper in parallel, supports optional argument to limit number of cores, otherwise uses all available cores | No | None |

### Reference Arguments

| Parameter | Definition | Required? |
| --------- | ---------- | --------- |
| `-r | --reference-sequence` | Reference FASTA for alignment | **Yes** |
| `-t | --template-sequence` | Template sequence for CRISPR editing | **Yes** |
| `-b | --reference-bed` | A BED file with the genomic location of the reference sequence. Used to get chromosome name and sequence start point. Only the first entry will be read; any line starting with 'browser', 'track' or '`#`' will be ignored | No |

### Input Arguments

> **NOTE**: These arguments are mutually exclusive
>
> **NOTE**: At least one of these arguments is required

| Parameter | Definition |
| --------- | ---------- |
| `-i | --input-file` | Provide a FASTQ file for aligning, this option can be specified multiple times for multiple FASTQs |
| `-l | --sample-list` | Provide a list of FASTQ files for aligning; there should be one FASTQ file per line |
| `-d | --fastq-directory` | Provide a directory where FASTQ files are located; files must end in .fastq or .fq (.gz is allowed afterwards) |

### Output Arguments

| Parameter | Definition | Required? | Default |
| --------- | ---------- | --------- | ------- |
| `-d | --output-directory` | Choose where to put all output files | No | 'output' |
| `-j | --project` | Give this project a name, will be used as the basename for output files | No | 'edityper' |
| `--bam` | Convert SAM files to BAM files, optionally use CSI indecies instead of BAI indeces with `--bam csi`; requires [SAMtools](https://github.com/samtools/samtools) | No | None |

### Read Group Arguments

> **NOTE**: Information provided will be applied to *all* read groups in *all* FASTQ files

| Parameter | Definition | Required? | Default |
| --------- | ---------- | --------- | ------- |
| `-rc | --read-center` | Name of the sequencing center | No | None |
| `-rl | --read-library` | Sequencing library | No | None |
| `-rp | --read-platform` | Platform used for sequencing | No | None |
| `-rs | --read-sample` | Sample being sequenced | No | None |

### Suppression Arguments

| Parameter | Definition |
| --------- | ---------- |
| `--suppress-sam` | Suppress SAM output |
| `--suppress-events` | Suppress events table output |
| `--suppress-classification` | Suppress read classification output |
| `--suppress-tables` | Suppress both events and read classification |
| `--suppress-plots` | Suppress locus and quality plots |

## Output Files

For each output table, all lines starting with `#` are header lines. All lines starting with `##` are extra information. See details below for specifics about each output.

| Output file | Extension |
| ----------- | --------- |
| Alignments in SAM/BAM format | `.sam | .bam` |
| Table of events by base | `.events.txt` |
| Classification of reads in tabular format | `.classification.txt` |
| Locus and alignment quality plots | `.pdf` |
| Summary of read classifications per input FASTQ file | `.summary.txt` |
| Read assignments table | `.assignments.txt` |

### SAM/BAM Output

SAM alignments are standard SAM files. They have been presorted in coordinate order with read groups attached. If BAM output, BAM indices will also be output in either BAI or CSI format. For more details, see SAM/BAM format specification from [HTSlib](http://www.htslib.org/).

### Events Table

The `.events` table shows a locus-by-locus overview of indels and mismatches in each FASTQ file. One table is generated per FASTQ file.

<!-- reference state, coverage, number of deletions, average deletion length, number of insertions, average insertion length, and base counts at each position in reference sequence. -->

| Column | Meaning |
| ------ | ------- |
| `POS` | Position in reference sequence |
| `REF` | Nucleotide in reference sequence at this position |
| `COV` | Coverage in FASTQ at position |
| `DEL` | Number of deletions starting at this position |
| `AVG_DEL` | Average length of deletions starting at this position |
| `DCOUNT` | Number of times this position is deleted |
| `INS` | Number of insertions starting at this position |
| `AVG_INS` | Average length of insertions starting at this position |
| `A` | Count of mismatched A's at this position |
| `T` | Count of mismatched T's at this position |
| `C` | Count of mismatched C's at this position |
| `G` | Count of mismatched G's at this position |

### Read Classifications

The `.classifications.txt` table shows a breakdown of indels and mismatches per read category. One table is generated per FASTQ file. The read categories are HDR, MIX, NHEJ, NO_EDIT, and DISCARD.

<!-- counts, indels, and mismatches for HDR, NHEJ, no editing, and discarded reads. In addition, this file shows SNP state and position, read counts, and alignment scoring information. -->

| Column | Meaning |
| ------ | ------- |
| `TAG` | Which classification category, one of HDR, MIX, NHEJ, NO_EDIT, and DISCARD |
| `COUNT` | How many reads fall in this classification category? |
| `PERC_COUNT` | What percentage of reads fall in this classification category. Excludes discarded reads from calculation |
| `INS_EVENTS` | Total number of insertion events; reported for HDR, MIX, and NHEJ only |
| `AVG_INS` | Average number of insertion events per read; reported for HDR, MIX, and NHEJ only |
| `STD_DEV_INS` | Standard deviation of the distribution of insertions; reported for HDR, MIX, and NHEJ only |
| `DEL_EVENTS` | Total number of deletion events; reported for HDR, MIX, and NHEJ only |
| `AVG_DEL` | Average number of deletion events per read; reported for HDR, MIX, and NHEJ only |
| `STD_DEV_DEL` | Standard deviation of the distribution of deletions; reported for HDR, MIX, and NHEJ only |
| `MISMATCH_EVENTS` | Total number of mismatch events; reported for HDR, MIX, and NHEJ only |
| `AVG_MIS` | Average number of mismatch events per read; reported for HDR, MIX, and NHEJ only |
| `STD_DEV_MIS` | Standard deviation of the distribution of mismatches; reported for HDR, MIX, and NHEJ only |
| `NO_INDELS` | Number of reads with no indels; reported for HDR and MIX only |
| `PERC_NO_INDELS` | Percentage of reads with no indels; reported for HDR and MIX only |
| `INS_ONLY` | Number of reads with only insertions; reported for HDR and MIX only |
| `PERC_INS_ONLY` | Percentage of reads with only insertions; reported for HDR and MIX only |
| `DEL_ONLY` | Number of reads with only deletions; reported for HDR and MIX only |
| `PERC_DEL_ONLY` | Percentage of reads with only deletions; reported for HDR and MIX only |
| `INDELS` | Number of reads with both insertions and deletions; reported for HDR and MIX only |
| `PERC_INDELS` | Percentage of reads with both insertions and deletions; reported for HDR and MIX only |

### Summary Table

The `.summary.txt` table shows total reads, unique reads, discarded reads, SNP information, no editing, HDR, NHEJ, and mismatch percentages by base. One table is generated for *all* FASTQ files.

| Column | Meaning |
| ------ | ------- |
| `FASTQ` | Name of FASTQ file |
| `TOTAL_READS` | Total number of reads in this FASTQ file |
| `TOTAL_NON_DISC` | Total number of reads in this FASTQ file, excluding discarded reads |
| `UNIQ_READS` | Total number of unique non-discarded reads |
| `DISCARDED` | Number of discarded reads |
| `SNP_POS` | SNP position |
| `REF_STATE` | Reference state |
| `TEMP_SNP` | Alternate state |
| `NO_EDIT` | Number of non-edited reads |
| `PERC_NO_EDIT` | Percentage of non-edited reads, excluding discarded reads |
| `HDR_CLEAN` | Number of clean HDR (not MIX) reads |
| `PERC_HDR_CLEAN` | Percentage of clean HDR (not MIX) reads, excluding discarded reads |
| `HDR_GAP` | Number of gapped HDR (MIX) reads |
| `PERC_HDR_GAP` | Percentage of gapped HDR (MIX) reads, excluding discarded reads |
| `NHEJ` | Number of NHEJ reads |
| `NHEJ_GAP` | Percentage of NHEJ reads, excluding discarded reads |
| `PERC_MIS_A` | Percentage of mismatches with an A compared to total number of mismatches |
| `PERC_MIS_T` | Percentage of mismatches with an T compared to total number of mismatches |
| `PERC_MIS_C` | Percentage of mismatches with an C compared to total number of mismatches |
| `PERC_MIS_G` | Percentage of mismatches with an G compared to total number of mismatches |

### Locus and Quality Plots

The locus plot shows events at each base along the reference and the number of supporting reads for each event. One PDF is generated per FASTQ file. The first page is scaled to total number of reads in the FASTQ file, the second is scaled to maximum number of supporting reads accross all events. Coverage at each base is also shown.

The quality plot shows the distribution of alignment quality scores.  One PDF is generated for *all* FASTQ files. The quality-score threshold for discarding reads is shown as a black bar.

### Assignments Table

The `.assignments.txt` table shows how each read was classified as well as the number of insertions, deletions, and mismatches for each read. One table is generated per FASTQ. This table is *only* generated when verbosity is set to 'debug' (`-v debug | --verbosty debug`) and `--suppress-tables` is **not** passed.

| Column | Meaning |
| ------ | ------- |
| `ReadID` | Read ID from the FASTQ file |
| `Label` | Classification assigned to this read |
| `NumDel` | Number of deletions in this read |
| `NumIns` | Number of insertions in this read |
| `NumMis` | Number of mismatches in this read |

<!-- ---

<img src='.nygc.jpeg' alt='New York Genome Center' height='100' width='100'> -->
