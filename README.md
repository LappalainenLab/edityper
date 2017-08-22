# CRISPRonto

A python program for aligning CRISPR-edited RNA-seq data. Written by Alexandre Yahi.

## Installation

To install, run the following command

```bash
python setup.py install
```

Please note, to install system-wide, you must have sudo acces. If you are not installing system-wide, please install within your `$PYTHONPATH`.

## Dependencies

CRISPRonto is compatible with Python 2.7+ and 3.3+; tt also depends on the following Python modules:
 - [Cython](http://cython.org/)
 - [matplotlib](http://matplotlib.org/)
 - [NumPy](http://www.numpy.org/)
 - [SciPy](https://www.scipy.org/)
 - [Biopython](http://biopython.org/)

Each of these modules is available on [PyPi](https://pypi.python.org/) and installable using [pip](https://pip.pypa.io/en/stable/)

## Basic Usage

To get a basic help message, simply run

```bash
CRISPRonto
```

## Arguments

### Alignment Arguments

| Parameter | Definition | Required? | Default |
| --------- | ---------- | --------- | ------- |
| `-p | --pvalue-threshold` | Set the p-value threshold for the alignment | No | 1 * 10<sup>-3</sup> |
| `-g | --gap-opening` | Set the gap opening penalty | No | 8 |
| `-e | --gap-extension` | Set the gap extension penalty | No | 1 |
| `-n | --num-cores` | Set the number of cores for multiprocessing | No | Number of cores available or number of FASTQ files provided, whichever is lower |

### Reference Arguments

| Parameter | Definition | Required? |
| --------- | ---------- | --------- |
| `-r | --reference-sequence` | Reference FASTA for alignment | **Yes** |
| `-t | --template-sequence` | Template sequence for CRISPR editing | **Yes** |

### Input Arguments

> **NOTE**: These arguments are mutually exclusive
>
> **NOTE**: At least one of these arguments is required

| Parameter | Definition |
| --------- | ---------- |
| `-i | --input-file` | Provide a single FASTQ file for aligning |
| `-l | --sample-list` | Provide a list of FASTQ files for aligning; there should be one FASTQ file per line |

### Output Arguments

| Parameter | Definition | Required? | Default |
| --------- | ---------- | --------- | ------- |
| `-d | --output-directory` | Choose where to put all output files | No | 'output' |
| `-j | --project` | Give this project a name, will be used as the basename for output files | No | 'crispronto'

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

For each output table, all lines starting with `#` are header lines. All lines starting with `##` are extra information.

| Output file | Extension |
| ----------- | --------- |
| Alignments in SAM format | `.sam` |
| Table of events by base | `.events` |
| Classification of reads in tabular format | `.classification` |
| Locus and alignment quality plots | `.pdf` |
| Summary of read classifications per input FASTQ file | `.summary` |

SAM alignments are standard SAM files. They have been presorted in coordinate order with read groups attached.

The `.events` table shows reference state, coverage, number of deletions, average deletion length, number of insertions, average insertion length, and base counts at each position in reference sequence.

The `.classifications` table shows counts, indels, and mismatches for HDR, NHEJ, no editing, and discarded reads. In addition, this file shows SNP state and position, read counts, and alignment scoring information.

The `.summary` table shows total reads, unique reads, discarded reads, SNP information, no editing, HDR, NHEJ, and mismatch percentages by base per FASTQ file.

The locus plot shows events at each base along the reference and the number of supporting reads for each event.

The quality plot shows the distribution of alignment quality scores for each FASTQ file.

<img src='.nygc.jpeg' alt='New York Genome Center' height='100' width='100'>
