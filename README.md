# crispr.py

A python program for aligning CRISPR-edited RNA-seq data. Written by Alexandre Yahi and Paul Hoffman.

## Basic Usage

To get a basic help message, simply run

```bash
./crispr.py
```

There are two subroutines as part of crispr.py: `CONFIG` and `ALIGN`. To view detailed help information about either subroutine, simply tell crispr.py to run the desired subroutine and pass the `-h` flag. For example:

```bash
./crispr.py CONFIG -h
```

for viewing help information about the `CONFIG` subroutine.

## Dependencies

crispr.py is compatible with Python 2.7 only. It also depends on the following Python modules:
 - [matplotlib](http://matplotlib.org/)
 - [NumPy](http://www.numpy.org/)
 - [SciPy](https://www.scipy.org/)

Each of these modules is available on [PyPi](https://pypi.python.org/) and installable using [pip](https://pip.pypa.io/en/stable/)

## `CONFIG`

The `CONFIG` subroutine creates a configuration file for the `ALIGN` subroutine. Using the `CONFIG` subroutine, one can provide alignment parameters, set reference and template arguments, choose an input FASTQ file or a sample list, set the output directory and project name, and provide extra read group information for SAM output. The configuration file is stored in [INI](https://en.wikipedia.org/wiki/INI_file) format. Please see the help information for more information

### Alignment Arguments

| Parameter | Definition | Required? | Default |
| --------- | ---------- | --------- | ------- |
| `-p | --pvalue-threshold` | Set the p-value threshold for the alignment | No | 1 * 10<sup>-3</sup> |
| `-g | --gap-opening` | Set the gap opening penalty | No | 8 |
| `-e | --gap-extension` | Set the gap extension penalty | No | 1 |

### Reference Arguments

| Parameter | Definition | Required? |
| --------- | ---------- | --------- |
| `-r | --reference-sequence` | Reference FASTA for alignment | **Yes** |
| `-t | --template-sequence` | Template sequence for CRISPR editing | **Yes** |

### Alignment Arguments

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
| `-j | --project` | Give this project a name, will be used as the basename for output files | No | 'crispr'

### Read Group Arguments

> **NOTE**: Information provided will be applied to *all* read groups in *all* FASTQ files

| Parameter | Definition | Required? | Default |
| --------- | ---------- | --------- | ------- |
| `-rc | --read-center` | Name of the sequencing center | No | None |
| `-rl | --read-library` | Sequencing library | No | None |
| `-rp | --read-platform` | Platform used for sequencing | No | None |
| `-rs | --read-sample` | Sample being sequenced | No | None |

## `ALIGN`

The `ALIGN` subroutine aligns CRISPR-edited RNA-seq data using the parameters stored in a configuration file. This file is easily generated using the `CONFIG` subroutine. It will generate several output files:

| Output file | Extension |
| ----------- | --------- |
| Alignments in SAM format | `.sam` |
| Table of insertion, deletion, and mismatch events | `.events` |
| Classification of reads in tabular format | `.classification` |
| Locus and alignment quality plots | `.pdf` |
| Summary of read classifications per input FASTQ file | `.summary` |

Any combination of these outputs can be suppressed, please see the help message for more information.

[<img src='.nygc.jpeg' alt='New York Genome Center' height='100' width='100'>](http://www.nygenome.org/)

