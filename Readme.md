# Genome Guesser

## Overview

**Genome Guesser** is a straightforward tool developed in Rust, designed to assist genomic researchers and bioinformaticians in identifying the reference genome used for variant calling and determining whether the variant positions in tabular variant files are based on 0-based or 1-based coordinates. This tool is currently in early development and not yet ready for production use.

## Key Features

- **Reference Genome Identification:** Quickly verifies a user-supplied reference genome could have been used for variant calling.
- **Coordinate System Detection:** Determines whether the variant positions are 0-based (used by bedvar) or 1-based (used by VCFs and MAFs), facilitating accurate data interpretation and analysis.
- **Compatibility:** Accepts tabular variant files with a header, including columns for chromosome (`chrom`), position (`pos`), and reference bases (`ref`). The tool intelligently identifies these columns regardless of order, provided they are named sensibly.

## Limitations

Identifying the exact reference genome for variant calling from a simple tabular is not always possible, especially with a limited set of variants and similar potential reference genomes. 
Consider a scenario where you aim to discern whether the variants in your file correspond to hg19 or hg38. Providing only one of these genomes to Genome Guesser and a mere handful of variants might lead the tool to affirm a match with the supplied genome. This outcome, while technically accurate since no discrepancies with the provided reference were detected might not reflect the complete picture. For a more robust verification, it's crucial to re-run Genome Guesser with each potential reference genome. In the given example, this approach would reveal that the five variants align perfectly with both hg19 and hg38, rendering the results inconclusive. 

Hence, the key takeaway is that to truly determine the basis of your variant file with respect to one genome or another, it's advisable to test against all conceivable references in Genome Guesser and evaluate the degree of concordance across them.


## Getting Started

### Installation

Currently, Genome Guesser is in early development. Installation instructions will be provided upon release.

### Usage

To use Genome Guesser, run the following command:

```shell
genomeguesser <variants> <genome>

Arguments:
  <variants>  TSV file with variants. Expects columns: 'Chrom', 'Pos', and 'Ref', but will also look for synonymous column names
  <genome>    FASTA file with genome to test

Options:
  -h, --help     Print help
  -V, --version  Print version
```