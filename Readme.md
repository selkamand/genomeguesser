# Genome Guesser

Genome Guesser is in early development and not yet ready for use.

A simple rust tool that checks whether tabular variant files are 0, or 1 based, and whether they align to user-specified reference genomes.


## Problem Statement

Ever stumbled accross a genomic variant file out in the wild, and wanted to confirm 

1) what reference genome was used to call these variants?
2) whether the variant positions are 0 or 1-based coordinates. VCFs and MAFs use 1based coordinates while bedvar uses 0-based

## Solution

Genome guesser can help confirm  variants were called against a particular reference genome - and whether coordinates are 0 or 1-based. This is accomplished by simplying pulling all SNVs with A/C/T/G reference bases and checking whether these reference bases match user-supplied reference genome assuming either 0 or 1-based positions (fasta).

Genome Guesser expects input file to have a header and include columns chrom, pos and ref. Column order doesn't matter, but they must be named something sensible so the tool stands a chance of guessing which is the chromosome, position and reference columns.

guess_genome -i variant_file.tsv -f <path_to_genomefasta>
