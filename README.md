# cellClassifier
Automatic classification of clusters from single-cell RNA-Seq experiments

```
Usage: cellClassifier.py [options] dbfile genesfile
```

## Introduction

This program reads a list of genes from `genesfile` and compares them against
the cell signatures contained in `dbfile` using Fisher's exact test. See the
Database section for a description of the format of the database file.

A signature matches if the P-value returned by the test is less than 0.05 (can 
be changed with the -p option). The test uses the number of genes in the database 
as the total number of genes, but this can be changed with -n. Cell types that 
match are printed to standard output (or to the file specified with the -o option). 

The output consists of four columns: tissue, cell type, P-value, list of genes
from input list matching the gene signature. 

`genesfile` is assumed to contain one gene identifier per line, or to be a 
tab-delimited file with identifiers in the first column. A different column can be 
specified with the -c option or using the syntax filename:column.

If -X is specified, the program switches to cellranger mode, suitable for parsing
differential expression files produced by cellranger. See the Cellranger Support
section for details.

## Options

Option   | Description
---------|-----------------
  -o O   | Write output to file O (default: standard output).
  -p P   | Set the P-value threshold for signature matching to P (default: 0.05).
  -c C   | Column containing gene names in input file (default: 1).
  -n N   | Total number of genes considered (default: number of genes in database).
  -X     | Enable cellranger mode. See the Cellranger Support section.
  -Xp P  | Set P-value threshold for cellranger file to P (default: 0.05).
  -Xfc F | Set fold change threshold for cellranger file to F (default: 2).

