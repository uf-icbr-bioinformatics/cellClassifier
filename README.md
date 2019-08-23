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

## Database format
The database file supplied as the first argument should be a
tab-delimited file with the following format:

- The first line should contain two fields, the first one being "# Genes:" and the second
one containing the total number of genes listed in the database.

- All successive lines should contain four fields:
  1. Tissue
  2. Cell type
  3. Number of genes in cell type signature
  4. Names of genes in signature, separated by a single comma.

For example, the following are the first two lines of a database file:

```
# Genes:        11541
Amniotic fluid  Amniotic fluid stem cell  4  CD44,ENG,NT5E,THY1
```

Four prebuilt databases are provided in the `databases` directory. The source for this databases is the [CellMarker](http://bio-bigdata.hrbmu.edu.cn/CellMarker/index.jsp) website. If you use any of these databases in your work, please cite the CellMarker paper: 

[CellMarker: a manually curated resource of cell markers in human and mouse](https://academic.oup.com/nar/advance-article/doi/10.1093/nar/gky900/5115823). Nucleic Acids Research. 2018. 

## Cellranger support
If `-X` is specified, the program assumes that the input is a 
`differential_expression.csv` file produced by cellranger. The file is comma-delimited
with two columns for gene identifier and gene name, followed by three columns for
each cluster containing average expression, fold change, and P-value, respectively.

In this mode, cellClassifier will first extract highly-expressed genes from each cluster
i.e. those genes with a fold change higher than the value specified with -Xfc (2 by default)
and a P-value smaller than the one specified with -Xp (0.05 by default). It will then create one
classifier for each set of genes.

The output of each classifier will be written to a file called `cc.clustN.csv`, where N is the
cluster number. For example, the classifications for the set of genes in the third cluster
will be written to cc.clust3.csv. The `cc` prefix can be changed with the -o option.
