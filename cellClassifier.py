#!/usr/bin/env python

import sys, csv
import os.path
import math
from scipy.stats import fisher_exact

# Utils

def findIntersection(g1, g2):
    return [g for g in g1 if g in g2]

def readNclusters(filename):
    """Read the header of a cellranger differential_expression file and return the number of clusters it contains."""
    with open(filename, "r") as f:
        hdr = f.readline().split(",")
        return (len(hdr) - 2) / 3

def safeInt(v):
    try:
        return int(v)
    except ValueError:
        sys.stderr.write("Error: the value `{}' should be an integer number.")
        sys.exit(1)
    
def safeFloat(v):
    try:
        return float(v)
    except ValueError:
        sys.stderr.write("Error: the value `{}' should be a number.")
        sys.exit(1)

def safeFilename(v):
    if os.path.isfile(v):
        return v
    else:
        sys.stderr.write("Error: file `{}' does not exist or is not readable.")
        sys.exit(1)

def fileAndColumn(v):
    p = v.rfind(":")
    if p > 0:
        try:
            x = int(v[p+1:])
            return v[:p], x - 1
        except:
            return v, 0
    return v, 0

def helpWanted(args):
    h = False
    for a in args:
        if h:
            return a
        elif a in ["-h", "--help"]:
            h = True
        else:
            h = False
    return h

class Classifier(object):
    """Class that stores a list of genes, to test them against the markers for different cell types."""
    filename = ""
    column = 0
    genes = []
    ngenes = 0
    pval = 0.05
    candidates = []
    outfile = None
    stream = sys.stdout

    def __init__(self):
        self.genes = []
        self.candidates = []
        
    def readGenesFromFile(self):
        with open(self.filename, "r") as f:
            c = csv.reader(f, delimiter='\t')
            for line in c:
                if self.column >= len(line):
                    continue
                g = line[self.column]
                if g not in self.genes:
                    self.genes.append(g)
        self.ngenes = len(self.genes)
        return self.genes

    def addGene(self, gene):
        if gene not in self.genes:
            self.genes.append(gene)
            self.ngenes += 1

    def classify(self, dbrow, totgenes):
        ng = int(dbrow[2])
        cg = dbrow[3].split(",")
        common = findIntersection(cg, self.genes)
        A = len(common)
        B = ng - A
        C = self.ngenes - A
        D = totgenes - A - B - C
        table = [[A, B], [C, D]]
        odds, pval = fisher_exact(table, alternative="greater")
        if pval < self.pval:
            score = - math.log(pval, 10.0) * math.log(A, 10.0)
            self.candidates.append([pval, score, dbrow[0], dbrow[1], A, ",".join(common)])

    def writeResults(self, sortby=0):
        rev = (sortby == 1)
        if self.outfile:
            self.stream = open(self.outfile, "w")
        try:
            self.candidates.sort(key=lambda a: a[sortby], reverse=rev)
            for c in self.candidates:
                self.stream.write("\t".join([c[2], c[3], str(c[0]), str(c[1]), str(c[4]), c[5]]) + "\n")
        finally:
            if self.outfile:
                self.stream.close()
        
class Manager(object):
    mode = "normal"             # or "cellranger"
    dbfile = None
    filename = None
    column = None
    outfile = None
    classifiers = []
    totgenes = None
    pval = 0.05
    sortby = 0                # 0=pval, 1=score
    Xpval = 0.05
    Xfc = 2

    def __init__(self):
        self.classifiers = []

    def parseArgs(self, args):
        w = helpWanted(args)
        if w:
            return self.usage(w)
        prev = ""
        for a in args:
            if prev == "-p":
                self.pval = safeFloat(a)
                prev = ""
            elif prev == "-n":
                self.totgenes = safeInt(a)
                prev = ""
            elif prev == "-o":
                self.outfile = a
                prev = ""
            elif prev == "-c":
                self.column = safeInt(a) - 1
                prev = ""
            elif prev == "-Xp":
                self.Xpval = safeFloat(a)
                prev = ""
            elif prev == "-Xfc":
                self.Xfc = safeFloat(a)
                prev = ""
            elif a in ["-n", "-p", "-o", "-c", "-Xp", "-Xfc"]:
                prev = a
            elif a == "-s":
                self.sortby = 1
            elif a == "-X":
                self.mode = "cellranger"
            elif self.dbfile is None:
                self.dbfile = safeFilename(a)
            elif self.filename is None:
                f, col = fileAndColumn(a)
                self.filename = safeFilename(f)
                if not self.column:
                    self.column = col
        if self.dbfile and self.filename:
            return True
        else:
            self.usage()
            return False

    def initialize(self):
        if self.mode == "normal":
            C = Classifier()
            self.classifiers.append(C)
            C.pval = self.pval
            C.filename = self.filename
            C.column = self.column
            C.outfile = self.outfile
        elif self.mode == "cellranger":
            if not self.outfile:
                self.outfile = "cc"
            nclusters = readNclusters(self.filename)
            for i in range(nclusters):
                c = Classifier()
                c.pval = self.pval
                c.outfile = "{}.clust{}.csv".format(self.outfile, i+1)
                self.classifiers.append(c)

    def readCellrangerGenes(self):
        nclust = len(self.classifiers)
        with open(self.filename, "r") as f:
            f.readline()        # skip header
            c = csv.reader(f, delimiter=",")
            for line in c:
                for i in range(nclust):
                    fc = float(line[i*3 + 3])
                    pv = float(line[i*3 + 4])
                    if pv < self.Xpval and abs(fc) >= self.Xfc:
                        self.classifiers[i].addGene(line[1])
        for i in range(nclust):
            sys.stderr.write("Cluster {} ({}): {} genes.\n".format(i+1, self.classifiers[i].outfile, self.classifiers[i].ngenes))
                
    def run(self):
        self.initialize()
        if self.mode == "cellranger":
            return self.runX()
        C = self.classifiers[0]
        C.readGenesFromFile()
        with open(self.dbfile, "r") as f:
            c = csv.reader(f, delimiter='\t')
            l1 = c.next()
            if self.totgenes is None:
                if l1[0] == "# Genes:":
                    self.totgenes = int(l1[1])
                else:
                    raise "Error: malformed database."
            for line in c:
                C.classify(line, self.totgenes)
        # When done, print results
        C.writeResults(self.sortby)

    def runX(self):
        self.readCellrangerGenes()
        with open(self.dbfile, "r") as f:
            cr = csv.reader(f, delimiter='\t')
            l1 = cr.next()
            if self.totgenes is None:
                if l1[0] == "# Genes:":
                    self.totgenes = int(l1[1])
                else:
                    raise "Error: malformed database."
            for line in cr:
                for c in self.classifiers:
                    c.classify(line, self.totgenes)
        # When done, print results
        for c in self.classifiers:
            c.writeResults(self.sortby)
        
    def usage(self, what=None):
        sys.stdout.write("""cellClassifier.py - Classify cells based on highly expressed genes.
""")
        if what == "cellranger":
            sys.stdout.write("""
* CellRanger mode. In this mode, the program assumes that the input is a 
differential_expression.csv file produced by cellranger. The file is comma-delimited
with two columns for gene identifier and gene name, followed by three columns for
each cluster.

In this mode, cellClassifier will first extract highly-expressed genes from each cluster
i.e. those genes with a fold change higher than the value specified with -Xfc ({} by default)
and a P-value smaller than the one specified with -Xp ({} by default). It will then create one
classifier for each set of genes.

The output of each classifier will be written to a file called cc.clustN.csv, where N is the
cluster number. For example, the classifications for the set of genes in the third cluster
will be written to cc.clust3.csv. The `cc' prefix can be changed with the -o option.

""".format(self.Xfc, self.Xpval))
        elif what == "db":
            sys.stdout.write("""
* Database format. The database file supplied as the first argument should be a
tab-delimited file with the following format:

- The first line should contain two fields, the first one being "# Genes:" and the second
one containing the total number of genes listed in the database.

- All successive lines should contain four fields:
  1. Tissue
  2. Cell type
  3. Number of genes in cell type signature
  4. Names of genes in signature, separated by a single comma.

For example, the following are the first two lines of a database file:

# Genes:        11541
Amniotic fluid  Amniotic fluid stem cell  4  CD44,ENG,NT5E,THY1

""")
        elif what == "output":
            sys.stdout.write("""
For each entry in the signatures database, the program performs a Fisher exact test
comparing the genes in the signature with the provided ones. If the test is successful,
the entry is considered a candidate match. In addition to the P-value, the program
computes a match score using the following formula:

  S = -log10(P-value) * log10(I)

where I is the number of genes in the intersection between the provided list and the
genes in the signature. Higher values of the score indicate a better match. This score
favors the entries for which the intersection is a large number of genes; therefore
sorting based on this score (using the -s option) may provide a somewhat different 
ordering compared to the default one based on the P-value.

The output consists of a tab-delimited file with five columns: 

* Tissue
* Cell type
* P-value
* Score
* List of genes from input list matching the gene signature. 

Results are sorted in order of increasing P-value, or decreasing match score is `-s' is supplied.

""")
        else:
            sys.stdout.write("""
Usage: cellClassifier.py [options] dbfile genesfile

This program reads a list of genes from `genesfile' and compares them against
the cell signatures contained in `dbfile', looking for significant matches. It
will return a list of the cell type signatures that best match the provided list
with a P-value and a match score. Use `-h db' for a description of the format of 
the database file, and `-h output' for a description of the program's output.

A signature matches if the P-value returned by the enrichment test is less than 0.05 
(this limit can be changed with the -p option). The test uses the number of genes in the database 
as the total number of genes, but this can be changed with -n. Cell types that 
match are printed to standard output (or to the file specified with the -o option). 

`genesfile' is assumed to contain one gene identifier per line, or to be a 
tab-delimited file with identifiers in the first column. A different column can be 
specified with the -c option or using the syntax filename:column.

If -X is specified, the program switches to cellranger mode, suitable for parsing
differential expression files produced by cellranger. Use -h cellranger for details.

Options:

  -o O   | Write output to file O (default: standard output).
  -p P   | Set the P-value threshold for signature matching to P (default: {}).
  -c C   | Column containing gene names in input file (default: 1).
  -n N   | Total number of genes considered (default: number of genes in database).
  -s     | Sort results by score rather than P-value.
  -X     | Enable cellranger mode. See -h cellranger.
  -Xp P  | Set P-value threshold for cellranger file to P (default: {}).
  -Xfc F | Set fold change threshold for cellranger file to F (default: {}).

""".format(self.pval, self.Xpval, self.Xfc))
        sys.stdout.write("""(c) 2019, A. Riva, ICBR Bioinformatics Core, University of Florida.

""")                 
        
if __name__ == "__main__":
    M = Manager()
    if M.parseArgs(sys.argv[1:]):
        M.run()

