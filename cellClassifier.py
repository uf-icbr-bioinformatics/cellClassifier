#!/usr/bin/env python

import sys, csv
from scipy.stats import fisher_exact

def findIntersection(g1, g2):
    return [g for g in g1 if g in g2]

def readNclusters(filename):
    """Read the header of a cellranger differential_expression file and return the number of clusters it contains."""
    with open(filename, "r") as f:
        hdr = f.readline().split(",")
        return (len(hdr) - 2) / 3

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
            self.candidates.append([pval, dbrow[0], dbrow[1], len(common), ",".join(common)])

    def writeResults(self):
        if self.outfile:
            self.stream = open(self.outfile, "w")
        try:
            self.candidates.sort(key=lambda a: a[0])
            for c in self.candidates:
                self.stream.write("\t".join([c[1], c[2], str(c[0]), str(c[3]), c[4]]) + "\n")
        finally:
            if self.outfile:
                self.stream.close()
        
class Manager(object):
    mode = "normal"             # or "cellranger"
    dbfile = None
    filename = None
    column = 0
    outfile = None
    classifiers = []
    totgenes = None
    pval = 0.05
    Xpval = 0.05
    Xfc = 2

    def __init__(self):
        self.classifiers = []

    def parseArgs(self, args):
        if "-h" in args or "--help" in args:
            return self.usage()
        prev = ""
        for a in args:
            if prev == "-p":
                self.pval = float(a)
                prev = ""
            elif prev == "-n":
                self.totgenes = int(a)
                prev = ""
            elif prev == "-o":
                self.outfile = a
                prev = ""
            elif prev == "-c":
                self.column = int(a) - 1
                prev = ""
            elif prev == "-Xp":
                self.Xpval = float(a)
                prev = ""
            elif prev == "-Xfc":
                self.Xfc = float(a)
                prev = ""
            elif a in ["-n", "-p", "-o", "-c", "-Xp", "-Xfc"]:
                prev = a
            elif a == "-X":
                self.mode = "cellranger"
            elif self.dbfile is None:
                self.dbfile = a
            elif self.filename is None:
                self.filename = a
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
            if not self.totgenes:
                if l1[0] == "# Genes:":
                    self.totgenes = int(l1[1])
                else:
                    raise "Error: malformed database."
            for line in c:
                C.classify(line, self.totgenes)
        # When done, print results
        C.writeResults()

    def runX(self):
        self.readCellrangerGenes()
        with open(self.dbfile, "r") as f:
            cr = csv.reader(f, delimiter='\t')
            l1 = cr.next()
            if not self.totgenes:
                if l1[0] == "# Genes:":
                    ngenes = int(l1[1])
                else:
                    raise "Error: malformed database."
            for line in cr:
                for c in self.classifiers:
                    c.classify(line, self.totgenes)
        # When done, print results
        for c in self.classifiers:
            c.writeResults()
        
    def usage(self):
        sys.stdout.write("""cellClassifier.py - Classify cells based on high-expressed genes.

Usage: cellClassifier [options] dbfile genesfile

This program reads a list of genes from `genesfile' and compares them against
the cell signatures contained in `dbfile' using Fisher's exact test. Cell types
that match are printed to standard output (or to a specified file).

""")
        
if __name__ == "__main__":
    M = Manager()
    if M.parseArgs(sys.argv[1:]):
        M.run()


# def findCellType(dbfile, genes):
#     nmine = len(genes)
#     with open(dbfile, "r") as f:
#         c = csv.reader(f, delimiter='\t')
#         l1 = c.next()
#         if l1[0] == "# Genes:":
#             totgenes = int(l1[1])
#         else:
#             raise "Error: malformed database."

#         candidates = []
#         for line in c:
#             ng = int(line[2])
#             cg = line[3].split(",")
#             common = findIntersection(cg, genes)
#             A = len(common)
#             B = ng - A
#             C = nmine - A
#             D = totgenes - A - B - C
#             table = [[A, B], [C, D]]
#             odds, pval = fisher_exact(table, alternative="greater")
#             if pval < 0.05:
#                 candidates.append([pval, line[0], line[1], len(common), ",".join(common)])

#         candidates.sort(key=lambda a: a[0])
#         for c in candidates:
#             sys.stdout.write("\t".join([c[1], c[2], str(c[0]), str(c[3]), c[4]]) + "\n")
    
