#!/usr/bin/env python

import sys
import csv
import numpy as np
import scipy.stats

class ClassCollection(object):
    gene_names = []
    ngenes = 0
    tissues = []
    celltypes = []
    classes = []

    def __init__(self):
        self.gene_names = []
        self.classes = []
        self.celltypes = []
        self.tissues = []

    def add_tissue(self, tissue):
        if tissue not in self.tissues:
            self.tissues.append(tissue)

    def add_celltype(self, celltype):
        if celltype not in self.celltypes:
            self.celltypes.append(celltype)

    def readCentroids(self, filename):
        with open(filename, "r") as f:
            c = csv.reader(f, delimiter='\t')
            line = c.next()
            self.gene_names = line[2:]
            self.ngenes = len(self.gene_names)

            self.gene_names = [g.upper() for g in self.gene_names]

            for line in c:
                cc = CellClass(self, line)
                self.classes.append(cc)
        sys.stderr.write("{} cell classes\n{} tissues\n{} cell types\n{} genes\n".format(len(self.classes), len(self.tissues), len(self.celltypes), self.ngenes))

    def vectorsFromCellranger(self, cellranger):

        with open(cellranger, "r") as f:
            c = csv.reader(f)
            hdr = c.next()
            nclust = (len(hdr) - 2) / 3

        sys.stderr.write("{} clusters in input file.\n\n".format(nclust))
        for i in range(nclust):
            sys.stderr.write("** Classifying cluster {}\n".format(i+1))
            g = []                  # gene names
            v = []                  # from cellranger
            col = 2 + i * 3

            with open(cellranger, "r") as f:
                c = csv.reader(f)
                c.next()
                for line in c:
                    g.append(line[1])
                    v.append(float(line[col]))
                        
            result = []
            for cl in self.classes:
                (rho, p) = cl.eval(g, v)
                result.append((rho, cl.tissue, cl.celltype))

            result.sort(key=lambda r: r[0], reverse=True)
            sys.stderr.write("  => Tissue: {}, CellType: {} ({})\n\n".format(result[0][1], result[0][2], result[0][0]))
            outfile = "cl{}.classification.csv".format(i+1)
            with open(outfile, "w") as out:
                for r in result:
                    out.write("{}\t{}\t{}\n".format(*r))
                    

class CellClass(object):
    tissue = ""
    celltype = ""
    vector = None

    def __init__(self, parent, row):
        self.vector = {}
        self.tissue = row[0]
        self.celltype = row[1]
        parent.add_tissue(self.tissue)
        parent.add_celltype(self.celltype)
        gene_names = parent.gene_names
        rp = 2
        for gene in gene_names:
            self.vector[gene] = float(row[rp])
            rp += 1

    def eval(self, genes, values):
        ng = 0
        nz = 0
        z = np.zeros(len(genes))
        idx = 0
        for g in genes:
            if g in self.vector:
                z[idx] = self.vector[g]
                if z[idx] == 0:
                    nz += 1
                ng += 1
            idx += 1
        return scipy.stats.spearmanr(z, values)

if __name__ == "__main__":
    cc = ClassCollection()
    cc.readCentroids(sys.argv[1])
    cc.vectorsFromCellranger(sys.argv[2])
