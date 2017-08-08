#!usr/bin/python
import os
import csv


def parse_clusters(clusterfile):
    clusters = {}
    if os.path.isfile(clusterfile):
        clfile = open(clusterfile, 'r')
        clreader = csv.reader(clfile, delimiter = '\t')
        i = 0
        for line in clreader:
            if line:
                if line[0].startswith("Cluster"):
                    i += 1
                    clusters[i] = []
                else:
                    if  line[0] != '':
                        clusters[i].append(line[0])
    else:
        print('Not a file')
    return(clusters)

def parse_matrix(matrixfile):
    matrix = []
    if os.path.isfile(matrixfile):
        mfile = open(matrixfile, 'r')
        mreader = csv.reader(mfile, delimiter = '\t')
        next(mreader, None)
        for line in mreader:
            matrix.append(line)

    return(matrix)


def calc_indices(aList):
    #print(aList)
    list_mean = sum(aList) / float(len(aList))
    list_min = min(aList)
    list_max = max(aList)
    return(list_mean, list_min, list_max)


def process(clusters, matrix):
    for cluster in clusters:
        genes = []
        ratio = []
        trait_num = 0
        for trait in clusters[cluster]:
            trait_num += 1
            for row in matrix:
                for cotrait in clusters[cluster]:
                    if row[0] == trait and row[1] == cotrait and \
                            row[0] != row[1]: #Exclude traits matching themselves
                        genes.append(int(row[3]))
                        ratio.append(float(row[4]))
        genes_mean, genes_min, genes_max = calc_indices(genes)
        ratio_mean, ratio_min, ratio_max = calc_indices(ratio)
        print(cluster, trait_num, clusters[cluster])
        #print(cluster, genes_mean, genes_min, genes_max)
        #print(cluster, ratio_mean, ratio_min, ratio_max)
        

if __name__ == '__main__':
    clusterfile = './analysis/trait_clusters.txt'
    matrixfile = './analysis/gene_matrix.txt'
    clusters = parse_clusters(clusterfile)
    matrix = parse_matrix(matrixfile)
    process(clusters, matrix)
