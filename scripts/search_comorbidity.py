#!/usr/bin/python

import os
import csv


def read_matrix(matrix_file):
    matrix = []
    mfile = open(matrix_file, 'r')
    mreader = csv.reader(mfile, delimiter = '\t')
    next(mreader, None)
    for row in mreader:
        matrix.append(row)
    return(matrix)

def find_traits(query_traits, matrix):
    comorbid = []
    traits = {}
    results = []
    query_traits = query_traits.split(',')
    for trait in query_traits:
        for row in matrix:
            if trait == row[0] or trait == row[1]:
                comorbid.append(row)
    
    for row in comorbid:
        xtrait = row[0]
        ytrait = row[1]
        if xtrait not in traits:
            traits[xtrait] = get_genes(xtrait)
        if ytrait not in traits:
            traits[ytrait] = get_genes(ytrait)
        common_genes = list(set(traits[xtrait].keys()).intersection(traits[ytrait].keys()))
        common_list = common_genes[0]
        if len(common_genes) > 1:
            for i in range(1, len(common_genes)):
                common_list += ' ,' + common_genes[i]
        results.append([row[0], row[1], row[2], row[3], common_list, row[4], row[5]])
    return(results)


def get_genes(trait):
    trait_dict = {}
    trait_dir = os.path.join('../data/downloaded/test/', trait)
    eqtl_file = open(os.path.join(trait_dir, 'sig_SNP-gene_eqtls.txt'), 'r')
    ereader = csv.reader(eqtl_file, delimiter = '\t')
    next(ereader, None)
    for row in ereader:
        gene = row[3]
        snp = row[0]
        tissue = row[7]
        if gene not in trait_dict:
            trait_dict[gene] = {'tissues':[], 'snps':[]}
        if snp not in trait_dict[gene]['snps']:
            trait_dict[gene]['snps'].append(snp)
        if tissue not in trait_dict[gene]['tissues']:
            trait_dict[gene]['tissues'].append(tissue)

    return(trait_dict)



if __name__=='__main__':
    matrix_file = '../data/downloaded/test/gene_matrix.txt'
    query_traits = 'comprehensive_strength_index__muscle_measurement'
    matrix = read_matrix(matrix_file)
    results = find_traits(query_traits, matrix)
    rfile = open(os.path.join('../analysis', query_traits) + '_query.txt', 'w')
    rwriter = csv.writer(rfile, delimiter = '\t')
    rwriter.writerow(['Trait_X', 'Trait_Y', '#Trait_X_Genes', '#Common_Genes', \
                          'Common_Genes', 'Per_Overlap', 'P_value'])
    rwriter.writerows(results)
