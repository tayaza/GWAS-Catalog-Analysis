#usr/bin/python
import os 
import sqlite3
import csv


"""
Date: 22052017
Purpose: To compare the gene and eQTL SNP correlation matrices 
"""

def read_files(trait_file):
    records = []
    if not os.path.isfile(trait_file):
        print("Not a file!")
    else:
        t_file = open(trait_file, 'r')
        t_records = csv.reader(t_file, delimiter = '\t')
        if trait_file != './correlations/snps/correlation_matrix/traits.txt':
            next(t_records, None)
        for row in t_records:
            records.append(row)
        t_file.close()
    return records


def compare(snp_file, gene_file, trait_file):
    combined = []
    snp_matrix = read_files(snp_file)
    gene_matrix = read_files(gene_file)
    traits = read_files(trait_file)
    traits = [trait[0] for trait in traits]
    for trait in traits:
        snp_data = []
        gene_data = []
        
        for snp_rec in snp_matrix:
            snp_xtrait = snp_rec[0]
            snp_ytrait = snp_rec[1]
            if snp_xtrait == trait and snp_ytrait in traits and \
                    snp_rec not in snp_data:
                snp_data.append(snp_rec)
                for gene_rec in gene_matrix:
                    gene_xtrait = gene_rec[0]
                    gene_ytrait = gene_rec[1]
                    if snp_xtrait == gene_xtrait and snp_ytrait == gene_ytrait: 
                        print(snp_rec)
                        print(gene_rec)
                        print('\n')
                        combined.append(snp_rec + gene_rec[2:])

    print(len(combined))
                
    print(len(gene_matrix), len(snp_matrix))

    #TODO: fix output to give 625 traits
    matrix_file = open('./correlations/snp_gene_matrix1.txt', 'w')
    twriter = csv.writer(matrix_file, delimiter = '\t')
    theader = ('xtrait', 'ytrait', '#snps', 'common_snps', 'snp_ratio', \
                   '#genes', 'common_genes', 'gene_ratio')
    twriter.writerow(theader)
    twriter.writerows(combined)
    matrix_file.close()

if __name__ == "__main__":
    snp_file = './correlations/snps/correlation_matrix/all_matrix.txt'
    gene_file = '../analysis/gene_matrix.txt'
    trait_file = './correlations/snps/correlation_matrix/traits.txt'
    compare(snp_file, gene_file, trait_file)


