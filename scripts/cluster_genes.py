#! /usr/bin/env python

import sys
import csv
import os
import argparse

"""
Find common eGenes among traits and calculate the ratio of occurrence.
Takes in as input significant (FDR < 0.05)spatial eQTL SNP-gene interactions
for each GWAS Catalog trait.
"""

def resolve_output_fp(input_fp, output_fp):
    print('Resolving IO parameters..')
    resolved_output_fp = ''
    if output_fp == 'NA':
        resolved_output_fp = input_fp
    else:
        resolved_output_fp = output_fp
    if not resolved_output_fp.endswith('/'):
        resolved_output_fp = resolved_output_fp + '/'
    if not os.path.isdir(resolved_output_fp):
        os.mkdir(resolved_output_fp)
    print(resolved_output_fp)
    return resolved_output_fp
                        


def find_common_genes(input_fp):
    """
    Find common Genes from all trait associations downloaded from GWAS Catalogue.
    """
    trait_genes = {}
    all_genes = []
    common_genes = []
    snp_count = {}
    traits = {}
    matrix = []
    print('Extracting genes from eQTL interactions for...')
    _,_,t_files = next(os.walk(input_fp), (None, None, []))
    for trait_file in t_files:
        trait = trait_file[:len(trait_file)-4]
        print('\t' + trait)
        tfile = open(os.path.join(input_fp, trait_file), 'r')
        eqtls= csv.reader(tfile, delimiter = '\t') 
        next(tfile, None)
        for line in eqtls:
            genes = []
            if trait in trait_genes.keys():
                genes = trait_genes[trait]
            genes.append(line[3])
            trait_genes[trait] = genes
            all_genes.append(line[3])
        tfile.close()
    
    for trait in trait_genes:
        trait_genes[trait] = list(set(trait_genes[trait]))
    all_genes = list(set(all_genes))
    print(len(all_genes))

    done_genes = []
    """
    for snp in all_snps:
        occur = all_snps.count(snp)
        if occur > 1 and snp not in done_snps:
            done_snps.append(snp)
            for record in trait_snps:
                if snp == record[1] and record not in common_snps:
                    common_snps.append(record)
                    snp_count[snp] = occur
                    to_dict = []
                    if record[0] not in traits.keys():
                        to_dict.append(snp)
                        traits[record[0]] = to_dict
                    else:
                        to_dict = traits[record[0]]
                        to_dict.append(snp)
                        traits[record[0]] = to_dict
    """
    for trait in trait_genes.keys():
        gene_count = {}
        genes_total = len(trait_genes[trait])
        compare_traits = trait_genes.keys()
        if genes_total > 3:
            for trait_gene in trait_genes[trait]:
                for compare in compare_traits:
                    if trait_gene in trait_genes[compare]:
                        if compare not in gene_count.keys():
                            gene_count[compare] = 1
                        else:
                            gene_count[compare] += 1
                    #else:
                    #    gene_count[compare] = 0
            row = []
            row.append(trait)
            for t in gene_count:
                ratio = round(gene_count[t]/float(genes_total), 7)
                matrix.append([trait, t, genes_total, gene_count[t], ratio])

    """
    with open (output_fp + '/' + 'common_snps_count.txt', 'wb') as cluster_file:
        writer = csv.writer(cluster_file, delimiter = '\t')
        writer.writerow(['snp', 'count'])
        for snp in snp_count:
            writer.writerow([snp,snp_count[snp]])
    """

    with open ('gene_matrix.txt', 'w') as cluster_file:
        writer = csv.writer(cluster_file, delimiter = '\t')
        writer.writerow(['trait_x', 'trait_y', '#total_genes', '#common_snps', \
                             'ratio'])
        writer.writerows(matrix)
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = "")
    parser.add_argument("-i", "--input", required = True, \
                            help = "File of GWAS Catalog traits and SNPs") 

    parser.add_argument("-o", "--output", default = 'NA', \
                            help = "Output directory") 

    args = parser.parse_args()
    eqtl_dir = args.input
    output_fp = resolve_output_fp(args.input, args.output)
    find_common_genes(args.input)
