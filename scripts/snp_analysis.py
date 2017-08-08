#!usr/bin/python

import os
import sqlite3
import csv

"""
Retrieves spatial eQTL SNPs from the 'spatial_eqtls' database
and produces an output containing: 1) the genes which with the SNP regions
interact with 2) the distance of the interactions and 3) whether these
interactions are cis (i.e. < 1mb) or trans (i.e. > 1mb or intra-chromosomal)
"""


def extract_snps(snp_dir):
    trait_snps = {}
    if os.path.isdir(snp_dir):
        _,_,trait_files = next((os.walk(snp_dir)), (None,None,[]))
        for trait_file in trait_files:
            trait = trait_file[:len(trait_file)-4]
            #trait_snps[trait] = []
            tfile = open(os.path.join(snp_dir, trait_file), 'r')
            snps = tfile.readlines()
            snps = [snp.strip() for snp in snps]
            trait_snps[trait] = snps
    
    return trait_snps



def process(eqltDB, trait_snps):
    conn = sqlite3.connect(eqtlDB)
    conn.text_factory = str
    cur = conn.cursor()

    snp_dict = {}
    cur.execute("SELECT * FROM eqtls;")
    cur_data = cur.fetchall()

    for row in cur_data:
        snp = row[0]
        gene = row[3]
        distance  = row[13]
        interaction = row[12]
        if interaction == 'True':
            interaction = 'Cis'
        else:
            interaction = 'Trans'
        if snp not in snp_dict:
            snp_dict[snp] = {}
        snp_dict[snp][gene] = {'dist':distance, 'inter':interaction}
    
    snp_file = open('../analysis/eqtl_snps.txt', 'w')
    swriter = csv.writer(snp_file, delimiter = '\t')
    sheader = ['SNP', '#Genes', 'Gene', 'Interaction', 'SNP-Gene_distance', \
                   '#Traits', 'Traits']
    swriter.writerow(sheader)
    
    for snp in snp_dict:
        trait_list = []
        for trait in trait_snps:
            if snp in trait_snps[trait]:
                trait_list.append(trait)
        if trait_list:
            traits = trait_list[0]
        if len(trait_list) > 1:
            for i in range (1, len(trait_list)):
                traits += ' ,' + trait_list[i]
        for gene in snp_dict[snp]:
            swriter.writerow([snp, len(snp_dict[snp]), gene, \
                      snp_dict[snp][gene]['inter'], \
                      snp_dict[snp][gene]['dist'], \
                      len(trait_list), traits])
    snp_file.close()

if __name__ == '__main__':
    eqtlDB = './spatial_eqtls.db'
    snp_dir = '../data/downloaded/test'
    trait_snps = extract_snps(snp_dir)
    process(eqtlDB, trait_snps)
