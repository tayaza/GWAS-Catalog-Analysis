#!/usr/bin/python
import os
import csv
import xlrd
import pandas
import random

def get_egenes(trait_files):
    egenes = []
    trait_genes = {}
    for tfile in trait_files:
        trait = tfile[len(eqtl_dir):len(tfile)-23]
        trait_genes[trait] = []
        t_file = open(tfile, 'r')
        treader = csv.reader(t_file, delimiter = '\t')
        next(treader, None)
        for row in treader:
            if row[3] not in egenes:
                egenes.append(row[3])
            if row[3] not in trait_genes[trait]:
                trait_genes[trait].append(row[3])
    print("Number of eGenes without duplicates: " + str(len(egenes)))
    return(egenes, trait_genes)

def get_mgenes(mendeliome):
    mgenes = []
    mrows = []
    mfile = xlrd.open_workbook(mendeliome)
    shnames = mfile.sheet_names()
    msheet = mfile.sheet_by_index(0)
    for rx in range(msheet.nrows):
        cx = msheet.row_values(rx)
        if cx[0] not in mgenes:
            mgenes.append(cx[0])
        mrows.append([cx[0], cx[7]])
    print(len(mgenes))
    return(mgenes, mrows)
                     
def compare_genes(egenes, trait_genes, mgenes, mrows):
    ofile = open('../analysis/mendiliome.txt', 'w')
    common = []
    for mg in mgenes:
        if mg in egenes:
            common.append(mg)
    fishers = 2 * len(common) / float(len(egenes) + len(mgenes))
    print(len(common))
    print(fishers)
    ofile.write('#GWAS genes = ' + str(len(egenes)) + '\n')
    ofile.write('#Mendeliome genes = ' + str(len(mgenes)) + '\n')
    ofile.write('#common genes = ' + str(len(common)) + '\n')
    ofile.write('Fishers exact test = ' + str(fishers) + '\n\n')

    cwriter = csv.writer(ofile, delimiter = '\t')
    cwriter.writerow(['Gene', 'Panel', 'GWAS trait'])
    genes = {}
    for gene in common:
        genes[gene] = []
        for trait in trait_genes:
            if gene in trait_genes[trait]:
                genes[gene].append(trait)
        trait_list = genes[gene][0]
        if len(genes[gene]) > 1:
            for i in range(1, len(genes[gene])):
                trait_list += ', ' + genes[gene][i]
        genes[gene] = trait_list
        for row in mrows:
            if row[0] == gene:
                cwriter.writerow([str(gene), str(row[1]), genes[gene]])

def control(trait_genes):
    genepool = []
    controls = {}
    for trait in trait_genes:
        for gene in trait_genes[trait]:
            genepool.append(gene)
    print("Number of genes with redundancy: " + str(len(genepool)))
    for trait in trait_genes:
        gene_len = len(trait_genes[trait])
        controls[trait] = []
        if gene_len > 0:
            for i in range(0, gene_len):
                rand_gene_index = random.randint(0, len(genepool))
                if rand_gene_index >= len(genepool):
                    rand_gene_index = rand_gene_index - 1
                print(trait, str(len(controls[trait])+1)+ '/'+  str(gene_len), \
                          genepool[rand_gene_index], rand_gene_index)
                controls[trait].append(genepool[rand_gene_index])
                genepool.remove(genepool[rand_gene_index])
                #print('\t' + str(len(genepool)) + ' remaining')
    control_dir = '../data/controls'
    matrix = []
    if not os.path.isdir(control_dir):
        os.mkdir(control_dir)
    for trait in controls:
        gene_count = {}
        genes_total = len(controls[trait])
        compare_traits = controls.keys()
        if genes_total > 3:
            for trait_gene in controls[trait]:
                for compare in compare_traits:

                    if trait_gene in controls[compare]:
                        if compare not in gene_count.keys():
                            gene_count[compare] = 1
                        else:
                            gene_count[compare] += 1
                        
            row = []
            row.append(trait)
            for compare in gene_count:
                total_compare_genes = len(controls[compare])
                pvalue = 1 - (2*gene_count[compare] \
                    /float((genes_total + total_compare_genes)))
                ratio = round(gene_count[compare]/float(genes_total), 7)
                matrix.append([trait, compare, genes_total, gene_count[compare], ratio, pvalue])

        tfile = open(os.path.join(control_dir, trait+'_control.txt'), 'w')
        for gene in controls[trait]:
            tfile.write(gene+'\n')
        tfile.close()

    with open (os.path.join(control_dir, 'controls_matrix.txt'), 'w') as cluster_file:
        writer = csv.writer(cluster_file, delimiter = '\t')
        writer.writerow(['trait_x', 'trait_y', '#total_genes', '#common_snps', \
                             'ratio', 'pvalue'])
        writer.writerows(matrix)
    


if __name__ == '__main__':
    eqtl_dir = '../data/downloaded/test/'
    mendeliome = '../data/13059_2015_693_MOESM4_ESM.xls'
    trait_files = []
    if os.path.isdir(eqtl_dir):
        _,trait_dir,_ = next(os.walk(eqtl_dir), (None, [], None))
    trait_files = [os.path.join(eqtl_dir,tfile +'/sig_SNP-gene_eqtls.txt') \
                       for tfile in trait_dir]
    egenes, trait_genes = get_egenes(trait_files)
    mgenes, mrows = get_mgenes(mendeliome)
    compare_genes(egenes, trait_genes, mgenes, mrows)
    control(trait_genes)
