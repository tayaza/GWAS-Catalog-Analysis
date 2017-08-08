#!usr/bin/python
import csv

def find_genes():
    snp_dhs = {}
    with open(dhsfile) as dfile:
        reader = csv.reader(dfile, delimiter = '\t')
        next(reader, None)
        for line in reader:
            snp = line[0]
            snp_chr = line[1]
            snp_pos = line[2]
            dhs_id = line[3]
            cluster_id = line[7]
            celltype = line[10]
            dhs_tissue = line[11]
            disease = line[13]
            to_dict = []
            genes = []
            if snp not in snp_dhs.keys():
                to_dict.append((snp_chr, snp_pos, dhs_id, cluster_id, \
                                    celltype, dhs_tissue, disease))
            else:
                to_dict = snp_dhs[snp]
                to_dict.append((snp_chr, snp_pos, dhs_id, cluster_id, \
                                    celltype, dhs_tissue, disease))
            to_dict.append(genes)
            snp_dhs[snp] = to_dict

    with open(t2d_eqtl_file) as tfile:
        reader = csv.reader(tfile, delimiter = '\t')
        next(reader, None)
        for line in reader:
            snp = line[0]
            gene = line[4]
            if snp in snp_dhs.keys():
                to_dict = snp_dhs[snp][1]
                to_dict.append(gene)
                snp_dhs[snp][1] = to_dict
        
    with open(obesity_eqtl_file) as ofile:
        reader = csv.reader(ofile, delimiter = '\t')
        next(reader, None)
        for line in reader:
            snp = line[0]
            gene = line[4]
            if snp in snp_dhs.keys():
                to_dict = snp_dhs[snp][1]
                to_dict.append(gene)
                snp_dhs[snp][1] = to_dict

    return snp_dhs

def sort_tissues(snp_genes):
    tissues = {}
    novel_genes = []
    with open(dhs_tissues_file) as dfile:
        reader = csv.reader(dfile, delimiter = '\t')
        next(reader, None)
        for line in reader:
            cell_line = line[0]
            tissue = line[1]
            cell_lines = []
            if tissue not in tissues.keys():
                cell_lines.append(cell_line)
            else:
                cell_lines = tissues[tissue]
                cell_lines.append(cell_line)
            tissues[tissue] = cell_lines
    dfile.close()

    with open(novel_genes_file) as nfile:
        reader = csv.reader(nfile, delimiter = '\t')
        next(reader, None)
        for line in reader:
            gene = line[0]
            novel_genes.append(gene)
    nfile.close()

    with open('../analysis/all_snp_tissues.txt', 'wb') as wfile:
        writer = csv.writer(wfile, delimiter = '\t')
        writer.writerow(('Gene', 'SNP', 'SNP_Celltype', 'SNP_tissue', 'Disease'))
        for snp in snp_genes.keys():
            genes = snp_genes[snp][1]
            celltype = snp_genes[snp][0][4]
            disease = snp_genes[snp][0][6]
            for tissue in tissues.keys():
                if celltype in tissues[tissue]:
                    for gene in genes:
                        #if gene in novel_genes:
                        to_file = gene, snp, celltype, tissue, disease
                        writer.writerow(to_file)
    wfile.close()

if __name__== '__main__':
    dhsfile = '../analysis/t2d-obesity_dhs.txt' # A merge of T2D and obesity snp_dhs.txt
    t2d_eqtl_file = '../diabetes/results/dhs_results/match.txt'
    obesity_eqtl_file = '../obesity/results/codes3d_results/dhs_results/match.txt'
    dhs_tissues_file = '../analysis/dhs_cells_tissues.txt'	# A concordance of tissue origin of cell lines
    novel_genes_file = '../analysis/novel/codes3d_novel_associations.txt'

    snp_genes = find_genes()
    sort_tissues(snp_genes)
