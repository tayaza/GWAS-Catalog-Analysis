#usr/bin/python
import os 
import sqlite3
import csv



def getSNPs(traits_dir):
    traits = {}
    print('Getting trait names')
    if os.path.isdir(traits_dir):
        for trait_file in os.listdir(traits_dir):
            if os.path.isfile(traits_dir + '/' + trait_file) and \
            trait_file.endswith('.txt') and \
            trait_file != 'all_catalog_snps.txt' and \
            trait_file != 'error_log.txt':
                trait = (trait_file[:len(trait_file)-4])
                tfile = open(traits_dir + '/' + trait_file, 'r')
                snps = tfile.readlines()
                snps = [x.strip() for x in snps]
                traits[trait] = {}
                for snp in snps:
                    traits[trait][snp] = []
    else:
        print('There is no such directory as ', traits_dir)

    return traits

def queryDB(eqtlDB, traits):
    conn = sqlite3.connect(eqtlDB)
    conn.text_factory = str
    cur = conn.cursor()
    genes = {}
    gene_details = []
    gene_summary = []
    theader = ('SNP', 'SNP_Chromosome', 'SNP_Locus',
               'Gene_Name', 'Gene_Chromosome', 'Gene_Start', 'Gene_End',
               'Tissue', 'p-value', 'q-value', 'Cell_Lines',
               'GTEx_cis_p_Threshold', 'cis_SNP-gene_interaction',
               'SNP-gene_distance', 'Expression_Level_in_eQTL_Tissue')
    if not os.path.isdir('./sig_eqtls'):
        print('Creating directory for significant eQTLs')
        os.mkdir('./sig_eqtls')
    print('Writing summary for... ')
    for trait in traits:
        print('\t'+ trait)
        tfile = open('./sig_eqtls/' + trait + '.txt', 'w')
        twriter = csv.writer(tfile, delimiter = '\t')
        twriter.writerow(theader)
        for snp in traits[trait]:
            cur.execute("SELECT * FROM eqtls WHERE snp=?;", (snp,))
            qdata = cur.fetchall()
            if qdata:
                for row in qdata:
                    traits[trait][snp].append(row)
                    twriter.writerow(row)
                    gene = row[3]
                    interaction = ''
                    if row[12] == 'True':
                        interaction = 'cis'
                    else:
                        interaction = 'trans'
                    if not gene in genes.keys():
                        genes[gene] = {'snps': [], 'traits': [], 'tissues': [],
                                       'interactions': []}
                    genes[gene]['snps'].append(snp)
                    genes[gene]['traits'].append(trait)
                    genes[gene]['tissues'].append(row[7])
                    genes[gene]['interactions'].append(interaction)
        tfile.close()
    conn.close()

    print('Processing eGene matrix...')
    for gene in genes:
        snplist = genes[gene]['snps']
        s_snps = set(snplist)
        snps = snplist[0]
        if len(snplist) > 1:
            for i in range (1, len(snplist)):
                snps = snps + ', ' + snplist[i]

        traitlist = list(genes[gene]['traits'])
        s_traitlist = set(traitlist)
        traits = traitlist[0]
        if len(traitlist) > 1:
            for i in range (1, len(traitlist)):
                traits = traits + ', ' + traitlist[i]

        tissuelist = list(genes[gene]['tissues'])
        s_tissuelist = set(tissuelist)
        tissues = tissuelist[0]
        if len(tissuelist) > 1:
            for i in range (1, len(tissuelist)):
                tissues = tissues + ', ' + tissuelist[i]

        interlist = list(genes[gene]['interactions'])
        s_interlist = set(interlist)
        inters = interlist[0]
        inter = ''
        if len(set(genes[gene]['interactions'])) == 1:
            inter = list(set(genes[gene]['interactions']))[0]
        else:
            inter = 'both'
        if len(interlist) > 1:
            for i in range (1, len(interlist)):
                inters = inters + ', ' + interlist[i]

        gene_summary.append((gene, len(s_snps), len(s_traitlist), \
                                 len(s_tissuelist), inter))
        gene_details.append((gene, snps, traits, tissues, inters))

    print('Writing eGene files...')
    sfile = open('./gene_summary.txt', 'w')
    swriter = csv.writer(sfile, delimiter = '\t')
    swriter.writerow(('Gene', '#SNPs', '#Traits', '#Tissues', 'Interactions'))
    swriter.writerows(gene_summary)
    dfile = open('./gene_details.txt', 'w')
    dwriter = csv.writer(dfile, delimiter = '\t')
    dwriter.writerow(('Gene', 'SNPs', 'Traits', 'Tissues', 'Interactions'))
    dwriter.writerows(gene_details)
        


if __name__ == "__main__":
    eqtlDB = './spatial_eqtls.db'
    traits_dir = '../data/downloaded/test'
    traits = getSNPs(traits_dir)
    queryDB(eqtlDB, traits)

