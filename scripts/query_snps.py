#!usr/bin/python
import csv
import argparse
import os



def set_filepath(input_fp, output_fp):
    """ Set output directory filepath."""
    if output_fp == 'default':
        filepath = input_fp
        if ('/' in filepath):
            while not (filepath.endswith('/')):
                filepath = filepath[:len(filepath)-1]
            filepath = filepath[:len(filepath)-1]
        else:
            filepath = ''
    else:
        filepath = output_fp
    if not os.path.isdir(filepath):
        os.mkdir(filepath)

    return filepath
                                                
def parse_snps(snp_file):
    snps = []
    sfile = open (snp_file, 'rb')
    reader = csv.reader(sfile, delimiter = '\t')
    next(reader, None)
    for line in reader:
        snps.append(line[0])

    return snps

def query_data(snps, data_dir):
    dir_list = os.listdir(data_dir)
    eqtls = []
    failed_batch = []
    for result_dir in dir_list:
        summary = data_dir + '/' + result_dir + '/summary.txt'
        if os.path.isfile(summary):
            print '\t Checking batch ' + result_dir
            sfile = open(summary, 'rb')
            reader = csv.reader(sfile, delimiter = '\t')
            next(reader, None)
            for snp in snps:
                for line in reader:
                    if snp == line[0]:
                        eqtls.append(line)
        else:
            failed_batch.append(result_dir)
            print result_dir
    outfile = open(output_fp + '/codes3d_results.txt', 'wb')
    writer = csv.writer(outfile, delimiter = '\t')
    header = 'SNP', 'SNP_Chromosome', 'SNP-Locus', 'Gene_Name', \
                'Gene_Chromosome', 'Gene_Start', 'Gene_End', 'Tissue', 'p-value', \
                'q-value', 'Cell_Lines', 'GTEx_cis_p_Threshold', \
                'cis_SNP_Gene_Interaction', 'SNP-gene_Distance', \
                'Expression_Level_in_eQTL_Tissue', 'Max_Expressed_Tissue', \
                'Maximum_Expression_Level', 'Min_Expressed_Tissue', \
                'Min_Expressed_Level'
    writer.writerow(header)
    writer.writerows(eqtls)
    outfile.close()
    failed_batches = open(output_fp + '/failed_batches.txt', 'wb')
    writer = csv.writer(failed_batches, delimiter = '\n')
    writer.writerows(failed_batch)
    


if __name__== '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', required = True, \
                            help = 'A text file with a list of SNPs')
    parser.add_argument('-o', '--output', default = 'default', \
                            help = 'The filepath for output')
    args = parser.parse_args()
    
    data_dir = '/mnt/3dgenome/projects/tfad334/all-gwas/data/results' # the directory containing raw GTEx data
    output_fp = set_filepath(args.input, args.output)
    snps = parse_snps(args.input)
    query_data(snps, data_dir)
    
