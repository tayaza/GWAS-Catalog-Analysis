#usr/bin/python
import os 
import sqlite3
import csv



def createDB():
    conn = sqlite3.connect('./spatial_eqtls.db')
    conn.text_factory = str
    cur = conn.cursor()
    
    cur.execute("DROP TABLE IF EXISTS eqtls")
    cur.execute(""" 
        CREATE TABLE eqtls 
        (snp TEXT, snp_chr INTEGER, snp_locus INTEGER, 
            gene TEXT, gene_chr INTEGER, gene_start INTEGER, gene_end INTEGER,
            cell_line TEXT, interaction TEXT, 'snp-gene_distance' INTEGER, 
            tissue TEXT, pvalue REAL)
        """)

    cur.execute("DROP TABLE IF EXISTS summary")
    cur.execute(""" 
        CREATE TABLE summary 
        (snp TEXT, snp_chr INTEGER, snp_locus INTEGER, 
            gene TEXT, gene_chr INTEGER, gene_start INTEGER, gene_end INTEGER,
            tissue TEXT, pvalue REAL, qvalue REAL, cell_lines TEXT, 
            gtex_cis_p_threshold REAL, cis TEXT, 'snp-gene_distance' INTEGER, 
            expression REAL)
        """)
    conn.close()

def popDB(dir):
    conn = sqlite3.connect('./spatial_eqtls.db')
    conn.text_factory = str
    cur = conn.cursor()

    dirs = []
    all_eqtls = []
    if os.path.isdir(dir):
        for d in os.listdir(dir):
            if os.path.isdir(dir + '/' + d):
                dirs.append(dir + '/' + d)

                
    for d in dirs:
        eqtls = []
        dfile = open(d + '/eqtls.txt', 'r')
        ddata = csv.reader(dfile, delimiter = '\t')
        next(ddata, None)
        for row in ddata:
            eqtls.append(row)
            all_eqtls.append(row)
        print(d, len(eqtls))
    for rec in all_eqtls:
        cur.executemany("""INSERT INTO eqtls  VALUES 
            (?,?,?,?,?,?,?,?,?,?,?,?);""", (rec,))

    """
    for d in dirs:
        eqtls = []
        dfile = open(d + '/summary.txt', 'r')
        ddata = csv.reader(dfile, delimiter = '\t')
        next(ddata, None)
        for row in ddata:
            if float(row[9]) <= 0.05: # Get sig eqtl interactions
                eqtls.append(row[:15])
                all_eqtls.append(row[:15])
        print(d, len(eqtls))
    for rec in all_eqtls:
        cur.executemany("""#INSERT INTO eqtls  VALUES 
            #(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?);
            #""", (rec,))
    
    conn.commit()
    conn.close()



if __name__ == "__main__":
    
    dir = '../data/results'
    createDB()
    popDB(dir)

