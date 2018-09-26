import numpy as np
import os, gzip

# directories
sc_hiv_root = '%s/work/CRG/projects/sc_hiv'%(os.getenv('HOME'))
matrices_dir = '%s/data/matrices'%(sc_hiv_root)
out_dir = '%s/monocle'%(matrices_dir)

# prepare meta-data of the cells
labels = ['Jurkat']*6 + ['J-Lat+DMSO']*30 + ['J-Lat+SAHA']*60

# analyze the two samples in the matrices directory
sample_names = ['P2449', 'P2458']
for sample_name in sample_names :

    # get matrix file name
    matrix_fname = '%s/%s.tsv.gz'%(matrices_dir, sample_name)

    # variables init
    genes = []

    # open and parse the matrix
    with gzip.open(matrix_fname, 'r') as f :

        # read the file line by line
        for lineno, line in enumerate(f) :

            # tab-separated file
            curatedline = line.strip('\n').split('\t')

            # get header and store cell names
            if lineno==0 :
                cells = curatedline[1:]
            else :
                # get gene names
                genes.append(curatedline[0])

    # write the phenoData file
    phenodata_fname = '%s/%s.pd.tsv'%(out_dir, sample_name)
    with open(phenodata_fname, 'w') as fout :
        fout.write('cell\tlabel\n')
        for cell, label in zip(cells, labels) :
            fout.write('%s\t%s\n'%(cell, label))

    # write the featureData file
    featuredata_fname = '%s/%s.fd.tsv'%(out_dir, sample_name)
    with open(featuredata_fname, 'w') as fout :
        fout.write('gene\tgene_short_name\n')
        for gene in genes :
            fout.write('%s\t%s\n'%(gene, gene))#.split('.')[0]))
