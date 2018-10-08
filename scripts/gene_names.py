import os, gzip
from pyensembl import EnsemblRelease

# genome
release = EnsemblRelease(93)

# directories
sc_hiv_root = '%s/work/CRG/projects/sc_hiv'%(os.getenv('HOME'))
matrices_dir = '%s/data/matrices'%(sc_hiv_root)
out_dir = '%s'%(matrices_dir)
sample_name = 'P2449'

# get matrix file name and out file name
matrix_fname = '%s/%s.tsv.gz'%(matrices_dir, sample_name)
out_fname = '%s/gene_annotations.tsv'%(matrices_dir)

# open and parse the matrix
with gzip.open(matrix_fname, 'r') as fin, open(out_fname, 'w') as fout :

    # read the file line by line
    for lineno, line in enumerate(fin) :

        # tab-separated file
        curatedline = line.strip('\n').split('\t')

        # get header and store cell names
        if lineno==0 :
            continue
        else :
            # get gene names
            gene_id = curatedline[0]
            try :
                gene_name = release.gene_by_id(gene_id.split('.')[0]).gene_name
            except ValueError :
                print "Warning: gene %s not found"%(gene_id)
            fout.write('%s\t%s\n'%(gene_id, gene_name))
