import pandas as pd
import os, sys

# get parameters from command line
fastq_dir = "../data/fastqs/postprocess"
metadata_dir = "../data/metadata"
matrices_dir = "../data/matrices"
name_mapping_1 = "%s/samplesheet.csv"%(matrices_dir)
name_mapping_2 = "%s/FILION_02-name-mapping.txt"%(metadata_dir)

# latent cells were called differently in the two samples
sample_mapping = {'J-Lat+DMSO' : 'J-LatA2+DMSO',
                  'J-Lat+SAHA' : 'J-LatA2+SAHA',
                  'Jurkat'     : 'Jurkat+DMSO'}

# parse the name mapping file for FILION_01
mapping = {}
with open(name_mapping_1, 'r') as f :
    for i, line in enumerate(f) :
        # skip first line
        if i == 0 : continue
        # now parse
        name, status, label = line.strip('\n').split('\t')
        mapping[name] = {'name' : name,
                         'sample' : sample_mapping[label]}

# parse the name mapping file for FILION_02
with open(name_mapping_2, 'r') as f :
    for line in f :
        cell, fname, name, sample = line.strip('\n').split()
        mapping[fname] = {'name' : name,
                          'sample' : sample}

# dive in the fastq dir
df = pd.DataFrame()
for d in os.listdir(fastq_dir) :
    full_dir_name = "%s/%s"%(fastq_dir, d)
    if os.path.isdir(full_dir_name) :
        data = pd.read_csv("%s/quant.genes.sf"%(full_dir_name), sep="\t")
        df['Name'] = data['Name']
        df[mapping[d]['name']] = data['NumReads']

# save data frame
df = df.set_index('Name')
df.to_csv("%s/exprMatrix.tsv"%(matrices_dir), sep="\t", index_label='Name')

# now write the sample sheet
with open('%s/sampleSheet.tsv'%(metadata_dir), 'w') as f :
    for name in list(mapping.keys()) :
        f.write('%s\t%s\n'%(mapping[name]['name'], mapping[name]['sample']))
