import pandas as pd
import os, sys

# check for proper invocation
if len(sys.argv) < 2 :
    print("Usage: make_expression_matrix.py <fastq_dir>", file=sys.stderr)
    sys.exit(1)

# get parameters from command line
fastq_dir = sys.argv[1]

# dive in the fastq dir
df = pd.DataFrame()
for d in os.listdir(fastq_dir) :
    full_dir_name = "%s/%s"%(fastq_dir, d)
    if os.path.isdir(full_dir_name) :
        data = pd.read_csv("%s/quant.genes.sf"%(full_dir_name), sep="\t")
        df['Name'] = data['Name']
        df[d] = data['NumReads']

# save data frame
df = df.set_index('Name')
df.to_csv("%s/exprMatrix.tsv"%(fastq_dir), sep="\t", index_label='Name')
