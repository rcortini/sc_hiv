import pandas as pd

# load the sample sheet
datadir = "../data/"
fastqs_dir = "%s/fastq"
samplesheet = pd.read_excel('%s/metadata/FILION_02.xls'%(datadir), skiprows=2)

# prepare the names of the cells in the 96 well plate
alpha = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
cells = []
for i in range(1,13) :
     for a in alpha :
         cells.append('%s%d'%(a, i))

# extract the names of the multiplex indices and of the sample names
multiplex_indices = samplesheet['MULTIPLEX INDEX']
sample_names = samplesheet['SAMPLE NAME']
file_names = samplesheet['Unnamed: 9']

# mapping of the plate number to the external plate number
plates = {'P1' : 'P2769',
          'P2' : 'P2770',
          'P3' : 'P2771'}

# iterate on both at the same time, and make some horrendous parsing based on
# the specific characteristics of what is written in the sample sheet
# with open("../data/matrices/samplesheet_2.tsv", "w") as f :
    # f.write("cellnames\tlabel\n")
i = 0
for index, sample_name, file_name in zip(multiplex_indices, sample_names,
                                         file_names) :
    plate_id = sample_name[:2]
    name = plates[plate_id] + "_" + index[:9]
    sample_description = sample_name.split(':')[1].lstrip(' ')
    cell_type = sample_description.split(',')[0]
    treatment = sample_description.split(' ')[1]
    cell_id = cell_type + "+" + treatment
    print("%s %s %s %s"%(cells[i%96], file_name, name, cell_id))
    i += 1
