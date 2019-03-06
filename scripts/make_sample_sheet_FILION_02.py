import pandas as pd

# load the sample sheet
samplesheet = pd.read_excel('../data/FILION_02/FILION_02.xls', skiprows=2)

# extract the names of the multiplex indices and of the sample names
multiplex_indices = samplesheet['MULTIPLEX INDEX']
sample_names = samplesheet['SAMPLE NAME']

# mapping of the plate number to the external plate number
plates = {'P1' : 'P2769',
          'P2' : 'P2770',
          'P3' : 'P2771'}

# iterate on both at the same time, and make some horrendous parsing based on
# the specific characteristics of what is written in the sample sheet
with open("../data/matrices/samplesheet_2.tsv", "w") as f :
    f.write("cellnames\tlabel\n")
    for index, sample_name in zip(multiplex_indices, sample_names) :
        plate_id = sample_name[:2]
        name = plates[plate_id] + "_" + index[:9]
        sample_description = sample_name.split(':')[1].lstrip(' ')
        cell_type = sample_description.split(',')[0]
        treatment = sample_description.split(' ')[1]
        cell_id = cell_type + "+" + treatment
        f.write("%s\t%s\n"%(name, cell_id))
