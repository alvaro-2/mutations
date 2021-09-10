import pandas as pd
import pyranges as pr


# Our database protein table
protein = pd.read_csv('datasets/protein.tsv', sep='\t')
protein_id = protein[['id_protein',	'uniprot_acc']].copy()

# Uniprot mappings
uniprot_data = pd.read_csv('datasets/uniprot_all_data_associated_to_mlo_preproc.tsv.txt', sep= '\t')
uniprot_data = uniprot_data[['id_protein', 'uniprot_acc', 'uniprot_other_accesions']]

# Mutations table
mutations = pd.read_csv('datasets/mutation.tsv', sep= '\t')

# PTMs from phosphositeplus 
ptms = pd.read_csv('ptms.csv').rename(columns= {'uniprot': 'uniprot_acc'})
#ptms = ptms.merge(uniprot_data).drop(columns='uniprot_acc') # se pierden ptms asi (seguro algunos uniprots quedaron obsoletos)
ptms = ptms.merge(protein_id).drop(columns='uniprot_acc')

ptms.duplicated().any() # True 1864
ptms.drop_duplicates(inplace= True)

# To use pyranges: change cols names (this is for mapping the mutations)
#mutations.rename(columns={'id_protein': 'Chromosome', 'start_aa': 'Start',	'end_aa': 'End'}, inplace= True)
# Create the pyranges object of mutations
#df_py = pr.PyRanges(mutations)

# Type table
type_ptm = pd.DataFrame(ptms.type.unique(), columns=['type'])
type_ptm["id_type"] = range(1, len(ptms.type.unique())+1)

# Class table
class_ptm = pd.DataFrame(ptms["clase"].unique(), columns=['class'])
class_ptm["id_class"] = range(1, len(ptms["clase"].unique())+1)

# ptm table
ptm = ptms.merge(type_ptm)
ptm = ptm.merge(class_ptm, left_on='clase', right_on='class').drop(columns= ['type', 'clase', 'class']).rename(columns= {'pos': 'poas_aa'})