# parse_uniprot.py

import pandas as pd
import numpy as np
import re

def seq3to1(seq):
    protein_letters_3to1 = {
        'Ala': 'A',
        'Cys': 'C',
        'Asp': 'D',
        'Glu': 'E',
        'Phe': 'F',
        'Gly': 'G',
        'His': 'H',
        'Ile': 'I',
        'Lys': 'K',
        'Leu': 'L',
        'Met': 'M',
        'Asn': 'N',
        'Pro': 'P',
        'Gln': 'Q',
        'Arg': 'R',
        'Ser': 'S',
        'Thr': 'T',
        'Val': 'V',
        'Trp': 'W',
        'Tyr': 'Y',
        'Asx': 'B',
        'Xaa': 'X',
        'Glx': 'Z',
        'Xle': 'J',
        'Sec': 'U',
        'Pyl': 'O',
        'Ter': '*'
    }
    seq2 = re.findall('[A-Z][a-z]{2}', seq)  
    return "".join(protein_letters_3to1.get(aa, "X") for aa in seq2)


# Processed proteins table with id_protein
protein = pd.read_csv('../db_tables/protein.tsv', sep='\t').drop(columns='disorder_content')
protein = protein[['id_protein', 'uniprot_acc', 'gene_name', 'gene_name_synonyms','length', 'sequence']]
humsavar = pd.read_csv('../raw_data/humsavar.tsv', sep='\t')

#len(protein.uniprot_acc.unique())
#protein.duplicated().any()
humsavar.info()
humsavar.change.isnull().any() # False

# How many unique proteins in humsavar dataset
len(humsavar.uniprot_acc.unique()) # 12891

# How many unique LLPS proteins in our dataset are in humsavar dataset
len((humsavar.uniprot_acc[humsavar.uniprot_acc.isin(protein.uniprot_acc)]).unique()) #2892

'''
# Cuantas de esas ya estan en clinvar?
mutations_clinvar = pd.read_csv('tablas_fer/db_tables/mutation.tsv', sep='\t')
mutations_clinvar.columns
clinvar =  mutations_clinvar[['id_mutation', 'snp_id', 'id_protein']].copy()
clinvar
len((clinvar.id_protein[clinvar.id_protein.isin(protein.id_protein)]).unique()) # 2499, casi todas ya estan en clinvar
'''

# Merge with LLPS human proteins
uniprot_llps = protein.merge(humsavar, left_on='uniprot_acc', right_on='uniprot_acc', how= 'left')
uniprot_llps.duplicated().any() # False
# Keep entries with non-null change col values
uniprot_llps = uniprot_llps[uniprot_llps.change.notnull()].copy()

# Same with isin() method
#humsavar[humsavar.uniprot_acc.isin(protein.uniprot_acc)]

len(uniprot_llps.uniprot_acc.unique()) # only 2892 LLPS proteins with protein change in UniProt humsavar dataset

uniprot_llps.gene_name_y.notnull().sum() # 21246
uniprot_llps.gene_name_x.notnull().sum() # 21246

#uniprot_llps.gene_name_y = uniprot_llps.gene_name_y.combine_first(uniprot_llps.gene_name_x) # keep this col
uniprot_llps.drop(columns= 'gene_name_y', inplace= True)
uniprot_llps.rename(columns= {'gene_name_x': 'gene_name'}, inplace= True)

# Separate mutations and positions in different cols
uniprot_llps[['id_protein', 'uniprot_acc', 'gene_name', 'snp_id', 'change', 'category', 'mim', 'disease']]

# Regex for extracting Missense mutations (as they are in uniprot's humsavar dataset)
uniprot_llps['aux'] = uniprot_llps.change[uniprot_llps.change.notnull()].map(lambda x: re.findall('^([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})$', x)).str[0]

uniprot_llps['from_aa'] = uniprot_llps.aux[uniprot_llps.aux.notnull()].map(lambda x: x[0])
uniprot_llps['to_aa'] = uniprot_llps.aux[uniprot_llps.aux.notnull()].map(lambda x: x[2])
uniprot_llps['start_aa'] = uniprot_llps.aux[uniprot_llps.aux.notnull()].map(lambda x: x[1])
uniprot_llps['end_aa'] = uniprot_llps.aux[uniprot_llps.aux.notnull()].map(lambda x: x[1])

uniprot_llps[['change', 'from_aa', 'to_aa',	'start_aa',	'end_aa']]

# Control: from_aa corresponding to sequence
#Paso los aa de tres letras a una
uniprot_llps['ctrl'] = False
#uniprot_llps['aa_ctrl'] = np.nan

for i in uniprot_llps.index:
    aa1 = uniprot_llps.from_aa[i]
    aa2 = uniprot_llps.to_aa[i]
    if str(aa1) != 'nan':
        uniprot_llps['from_aa'][i] = str(seq3to1(aa1))

        #Evaluo
        if int(uniprot_llps.start_aa[i]) <= int(uniprot_llps.length[i]): 
            if uniprot_llps.sequence[i][int(uniprot_llps.start_aa[i])-1] == uniprot_llps.from_aa[i]:
                uniprot_llps.ctrl[i] = True

    if str(aa2) != 'nan':
        uniprot_llps['to_aa'][i] = str(seq3to1(aa2))

uniprot_llps.ctrl.value_counts() # 280 False. 279 from TTN protein and 1 from MTHFD1 protein... Maybe isoforms?
uniprot_llps.drop(columns= ['aux', 'sequence', 'length', 'gene_name_synonyms'], inplace= True)

uniprot_llps[uniprot_llps.disease.notnull()] # 11926 entries with a disease annotated

# Add columns consequence and source
uniprot_llps['consequence'] = 'missense'
uniprot_llps['source'] = 'uniprot'

# Final table
uniprot_llps.info()
(uniprot_llps[uniprot_llps.disease.notnull()].disease.value_counts() > 10).sum() # 225
uniprot_llps[uniprot_llps.disease.notnull()].disease.value_counts()[:20]

# Save
#uniprot_llps.to_csv('datasets/uniprot_all_proteins_mutations.tsv.gz', sep='\t', index= False, compression= 'gzip')