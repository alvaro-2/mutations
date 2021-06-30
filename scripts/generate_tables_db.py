#%pip install cython
#%pip install pyranges


import pandas as pd
import numpy as np
import pyranges as pr


def format_snp(df, column):
    '''
    format an int snps column in a DataFrame containing -1 values.
    Returns: the snp column in str format ('rs1580653772' or 'nan')
    '''
    #a = df.column.replace(-1, 'nan')
    a = df[column]
    #a = a.apply(str)
    a = a.map(lambda x: 'rs' + str(int(x)) if not np.isnan(x) else x)
    df[column] = a


# Proteins from each LLPS database with their roles, mlos and dataset
# Each combination: uniprot, mlo, db, rol must be unique
database_entrada = pd.read_csv('database_entrada.csv')
database_entrada.drop(columns='organism', inplace= True)
database_entrada.drop_duplicates(inplace = True)
database_entrada.rename(columns= {'uniprot': 'uniprot_acc'}, inplace= True)
database_entrada.info()

# Drop rows with same uniprot, mlo, db and different rol
aa = database_entrada.drop_duplicates().groupby(['uniprot_acc', 'mlo', 'db']).size().reset_index(name='counts')
aa[aa.counts > 1]

# Check one of them
database_entrada[(database_entrada.uniprot_acc == 'O43781') & (database_entrada.mlo == 'Stress granule') & (database_entrada.db == 'phasepdb_rev')]
database_entrada[(database_entrada.uniprot_acc == 'P31483') & (database_entrada.mlo == 'Stress granule') & (database_entrada.db == 'phasepdb_rev')]
database_entrada[(database_entrada.uniprot_acc == 'Q92973') & (database_entrada.mlo == 'Stress granule') & (database_entrada.db == 'phasepdb_rev')]
database_entrada[(database_entrada.uniprot_acc == 'Q9UHD9') & (database_entrada.mlo == 'Stress granule') & (database_entrada.db == 'phasepdb_rev')]

database_entrada.drop([137, 141, 277, 204], inplace= True)

# ## Deal with mlos annotations

database_entrada.mlo = database_entrada.mlo.str.strip()
#len(database_entrada.mlo.unique()) # no blank spaces
#database_entrada.mlo.value_counts()

# ### Paraspeckle
(database_entrada.mlo == 'Paraspeckle').sum()
# Unify paraspeckle with Paraspeckle
database_entrada.replace('paraspeckle', 'Paraspeckle', inplace= True)
(database_entrada.mlo == 'Paraspeckle').sum()

# ### Sam68
(database_entrada.mlo == 'Sam68 nuclear bodies').sum()
(database_entrada.mlo == 'Sam68 nuclear bodies (SNBs)').sum()
(database_entrada.mlo == 'Sam68 nuclear body').sum()
database_entrada.replace(['Sam68 nuclear bodies', 'Sam68 nuclear bodies (SNBs)'], 'Sam68 nuclear body', inplace= True)
(database_entrada.mlo == 'Sam68 nuclear body').sum()

# ### PML body  
# **PhaSepDB**: The PML bodies are dynamic nuclear protein aggregates interspersed between chromatin. These punctate nuclear structures are call PML bodies because the PML gene is essential for their formation. are present in most mammalian cell nuclei and typically number 1 to 30 bodies per nucleus.  
# **DrLLPS**: PML nuclear bodies are annotetad in the nucleus. They are matrix-associated domains that recruit an astonishing variety of seemingly unrelated proteins.

(database_entrada.mlo == 'PML nuclear body').sum()
(database_entrada.mlo == 'PML body').sum()
database_entrada.replace('PML body', 'PML nuclear body', inplace= True)
(database_entrada.mlo == 'PML nuclear body').sum()

# ### Polycomb body
(database_entrada.mlo == 'Polycomb bodies').sum()
database_entrada.replace('Polycomb bodies', 'Polycomb body', inplace= True)
(database_entrada.mlo == 'Polycomb body').sum()

# ### Pre and postsynaptic density
database_entrada.replace('Pre and postsynaptic densities', 'Pre and postsynaptic density', inplace= True)
(database_entrada.mlo == 'Pre and postsynaptic density').sum()

# ### Nuclear speckle
(database_entrada.mlo == 'Nucleus speckles').sum() #phasepdb
(database_entrada.mlo == 'Nuclear speckle').sum() # drllps
(database_entrada.mlo == 'Nuclear speckles').sum() #phasepdb
(database_entrada.mlo == 'nuclear speckle').sum()
database_entrada.replace(['Nucleus speckles', 'Nuclear speckles', 'nuclear speckle'], 'Nuclear speckle', inplace= True)
(database_entrada.mlo == 'Nuclear speckle').sum()

# ### Heterochromatin
(database_entrada.mlo == 'heterochromatin').sum()
database_entrada.replace('heterochromatin', 'Heterochromatin', inplace= True)
(database_entrada.mlo == 'Heterochromatin').sum()

# ### Cytoplasmic ribonucleoprotein granule
(database_entrada.mlo == 'cytoplasmic ribonucleoprotein granule').sum()
database_entrada.replace('cytoplasmic ribonucleoprotein granule', 'Cytoplasmic ribonucleoprotein granule', inplace= True)

# ### Membrane cluster
(database_entrada.mlo == 'Membrane clusters').sum()
(database_entrada.mlo == 'membrane cluster').sum()
database_entrada.replace(['Membrane clusters', 'membrane cluster'], 'Membrane cluster', inplace= True)
(database_entrada.mlo == 'Membrane cluster').sum()

# ### Nuclear body
database_entrada.replace('nuclear body', 'Nuclear body', inplace= True)
(database_entrada.mlo == 'Nuclear body').sum()

# ### Nucleolus
database_entrada.replace('nucleolus', 'Nucleolus', inplace= True)
(database_entrada.mlo == 'Nucleolus').sum()
(database_entrada.mlo == 'Centrosome/Spindle pole body').sum() # keep this annotation

# EXPLODE:
# P-body, Stress granule
# P-body, GW body
# Set mlo col as list-like and explode() to separate list elements into separate rows
# before: 8178 rows
database_entrada = database_entrada.assign(mlo= database_entrada.mlo.str.split(',')).explode('mlo')
database_entrada.mlo = database_entrada.mlo.str.strip()
database_entrada.drop_duplicates(inplace= True)
database_entrada
# after: 8183 rows

# GW-body
database_entrada.replace('GW body', 'GW-body', inplace= True)
(database_entrada.mlo == 'GW-body').sum()

# Postsynaptic density
database_entrada.replace('postsynaptic density', 'Postsynaptic density', inplace= True)
(database_entrada.mlo == 'Postsynaptic density').sum()

# Cytoplasmic ribonucleoprotein granule
database_entrada.replace('cytoplasmic ribonucleoprotein granule', 'Cytoplasmic ribonucleoprotein granule', inplace= True)
(database_entrada.mlo == 'Cytoplasmic ribonucleoprotein granule').sum()

# Histone locus body
database_entrada.replace('Histone Locus body', 'Histone locus body', inplace= True)
(database_entrada.mlo == 'Histone locus body').sum()

# Stress granule
database_entrada.replace('Sress granule', 'Stress granule', inplace= True)
(database_entrada.mlo == 'Stress granule').sum()

database_entrada.drop_duplicates(inplace= True)
database_entrada.info()


# Load proteins, mutations, domains and regions tables

# protein table for our db. Same above but one protein by row
protein = pd.read_csv('db_tables/protein.tsv', sep='\t')

# DataFrame with unique id_protein col
id_protein = protein[['id_protein', 'uniprot_acc']].copy()

# only ClinVar mutations at the moment
mutations = pd.read_csv('../datasets/mutations.tsv.gz', sep='\t', compression='gzip') # comes from parse_clinvar.py
# ClinVar mutations with source dataset
mutations_with_source = pd.read_csv('../datasets/mutations_with_source.tsv.gz', sep='\t', compression='gzip')
# ClinVar PMIDs
pmid = pd.read_csv('var_citations.txt', sep='\t')

# Our llps' proteins domains and regions
disorder = pd.read_csv('disorder_lite.csv').rename(columns={'uniprot': 'uniprot_acc'})
low_complexity = pd.read_csv('low_complexity.csv').rename(columns={'uniprot': 'uniprot_acc'})

pfam = pd.read_csv('pfam.csv').rename(columns={'uniprot': 'uniprot_acc', 'tipo': 'pfam_name'})
pfam_map = pd.read_csv('pfam_map.csv')
pfam_map.duplicated().any()
pfam = pfam.merge(pfam_map)
pfam.duplicated().any()

# Add and unique integer ID fow low_complexity and disorder
low_complexity['id_lc'] = range(1, len(low_complexity)+1)
disorder['id_idr'] = range(1, len(disorder)+1)

##################### Generate the tables #########################

# consequence table
cf = mutations.consequence.value_counts()
consequence = pd.DataFrame({'id_consequence': range(1, len(cf)+1), 'consequence': cf.index})
consequence.to_csv('db_tables/consequence.tsv', sep='\t', index = False)

# mutation table  
# cols: *id_mutation, snp_id, chromosome, start_genomic, end_genomic, start_aa, end_aa, from_aa, to_aa, id_source, id_protein, nt_change*

# Subset by cols to keep for mutation db table
mutation = mutations[['id_protein', 'id_mutation', 'snpid', 'chromosome', 'start', 'stop', 'start_aa', 'end_aa', 'from', 'to', 'consequence', 'cambio_nt']].copy()
mutation.rename(columns={'snpid': 'snp_id', 'start': 'start_genomic', 'stop': 'end_genomic', 'from': 'from_aa', 'to': 'to_aa', 'cambio_nt': 'nt_change'}, inplace= True)

# %%
# Add IDs from consequence
mutation = mutation.merge(consequence)
mutation.drop(columns='consequence', inplace= True)

# Format snp_id col
format_snp(mutation, 'snp_id')

mutation = mutation[['id_mutation', 'snp_id', 'chromosome', 'start_genomic', 'end_genomic', 'start_aa','end_aa',
                    'from_aa', 'to_aa', 'id_protein', 'id_consequence', 'nt_change']].sort_values('id_mutation')

mutation.chromosome = mutation.chromosome.apply(str)
#mutation.duplicated().any()
#mutation.nt_change.str.len().max()
mutation.to_csv('db_tables/mutation.tsv', sep='\t', index = False)

# Para asignar los rangos debo tener:  
# Tabla de mutaciones con id_mutation, *id_protein(Chromosome), start_aa(Start), end_aa(End)*  
# Tablas de lc, idr y pfam con id unico

# Pfam Tables
#len(pfam.pfam_name.unique())
#len(pfam.pfam_acc.unique())
#pfam[['pfam_name', 'pfam_acc']].drop_duplicates()

# ## pfam_domain  
# cols: pfam_id, pfam_domain, por ej: PF00003 7tm_3

# Array with unique pfam domains
pf_domain = pfam.pfam_name.unique() # unique pfam domains (2939 for this set of proteins)
pfam_domain = pfam[['pfam_name', 'pfam_acc']].drop_duplicates()
pfam_domain.rename(columns={'pfam_acc': 'id_pfam', 'pfam_name': 'pfam_domain'}, inplace= True)
pfam_domain.to_csv('db_tables/pfam_domain.tsv', sep='\t', index= False)

# ## protein_has_pfam_domain  
# cols: id_protein, id_pfam, start, end, length

protein_has_pfam_domain = pfam.merge(id_protein) # agregar col id_protein
protein_has_pfam_domain['length'] = protein_has_pfam_domain.end - protein_has_pfam_domain.start + 1 # col length
protein_has_pfam_domain.drop(columns='pfam_name', inplace= True)
protein_has_pfam_domain = protein_has_pfam_domain.merge(pfam) # to add the col pfam_id
protein_has_pfam_domain.rename(columns={'pfam_acc': 'id_pfam'}, inplace= True)
protein_has_pfam_domain = protein_has_pfam_domain[['id_protein', 'id_pfam', 'start', 'end', 'length']].sort_values('id_protein')

#protein_has_pfam_domain.duplicated().any()
protein_has_pfam_domain.to_csv('db_tables/protein_has_pfam_domain.tsv', sep='\t', index= False)

# ## mutation_has_pfam_domain  
# cols: id_mutation, id_protein, id_pfam, start, end

# ### Pyranges  
# columnas obligatorias: *Chromosome	 Start	End*  
# Chromosome: id_protein    
# otras columnas con ids son opcionales y cualquier nombre  
#   
# por ejemplo df seria tabla de mutaciones  
# df = pr.PyRanges(df.rename(columns={'chromosome':'Chromosome','start_position':'Start','end_position':'End'}))  
#   
# df = pyrange de mutaciones (columnas: Chromosome, Start, End, id_mutacion)  
# low_c = pyrange de low complexity(columnas: Chromosome, Start, End, id_low, id_proteina)  
# data = df.join(low_c, strandedness=False, slack=1).drop(like="_b") # mutaciones lo junto con low_complex  
# strandedness=False no tener en cuenta el Strand  
# slack=1 coincidir los extremos. Importante  
# drop(like="_b") eliminar el Chromosome, Start, End de low_c (en pfam no hacer el drop)  
# data = data.df[[Chromosome, Start, End, id_mutacion, id_low, id_proteina]] # pasa de pyrange a dataframe

# df has pfam domains data
df = pfam.rename(columns={'pfam_name': 'pfam_domain'}).merge(pfam_domain)
df = df.merge(id_protein) # mapping uniprot_acc - id_protein
df.drop(columns='uniprot_acc', inplace= True)
df.rename(columns={'id_protein': 'Chromosome', 'start': 'Start', 'end': 'End'}, inplace= True)

# Create the pyranges object of pfam domains
df_py = pr.PyRanges(df)

aux = mutation[['start_aa', 'end_aa', 'id_mutation', 'id_protein']].copy()
aux.rename(columns={'id_protein': 'Chromosome', 'start_aa': 'Start', 'end_aa': 'End'}, inplace= True)

# Pyranges object of mutations
aux_py = pr.PyRanges(aux)

# Join both pyranges object: this assings mutations to pfam domains
pfam_py = df_py.join(aux_py, strandedness= False, slack= 1)  # strandedness= False doesnt take count of the chain strand; slack= 1 include bound
#pfam_py.head() # Start and End are from the pfam domain in that protein (a protein may have the same pfam domain repeated at different positions along its sequence).
                # Start_b and End_b are from the mutation in this case

# Pyranges to DataFrame
mutation_has_pfam_domain = pfam_py.df[['id_mutation', 'Chromosome', 'id_pfam', 'Start', 'End']] # cols to keep
mutation_has_pfam_domain.rename(columns={'Chromosome': 'id_protein', 'Start': 'start', 'End': 'end'}, inplace= True)

#mutation_has_pfam_domain.head()

# control
#mutation[mutation.id_mutation == 23987]
# control
#pfam_domain[pfam_domain.id_pfam == 814] 
#mutation_has_pfam_domain[mutation_has_pfam_domain.id_pfam == 814] # ok!

# Save
mutation_has_pfam_domain.to_csv('db_tables/mutation_has_pfam_domain.tsv', sep='\t', index= False)

# # low-complexity Tables

# ## low_complexity  
# cols: id_lc, start, end, length, id_protein

# Add length col 
low_complexity['length'] = low_complexity.end - low_complexity.start + 1 
# Add id_proteins
low_complexity.rename(columns={'uniprot': 'uniprot_acc'}, inplace= True)
low_complexity = low_complexity.merge(id_protein)
low_complexity.drop(columns='uniprot_acc', inplace= True)
# Save
low_complexity.to_csv('db_tables/low_complexity.tsv', sep='\t', index= False)

# ## mutation_has_low_complexity  
# cols: id_mutation, id_lc

# Table for LC data
lc_has = low_complexity.copy()
lc_has.rename(columns={'id_protein': 'Chromosome', 'start': 'Start', 'end': 'End'}, inplace= True)

# Auxiliar table for mutations
aux_lc = mutation[['start_aa', 'end_aa', 'id_mutation', 'id_protein']].copy()
aux_lc.rename(columns={'id_protein': 'Chromosome', 'start_aa': 'Start', 'end_aa': 'End'}, inplace= True)

# Create the Pyranges objects
lc_has_py = pr.PyRanges(lc_has)
aux_lc_py = pr.PyRanges(aux_lc)

# Join both pyranges object: this assings mutations to low-complexity regions
lc_py = aux_lc_py.join(lc_has_py, strandedness= False, slack=1).drop(like="_b") # strandedness= False doesnt take count of the chain strand;
                                                                                # slack= 1 include bounds; drop(like="_b"): delete those cols (redudants)

# Pyrange to DataFrame
mutation_has_low_complexity = lc_py.df[['id_mutation', 'id_lc']] # cols to keep

# Control
#low_complexity[low_complexity.id_lc == 5240]
#protein.iloc[8]
#mutation[mutation.id_mutation == 23989] # It's allright! Mutation in aa 1133, which belongs to the low-complexity region between 1128 - 1139 in that protein

# Save
mutation_has_low_complexity.to_csv('db_tables/mutation_has_low_complexity.tsv', sep='\t', index= False)

# # Disorder Tables

# ## disorder_region  
# cols: id_idr, start, end, length, id_protein

# Add length col 
disorder['length'] = disorder.end - disorder.start + 1 
disorder_region = disorder.rename(columns={'uniprot': 'uniprot_acc'}).merge(id_protein).sort_values('id_protein')
disorder_region.drop(columns='uniprot_acc', inplace= True)
# Save
disorder_region.to_csv('db_tables/disorder_region.tsv', sep='\t', index= False)

# ## mutation_has_disorder_region  
# cols: id_mutation, id_idr

# Auxiliar table for mutations from low-complexity is the same for disorder. id-protein, start and end of the mutation
aux_idr = aux_lc

# Table for IDRs data
idr_has = disorder_region.copy()
idr_has.rename(columns={'id_protein': 'Chromosome', 'start': 'Start', 'end': 'End'}, inplace= True)

# Create the Pyranges objects
idr_has_py = pr.PyRanges(idr_has)
aux_idr_py = pr.PyRanges(aux_idr)

# Join both pyranges object: this assings mutations to pfam domains
idr_py = aux_idr_py.join(idr_has_py, strandedness= False, slack=1).drop(like="_b") # strandedness= False doesnt take count of the chain strand;
                                                                                   # slack= 1 include bounds; drop(like="_b"): delete those cols (redudants)

# Pyrange to DataFrame
mutation_has_disorder_region = idr_py.df[['id_mutation', 'id_idr']] # cols to keep

# Control
#mutation[mutation.id_mutation == 24001]
#id_protein[id_protein.id_protein == 12]
#disorder[disorder.id_idr == 4339] # It's Ok. A point mutation in position 70 in the idr region between 48-75

# Save
mutation_has_disorder_region.to_csv('db_tables/mutation_has_disorder_region.tsv', sep='\t', index= False)
 
# # Rol table  
# cols: id_rol, rol

#database_entrada.rol.unique()
#database_entrada.rol.value_counts()

rol = pd.DataFrame({'rol': database_entrada.rol.value_counts().index, 'id_rol': range(1, len(database_entrada.rol.value_counts())+1)})
# Save
rol.to_csv('db_tables/rol.tsv', sep='\t', index= False)

# # dataset table  
# cols: id_dataset, dataset
#database_entrada.db.value_counts()

dataset = pd.DataFrame({'dataset': database_entrada.db.value_counts().index, 'id_dataset': range(1, len(database_entrada.db.value_counts())+1)})
# Save
dataset.to_csv('db_tables/dataset.tsv', sep='\t', index= False)

# # MLOs tables  
# cols: id_mlo, mlo
#database_entrada.info()
#database_entrada.mlo.value_counts()
#len(database_entrada.uniprot_acc.unique()) # OK

mlo = pd.DataFrame({'mlo': database_entrada.mlo.value_counts().index, 'id_mlo': range(1, len(database_entrada.mlo[database_entrada.mlo.notnull()].unique())+1)})
# Save
mlo.to_csv('db_tables/mlo.tsv', sep='\t', index= False)

# ## protein_has_mlo  
# cols: id_protein, id_mlo, id_rol, id_database

#len(database_entrada.uniprot_acc.unique())
protein_has_mlo = database_entrada.copy()

# Add id_protein
protein_has_mlo = protein_has_mlo.merge(id_protein)
# Add id_mlo
protein_has_mlo = protein_has_mlo.merge(mlo, how= 'left')
# Add id_rol and id_database
protein_has_mlo = protein_has_mlo.merge(rol)
protein_has_mlo = protein_has_mlo.rename(columns={'db': 'dataset'}).merge(dataset).sort_values('id_protein')
protein_has_mlo.drop(columns=['uniprot_acc', 'mlo', 'rol', 'dataset'], inplace= True)

protein_has_mlo[protein_has_mlo.duplicated()] # OK
protein_has_mlo['id_proteinmlo'] = range(1, len(protein_has_mlo)+1)

#protein_has_mlo.duplicated().any()
# Save
protein_has_mlo.to_csv('db_tables/protein_has_mlo.tsv', sep='\t', index= False)

# # source
source = pd.DataFrame({'id_source': [1,2,3], 'source': ['clinvar', 'disgenet', 'uniprot']})
# Save
source.to_csv('db_tables/source.tsv', sep='\t', index = False)


# # mutation_has_source
mutation_has_source = mutations_with_source.copy()
mutation_has_source.rename(columns={'variationid': 'id_insource'}, inplace= True)
mutation_has_source.drop_duplicates(inplace= True)
mutation_has_source = mutation_has_source.merge(source).drop(columns='source')
# Save
mutation_has_source.to_csv('db_tables/mutation_has_source.tsv', sep='\t', index= False)


# # citation_source
pmid.columns = pmid.columns.str.lower().str.replace(' ',"_").str.replace("-",'_').str.replace('/','_')
pmid = pmid[['variationid', 'citation_source', 'citation_id']].copy()

#pmid.variationid.isnull().any() # False. This is the ClinVar ID
pmid.variationid.drop_duplicates(inplace= True)

citation_source = pd.DataFrame({'name': pmid.citation_source.value_counts().index, 'id_citation_source': range(1, len(pmid.citation_source.unique())+1) })
# Save
citation_source.to_csv('db_tables/citation_source.tsv', sep='\t', index= False)


# # mutation_has_citation
mutation_has_citation = pmid.rename(columns={'citation_id': 'id_citation'})
mutation_has_citation = mutation_has_citation.merge(citation_source.rename(columns={'name': 'citation_source'})).drop(columns= 'citation_source')
mutation_has_citation = mutation_has_citation.merge(mutations_with_source).drop(columns=['source', 'variationid']).drop_duplicates()

# Save
mutation_has_citation.to_csv('db_tables/mutation_has_citation.tsv', sep='\t', index= False)