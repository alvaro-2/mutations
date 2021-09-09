import pandas as pd

# Proteins from each LLPS database with their roles, mlos and dataset
# Each combination: uniprot, mlo, db, rol must be unique

#uniprot	 mlo db	rol	reviewed

mlo_db_rol = pd.read_csv('../raw_data/tablas_disphase_30-08/dataset_entrada_mlos.csv')
mlo_db_rol.drop_duplicates(inplace = True)
mlo_db_rol.rename(columns= {'uniprot': 'uniprot_acc'}, inplace= True)
#mlo_db_rol.info()

# ## Deal with mlos annotations

mlo_db_rol.mlo = mlo_db_rol.mlo.str.strip()
#len(database_entrada.mlo.unique()) # no blank spaces
#database_entrada.mlo.value_counts()

# ### Paraspeckle
(mlo_db_rol.mlo == 'Paraspeckle').sum()
# Unify paraspeckle with Paraspeckle
mlo_db_rol.replace('paraspeckle', 'Paraspeckle', inplace= True)
(mlo_db_rol.mlo == 'Paraspeckle').sum()

# ### Sam68
(mlo_db_rol.mlo == 'Sam68 nuclear bodies').sum()
(mlo_db_rol.mlo == 'Sam68 nuclear bodies (SNBs)').sum()
(mlo_db_rol.mlo == 'Sam68 nuclear body').sum()
mlo_db_rol.replace(['Sam68 nuclear bodies', 'Sam68 nuclear bodies (SNBs)'], 'Sam68 nuclear body', inplace= True)
(mlo_db_rol.mlo == 'Sam68 nuclear body').sum()

# ### PML body  
# **PhaSepDB**: The PML bodies are dynamic nuclear protein aggregates interspersed between chromatin. These punctate nuclear structures are call PML bodies because the PML gene is essential for their formation. are present in most mammalian cell nuclei and typically number 1 to 30 bodies per nucleus.  
# **DrLLPS**: PML nuclear bodies are annotetad in the nucleus. They are matrix-associated domains that recruit an astonishing variety of seemingly unrelated proteins.

(mlo_db_rol.mlo == 'PML nuclear body').sum()
(mlo_db_rol.mlo == 'PML body').sum()
mlo_db_rol.replace('PML body', 'PML nuclear body', inplace= True)
(mlo_db_rol.mlo == 'PML nuclear body').sum()

# ### Polycomb body
(mlo_db_rol.mlo == 'Polycomb bodies').sum()
mlo_db_rol.replace('Polycomb bodies', 'Polycomb body', inplace= True)
(mlo_db_rol.mlo == 'Polycomb body').sum()

# ### Pre and postsynaptic density
mlo_db_rol.replace('Pre and postsynaptic densities', 'Pre and postsynaptic density', inplace= True)
(mlo_db_rol.mlo == 'Pre and postsynaptic density').sum()

# ### Nuclear speckle
(mlo_db_rol.mlo == 'Nucleus speckles').sum() #phasepdb
(mlo_db_rol.mlo == 'Nuclear speckle').sum() # drllps
(mlo_db_rol.mlo == 'Nuclear speckles').sum() #phasepdb
(mlo_db_rol.mlo == 'nuclear speckle').sum()
mlo_db_rol.replace(['Nucleus speckles', 'Nuclear speckles', 'nuclear speckle'], 'Nuclear speckle', inplace= True)
(mlo_db_rol.mlo == 'Nuclear speckle').sum()

# ### Heterochromatin
(mlo_db_rol.mlo == 'heterochromatin').sum()
mlo_db_rol.replace('heterochromatin', 'Heterochromatin', inplace= True)
(mlo_db_rol.mlo == 'Heterochromatin').sum()

# ### Cytoplasmic ribonucleoprotein granule
(mlo_db_rol.mlo == 'cytoplasmic ribonucleoprotein granule').sum()
mlo_db_rol.replace('cytoplasmic ribonucleoprotein granule', 'Cytoplasmic ribonucleoprotein granule', inplace= True)

# ### Membrane cluster
(mlo_db_rol.mlo == 'Membrane clusters').sum()
(mlo_db_rol.mlo == 'membrane cluster').sum()
mlo_db_rol.replace(['Membrane clusters', 'membrane cluster'], 'Membrane cluster', inplace= True)
(mlo_db_rol.mlo == 'Membrane cluster').sum()

# ### Nuclear body
mlo_db_rol.replace('nuclear body', 'Nuclear body', inplace= True)
(mlo_db_rol.mlo == 'Nuclear body').sum()

# ### Nucleolus
mlo_db_rol.replace('nucleolus', 'Nucleolus', inplace= True)
(mlo_db_rol.mlo == 'Nucleolus').sum()
(mlo_db_rol.mlo == 'Centrosome/Spindle pole body').sum() # keep this annotation

# EXPLODE:
# P-body, Stress granule
# P-body, GW body
# Set mlo col as list-like and explode() to separate list elements into separate rows
# before: 8178 rows
mlo_db_rol = mlo_db_rol.assign(mlo=mlo_db_rol.mlo.str.split(',')).explode('mlo')
mlo_db_rol.mlo = mlo_db_rol.mlo.str.strip()
mlo_db_rol.drop_duplicates(inplace= True)
# after: 8183 rows

# GW-body
mlo_db_rol.replace('GW body', 'GW-body', inplace= True)
(mlo_db_rol.mlo == 'GW-body').sum()

# Postsynaptic density
mlo_db_rol.replace('postsynaptic density', 'Postsynaptic density', inplace= True)
(mlo_db_rol.mlo == 'Postsynaptic density').sum()

# Cytoplasmic ribonucleoprotein granule
mlo_db_rol.replace('cytoplasmic ribonucleoprotein granule', 'Cytoplasmic ribonucleoprotein granule', inplace= True)
(mlo_db_rol.mlo == 'Cytoplasmic ribonucleoprotein granule').sum()

# Histone locus body
mlo_db_rol.replace('Histone Locus body', 'Histone locus body', inplace= True)
(mlo_db_rol.mlo == 'Histone locus body').sum()

# Stress granule
mlo_db_rol.replace('Sress granule', 'Stress granule', inplace= True)
(mlo_db_rol.mlo == 'Stress granule').sum()

mlo_db_rol.drop_duplicates(inplace= True)

print(f'MLOs database and rol found {mlo_db_rol.shape[0]} rows')
print(f'unique MLOs {len(set(mlo_db_rol["mlo"].tolist()))}')
#mlo_db_rol.info()
mlo_db_rol.to_csv('../raw_data/mlo_db_rol_cleaned.tsv', sep="\t", index=False)
