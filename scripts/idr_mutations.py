from os import read
import pandas as pd
import matplotlib.pyplot as plt
# To do: num de mutaciones en zonas desord vs ordenadas; las PTMs también
# 

protein = pd.read_csv(
    '../db_tables/protein.tsv',
    sep='\t',
    usecols=['id_protein', 'length', 'disorder_content']
)

'''
protein = pd.read_csv(
    '../raw_data/uniprot_all_data_associated_to_mlo_preproc.tsv',
    sep= "\t",
    usecols= ['id_protein', 'uniprot_acc', 'uniprot_other_accesions']
)

id_protein = protein[['id_protein', 'uniprot_other_accesions']].rename(columns={'uniprot_other_accesions': 'uniprot_acc'})
id_protein['uniprot_acc'] = id_protein['uniprot_acc'].str.split(';')
id_protein = id_protein.explode('uniprot_acc')
id_protein['uniprot_acc'] = id_protein['uniprot_acc'].str.strip()
id_protein = id_protein[~id_protein['uniprot_acc'].isnull()]
id_protein = id_protein[id_protein['uniprot_acc'] != '']
id_protein = pd.concat([protein[['id_protein', 'uniprot_acc']], id_protein])
'''

disorder_region = pd.read_csv('../db_tables/disorder_region.tsv', sep='\t')
mutation = pd.read_csv('../db_tables/mutation.tsv', sep='\t')
mutation_has_disorder = pd.read_csv('../db_tables/mutation_has_disorder_region.tsv', sep='\t')

# sumar los length de idrs
dis_aa = disorder_region.length.sum() # 574352 aa desordenados
# los de las protein total
total_aa = protein.length.sum() # 3718008 totales
# aa ordenados
ord_aa = total_aa - dis_aa # 3143656

# fraccion de desorden
dis_aa/total_aa #0.15

# mutaciones en idrs
mut_idr = len(mutation_has_disorder)

# mutaciones por idrs aa
mut_idr/dis_aa # 0.44
# mutaciones por aa totales
#mut_idr/total_aa # 0.06

# mutaciones en region No idr
mut_no_idr = len(mutation[~mutation.id_mutation.isin(mutation_has_disorder.id_mutation)])
mut_no_idr/total_aa # 0.37


fig, ax = plt.subplots(figsize=(6,5))
ax.bar(["ordered aa", "disordered aa"] , [mut_no_idr/total_aa, mut_idr/dis_aa], color=["lightgray", "lightblue"])
ax.set_ylabel("Proportion")
ax.set_title("Mutations content\nOrder vs Disorder", fontweight="bold")
plt.show()

# COnsiderar solo 1 mut por aa
# subset de una sola mutacion por aa
mut_aa_unique = mutation[['start_aa', 'end_aa', 'id_protein']].drop_duplicates() # esto es tanto ordenados como desordenados
# eliminar mutaciones en rangos
mut_aa_unique = mut_aa_unique[mut_aa_unique.end_aa - mut_aa_unique.start_aa == 0]

# Ahora sumar todos los aa en estas proteinas
mut_aa_unique.id_protein.unique() # 5402
total_unique_aa = protein[protein.id_protein.isin(mut_aa_unique.id_protein.unique())].length.sum() # 3514830

# Sacar el nro de ordenados
ord_unique_aa = total_unique_aa - dis_aa
# coeficiente
len(mut_aa_unique)/total_unique_aa #0.348


# ahora ver idrs
idr = mutation_has_disorder.merge(disorder_region)
idr[idr.id_protein.isin(mut_aa_unique.id_protein.unique())] # estan todas

# Plot
fig, ax = plt.subplots(figsize=(6,8))
ax.bar(
    ["ordered aa", "disordered aa"] ,
    [len(mut_aa_unique)/ord_unique_aa, mut_idr/dis_aa],
    color=["lightgray", "lightblue"],
    width= 0.5
)
ax.set_ylabel("Proportion mutations/aminoacids", fontweight= "bold")
ax.set_title("Mutations content\nOrder vs Disorder", fontweight="bold")
plt.show()

