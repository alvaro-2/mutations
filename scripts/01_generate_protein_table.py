import pandas as pd
import numpy

#uniprot_acc mappping all data (gene_name, gene_id, hgnc_id, etc) of proteins related to mlo
protein = pd.read_csv('../raw_data/uniprot_associated_to_mlo.tsv', sep="\t")

# Mobidb disorder content
mobidb = pd.read_csv('../raw_data/dc_mobidb_lite.csv').rename(columns={'uniprot': 'uniprot_acc', 'dc': 'disorder_content'})

#add disorder content
protein = protein.merge(mobidb, how= 'left')
protein['id_protein'] = range(1, len(protein)+1)

protein = protein[['id_protein', 'uniprot_acc', 'uniprot_status', 'length', 'uniprot_name', 'hgnc_id', 'gene_id', 'gene_name', 'disorder_content', 'gene_name_synonyms', 'protein_names', 'sequence']]

protein.to_csv('../db_tables/protein.tsv', sep='\t', index= False)
'''
print(numpy.nanmax(protein['uniprot_status'].str.len())) #10
print(numpy.nanmax(protein['sequence'].str.len())) #32759
print(numpy.nanmax(protein['hgnc_id'].str.len())) #72
print(numpy.nanmax(protein['gene_name'].str.len())) #75
print(numpy.nanmax(protein['gene_id'].str.len())) #74
print(numpy.nanmax(protein['gene_name_synonyms'].str.len())) #280
print(numpy.nanmax(protein['protein_names'].str.len())) #865
print(numpy.nanmax(protein['uniprot_name'].str.len())) #16
print(numpy.nanmax(protein['uniprot_acc'].str.len())) #10
'''