import pandas as pd
import numpy as np


# ############################## LOAD DATASETS: ##############################    
# %% Disgenet data: variant - disease association curated; Summary of Curated VDAs
# this dataset comes from disgenet downloads (https://www.disgenet.org/downloads; Octubre 2020)
vda = pd.read_csv('../raw_data/curated_variant_disease_associations.tsv', sep='\t')
vda.columns = vda.columns.str.lower().str.replace(' ',"_").str.replace("-",'_').str.replace('/','_')
vda = vda.rename(columns={'snpid':'snp_id'})
# Add a generic id for each record
vda['id_disgenet'] = range(1, len(vda)+1)
# %%
# Cross references: Mappings UniProts. uniprots id with entrez gene id (https://www.disgenet.org/downloads, Mappings)
uniprots = pd.read_csv('../raw_data/mapa_geneid_4_uniprot_crossref.tsv.gz', sep='\t', compression='gzip').rename(columns={'UniProtKB':'uniprot_acc', 'GENEID':'gene_id'})

# The file contains the mappings of DisGeNET variants (dbSNP Identifiers)
# to NCBI Entrez identifiers according to dbSNP database (https://www.disgenet.org/downloads, Variant-Gene Mappings File)
variant_gene = pd.read_csv('../raw_data/variant_to_gene_mappings.tsv.gz', sep='\t', compression='gzip')
variant_gene = variant_gene.rename(columns={'snpId':'snp_id', 'geneId':'gene_id', 'geneSymbol': 'gene_name', 'sourceId':'source_id'})

# %% From ClinVar and COSMIC parsed data:
# All mutations table
mutation = pd.read_csv('../db_tables/mutation.tsv', sep= '\t') # es clinvar + cosmic + uniprot
mutation_has_source = pd.read_csv('../db_tables/mutation_has_source.tsv', sep='\t')
source = pd.read_csv('../db_tables/source.tsv', sep='\t')

# Create an auxiliar table
mutation_ids = mutation.merge(mutation_has_source).drop(
    columns=['chromosome','start_genomic', 'end_genomic', 'start_aa', 'end_aa', 'notation_cds', 'notation_aa', 'id_consequence', 'id_insource']
)
# %%
# Our database protein table
protein = pd.read_csv('../db_tables/protein.tsv', sep='\t')
protein_id = protein[['id_protein',	'uniprot_acc']].copy()

# %% Preprocessing

print(f'unique snps in vda dataset: {len(vda.snp_id.unique())}')
vda.snp_id.isnull().any() # all entries have a snp
vda.snp_id.value_counts()
#vda[vda.snp_id == 'rs3184504']  # cambian los disease. Un snp puede tener varios diseases asociados
# %% Generete diseases table
disgenet_disease = vda[
    ['id_disgenet', 'snp_id', 'diseaseid', 'diseasename', 'diseasetype', 'diseaseclass', 'diseasesemantictype', 'score', 'ei', 'yearinitial', 'yearfinal', 'dsi', 'dpi', 'nofpmids', 'source']
    ].copy()
# Add an unique identifier for diseases in disgenet
#disgenet_disease['disgenet_disease_id'] = range(1, len(disgenet_disease)+1)

# %% Subset without diseases
vda_gene = vda[['id_disgenet', 'snp_id', 'chromosome', 'position', 'source']].copy()
# vda_gene.drop_duplicates() # 179589 unique
# Check the chromosomes
#vda_gene.chromosome.value_counts() # Ok, no ensembled chromosomes
# %% Cross references data
# Add uniprot accessions
variant_gene = variant_gene.merge(uniprots)
# Drop source_id col (if the entry comes from VEP or dbsnp)
variant_gene = variant_gene.drop(columns= 'source_id').drop_duplicates()
len(variant_gene.snp_id.unique()) # 158086
# %% Add the proteins data
# Subset here only our LLPS set of proteins (add the id_protein)
variant_gene_llps = variant_gene.merge(protein_id).drop(columns= 'uniprot_acc')
len(variant_gene_llps.id_protein.unique()) # 4187
# %% With this, subset vda dataset by the snp
#vda_1 = vda_gene.merge(variant_gene_llps, on= 'snp_id', how= 'left').drop_duplicates()
vda_1 = vda_gene.merge(variant_gene_llps, on= 'snp_id').drop_duplicates()
# %%
len(vda_1.snp_id.unique()) # 53402
len(vda.snp_id.unique()) # 168051
len(vda_1.id_protein.unique()) # 3818

# %%
# Get SNPs that already exist in the others db (clinvar & cosmic)
snps = ((mutation_ids[mutation_ids.snp_id.notnull()]).snp_id.unique()) # array of unique snps; 165627

# %% SNPs unicos en disgenet: aquellos que no esten en el array snps
snps_unique_disgenet = variant_gene_llps[~variant_gene_llps.snp_id.isin(snps)]
len(snps_unique_disgenet.snp_id.unique()) # 34966
# Dos cosas: 
# 1. agregar estas mutaciones nuevas provenientes de disgenet (34966)
# 2. agregar el source de disgenet (2) y el id_insource a las que ya estan en mutation
# %% 2
# SNPs de disgenet que ya estan en mutation table: agregar el source disgenet (2)
#mutation[mutation.snp_id.isin(snps_unique_disgenet.snp_id.unique())] # vacio, esta ok
mutation[mutation.snp_id.isin(vda_1.snp_id.unique())]
# ok, de aqui me tengo que traer el id de la mutation
snps_also_in_disgenet = mutation[mutation.snp_id.isin(vda_1.snp_id.unique())][['id_mutation', 'snp_id']]
# len(snps_also_in_disgenet.snp_id.unique()) # 25820
# len(vda_1.snp_id.unique()) # 53402

# Data to add in mutation_has_source table
to_add = vda_1.merge(snps_also_in_disgenet, on= 'snp_id')[['id_disgenet', 'id_mutation']]
to_add["id_source"] = 2
#to_add.duplicated().sum() # True
# each row must be unique
to_add.drop_duplicates(inplace= True)
to_add.rename(columns={'id_disgenet': 'id_insource'}, inplace= True)

# %% Concat both tables
mutation_has_source = pd.concat([mutation_has_source, to_add], ignore_index= True).sort_values(by= 'id_mutation')
#mutation_has_source.duplicated().any() # False, ok

######################################################################################################
# %% Ahora agregar las de UNIPROT
# ft_id: unique and stable feature identifier (id_insource)
uniprot_variants = pd.read_csv('../raw_data/uniprot_all_proteins_mutations.tsv.gz', sep='\t')
# %% Generate a separated table for diseases. The id is ft_id col
uniprot_variants_disease = uniprot_variants[
    ['ft_id', 'id_protein', 'snp_id', 'category',
    'mim', 'disease', 'consequence', 'source']
].copy()
# Another table for genes and mutations
uniprot_variants_gene = uniprot_variants[
    ['ft_id', 'id_protein', 'gene_name', 'snp_id', 'change', 'category',
    'from_aa', 'to_aa', 'start_aa', 'end_aa', 'consequence', 'source']
].copy()
# %% Comprobar los SNPs que no estan en mutation.tsv
snps_unique_uniprot = uniprot_variants[~uniprot_variants.snp_id.isin(snps)] # 13027

# SNPs que ya estan en mutations.tsv
# Use dopna() to skip nans in snp_id
snps_also_in_uniprot = mutation[mutation.snp_id.isin(uniprot_variants_gene.snp_id.dropna())][['id_mutation', 'snp_id']] # 8984
# %% With this, add id_mutation to uniprot_variants
uniprot_to_add = uniprot_variants_gene.merge(snps_also_in_uniprot, on= 'snp_id')[['ft_id', 'id_mutation']]
# add source 3 for uniprot
uniprot_to_add["id_source"] = 3
uniprot_to_add.drop_duplicates(inplace= True)
uniprot_to_add.rename(columns={'ft_id': 'id_insource'}, inplace= True) # esto es mutation_has_source
# %% Agrego esas mutations que tambien tienen source de uniprot
mutation_has_source = pd.concat([mutation_has_source, uniprot_to_add], ignore_index= True).sort_values(by= 'id_mutation')
#mutation_has_source.duplicated().any() # False, ok

# %% Now, add mutations unique in disgenet and uniprot to mutation table
# Para agregar las de disgenet me esta faltando el id_consequence (ver)

# Agrego las de uniprot (son todas missense)
mutation_to_add = pd.DataFrame(columns= mutation.columns)
# %% Add notation_aa_col
snps_unique_uniprot["notation_aa"] = "p." + snps_unique_uniprot.from_aa.astype(str) + snps_unique_uniprot.start_aa.astype(str) + snps_unique_uniprot.to_aa.astype(str)
# %% Format
snps_unique_uniprot.replace({"consequence": {"missense": 1}}, inplace= True)
snps_unique_uniprot.rename(columns= {'consequence': 'id_consequence'}, inplace= True)
# %% Add mutations
mutation_to_add = pd.concat([mutation_to_add, snps_unique_uniprot[
    ['snp_id', 'start_aa', 'end_aa', 'notation_aa', 'id_protein', 'id_consequence']
]])
# Add a unique id for each new mutation
mutation_to_add.id_mutation = range(len(mutation)+1, len(mutation)+len(mutation_to_add)+1)
mutation = pd.concat([mutation, mutation_to_add], ignore_index= True)
#mutation.duplicated().any() # False, Ok

# %% Ahora agregar en mutation_has_source
# Traer el ft_id (id_insource)
mutation_has_source_to_add = snps_unique_uniprot.merge(mutation_to_add)
mutation_has_source_to_add = snps_unique_uniprot.merge(mutation_to_add)[['id_mutation', 'ft_id', 'source']].drop_duplicates()
mutation_has_source_to_add.replace({'source': {"uniprot": 3}}, inplace= True)
mutation_has_source_to_add.rename(columns={'ft_id': 'id_insource', 'source': 'id_source'}, inplace= True)
# %% Add to mutation_has_source table
mutation_has_source = pd.concat([mutation_has_source, mutation_has_source_to_add], ignore_index= True).sort_values(by= 'id_mutation')

# Hasta aca actualice: mutation.tsv; mutation_has_source.tsv
# %% Ok, faltaria agregar las unique de disgenet...
# %% save
mutation.to_csv('../db_tables/mutation_new.tsv.gz', sep= '\t', index= False, compression='gzip')
mutation_has_source.to_csv('../db_tables/mutation_has_source_new.tsv', sep= '\t', index= False)
