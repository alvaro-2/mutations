import pandas as pd
import numpy as np


# ############################## LOAD DATASETS: ##############################  
# ## Disgenet data:  
# * variant - disease association curated; Summary of Curated VDAs   
# ### Cross references  
# * mapa_geneid_4_uniprot_crossref.tsv.gz
# 
# ## From ClinVar and COSMIC parsed data:  
# * mutation.tsv table 
# * mutation_has_source.tsv table  
# * source.tsv table 
# ## From our database:  
# * protein.tsv table

# %%
# this dataset comes from disgenet downloads (https://www.disgenet.org/downloads; Octubre 2020)
vda = pd.read_csv('../raw_data/curated_variant_disease_associations.tsv', sep='\t')
vda.columns = vda.columns.str.lower().str.replace(' ',"_").str.replace("-",'_').str.replace('/','_')
vda = vda.rename(columns={'snpid':'snp_id'})
# Add a generic id for each record
vda['id_disgenet'] = range(1, len(vda)+1)
# %%
# Mappings UniProts: uniprots id with entrez gene id (https://www.disgenet.org/downloads, Mappings)
uniprots = pd.read_csv('../raw_data/mapa_geneid_4_uniprot_crossref.tsv.gz', sep='\t', compression='gzip').rename(columns={'UniProtKB':'uniprot_acc', 'GENEID':'gene_id'})

# The file contains the mappings of DisGeNET variants (dbSNP Identifiers)
# to NCBI Entrez identifiers according to dbSNP database (https://www.disgenet.org/downloads, Variant-Gene Mappings File)
variant_gene = pd.read_csv('../raw_data/variant_to_gene_mappings.tsv.gz', sep='\t', compression='gzip')
variant_gene = variant_gene.rename(columns={'snpId':'snp_id', 'geneId':'gene_id', 'geneSymbol': 'gene_name', 'sourceId':'source_id'})

# %%
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

# # Preprocessing
# ## rsid_data
# %%
rsid_data.head()
# drop that column
rsid_data.drop(columns= 'transcript_consequences', inplace= True)
rsid_data.snp_id.isnull().any() # all entries have a snp
((rsid_data.end_genome - rsid_data.start_genome + 1) > 1).any() # not all snps are single position
rsid_data[((rsid_data.end_genome - rsid_data.start_genome + 1) > 1)]

# %%
# Delete entries whose chromosome be like CHR_HG2334_PATCH
rsid_data.chromosome.value_counts()
# Drop those ensembled chromosomes
rsid_data = rsid_data[~rsid_data.chromosome.str.startswith('CHR')]

# %%
# Generate coordinates table
coordinates = rsid_data[['snp_id', 'chromosome', 'start_genome', 'end_genome']].copy()

# %% allele_data
#  A CoDing Sequence (CDS) is a region of DNA or RNA whose sequence determines the sequence of amino acids in a protein
allele_data.type.value_counts()
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
#vda_gene.head()
# vda_gene.drop_duplicates() # 179589 unique
# Check the chromosomes
#vda_gene.chromosome.value_counts() # Ok, no ensembled chromosomes
# %% Cross references data
# Add uniprot accessions
variant_gene = variant_gene.merge(uniprots)
# Drop source_id col (if the entry comes from VEP or dbsnp)
variant_gene = variant_gene.drop(columns= 'source_id').drop_duplicates()
variant_gene.head()
len(variant_gene.snp_id.unique())
# %%
variant_gene.snp_id.value_counts()
#variant_gene[variant_gene.snp_id == 'rs17119346']
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
# From allele_data: only get snps from LLPS proteins
allele_data_1 = allele_data[allele_data.snp_id.isin(vda_1.snp_id.unique())] 
#  A CoDing Sequence (CDS) is a region of DNA or RNA whose sequence determines the sequence of amino acids in a protein
allele_data_1 = allele_data_1[['snp_id', 'allele_string','allele_alt', 'type', 'gene_name',	'codons', 'cds_start', 'cds_end',  'from_to_aa', 'impact', 'aa_start', 'aa_end', 'consequence']]
allele_data_1 # 10727 rows
# %%
# Just keep those entries with changes in from_to_aa col
allele_data_1 = allele_data_1[allele_data_1.from_to_aa != "-"]
# Add genomic coordinates
allele_data_1 = allele_data_1.merge(coordinates, on= 'snp_id')
allele_data_1 # 5397 rows
# %%
len(allele_data.snp_id.unique()) # 5665
len(allele_data_1.snp_id.unique()) #3580; this is the subset with only llps proteins
# %%
# Get SNPs that already exist in the others db (clinvar & cosmic)
snps = ((mutation_ids[mutation_ids.snp_id.notnull()]).snp_id.unique()) # array of unique snps; 165627
# %%
# VER
#vda_1[vda_1.snp_id.isin(snp_clinvar)] # 27753 rows
# %%
# If the snp already exist in mutation table, pass it
subset_mutation_disgenet = allele_data_1[allele_data_1.snp_id.isin(snps)].copy() # 5184
#subset_clinvar_disgenet["source"] = "clinvar;disgenet"
subset_mutation_disgenet
# %%
subset_only_disgenet = allele_data_1[~allele_data_1.snp_id.isin(snps)].copy()
subset_only_disgenet["source"] = "disgenet"
subset_only_disgenet

#%%
#subset_clinvar_disgenet
# %%
#mutation[['id_mutation', 'snp_id']]



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

# %%
len(snps_also_in_disgenet.snp_id.unique())
len(vda_1.snp_id.unique())
# Data to add in mutation_has_source table
to_add = vda_1.merge(snps_also_in_disgenet, on= 'snp_id')[['id_disgenet', 'id_mutation']]
to_add["id_source"] = 2
#to_add.duplicated().sum() # True
# each row must be unique
to_add.drop_duplicates(inplace= True)
to_add.rename(columns={'id_disgenet': 'id_insource'}, inplace= True)

# %% Concat both tables
mutation_has_source = pd.concat([mutation_has_source, to_add], ignore_index= True).sort_values(by= 'id_mutation')
# %%
mutation_has_source.duplicated().any() # False, ok
# %% Subseteo mutation_has_source por el id_mutation de los snps en disgenet
# to_update = mutation_has_source[mutation_has_source.id_mutation.isin(snps_also_in_disgenet.id_mutation)]
# to_update.duplicated().any() # False, ok

# %% Agrego el id de source de disgenet (2).
# xx = to_update.append(
#     to_update.assign(**{'id_source': 2}), ignore_index= True
#     ).sort_values(by= 'id_mutation')

# not finished yet