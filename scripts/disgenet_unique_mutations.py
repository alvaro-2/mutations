import pandas as pd
import os
os.chdir('C:\\Users\\User\\Documents\\mutations\\scripts')

# Este script es para identificar las consequences en las mutaciones unicas en disgenet

# %% this dataset comes from disgenet downloads (https://www.disgenet.org/downloads; Octubre 2020)
vda = pd.read_csv('../raw_data/curated_variant_disease_associations.tsv', sep='\t')
vda.columns = vda.columns.str.lower().str.replace(' ',"_").str.replace("-",'_').str.replace('/','_')
vda = vda.rename(columns={'snpid':'snp_id'})
# %% Mappings UniProts: uniprots id with entrez gene id (https://www.disgenet.org/downloads, Mappings)
uniprots = pd.read_csv('../raw_data/mapa_geneid_4_uniprot_crossref.tsv.gz', sep='\t', compression='gzip').rename(columns={'UniProtKB':'uniprot_acc', 'GENEID':'gene_id'})

# The file contains the mappings of DisGeNET variants (dbSNP Identifiers)
# to NCBI Entrez identifiers according to dbSNP database (https://www.disgenet.org/downloads, Variant-Gene Mappings File)
variant_gene = pd.read_csv('../raw_data/variant_to_gene_mappings.tsv.gz', sep='\t', compression='gzip')
variant_gene = variant_gene.rename(columns={'snpId':'snp_id', 'geneId':'gene_id', 'geneSymbol': 'gene_name', 'sourceId':'source_id'})

# %% Clinvar mutations table
mutation_clinvar = pd.read_csv('../db_tables/mutation.tsv', sep= '\t')
# %% Our database protein table
protein = pd.read_csv('../db_tables/protein.tsv', sep='\t')
protein_id = protein[['id_protein',	'uniprot_acc']].copy()
# %% Dataset con consultas de VEP. one row per snp
cols1 = ['snp_id', 'allele_string', 'start_genome', 'end_genome', 'chromosome', 'assembly', 'most_severe_consequence', 'transcript_consequences']
rsid_data = pd.read_csv('../raw_data/rsid_data.txt', sep= '\t', names= cols1, skiprows= 1)
# %% In this dataset one row per allele, i.e., for C/A/G/T will be three rows for that snp
cols2 = ['snp_id', 'allele_string', 'type', 'ensembl_gene', 'allele_alt', 'from_to_aa', 'cdna_start', 'cdna_end', 'codons', 'impact', 'gene_name', 'cds_start', 'cds_end', 'aa_start', 'aa_end', 'consequence']
allele_data = pd.read_csv('../raw_data/allele_data.txt', sep= '\t', names= cols2, skiprows= 1)

# %% Preprocessing
rsid_data.head()
# drop that column
rsid_data.drop(columns= 'transcript_consequences', inplace= True)
rsid_data.snp_id.isnull().any() # all entries have a snp
((rsid_data.end_genome - rsid_data.start_genome + 1) > 1).any() # not all snps are single position
rsid_data[((rsid_data.end_genome - rsid_data.start_genome + 1) > 1)]

# %% Delete entries whose chromosome be like CHR_HG2334_PATCH
rsid_data.chromosome.value_counts()
# %% Drop those ensembled chromosomes
rsid_data = rsid_data[~rsid_data.chromosome.str.startswith('CHR')]
# %%
rsid_data.chromosome.value_counts()
rsid_data.columns
# %% Generate coordinates table
coordinates = rsid_data[['snp_id', 'chromosome', 'start_genome', 'end_genome']].copy()
# %% allele_data
#  A CoDing Sequence (CDS) is a region of DNA or RNA whose sequence determines the sequence of amino acids in a protein
allele_data.type.value_counts()
print(f'unique snps in vda dataset: {len(vda.snp_id.unique())}')
vda.snp_id.isnull().any() # all entries have a snp
vda.snp_id.value_counts()
# %%
