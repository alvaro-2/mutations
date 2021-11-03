import pandas as pd


pfam = pd.read_csv('../raw_data/pfam_domains.csv').rename(
    columns={'uniprot': 'uniprot_acc', 'pfam_acc': 'id_pfam', 'domain': 'pfam_domain'}
)
# Add this block of code to t_pfam_proteinpfam_mutationpfam() function. Script 11
plds = pd.read_csv('../raw_data/pld.csv')
#plds.duplicated().any() # False, ok
#pfam.duplicated().any() # False, ok

plds.rename(
    columns= {'uniprot': 'uniprot_acc'}, inplace= True
)

# Create a generic id for each PLD region
id_pfam = [ "PLD" + str(i) for i in range(1, len(plds)+1) ]

# Add id_pfam and pfam_domain cols
plds["id_pfam"] = id_pfam
plds["pfam_domain"] = "PLD"

plds.drop(
    columns=['score', 'seq'], inplace= True
)

# Concat both tables
pfam = pd.concat([pfam, plds]).drop_duplicates()

# Generate pfam_domain table
pfam_domain = pfam[['id_pfam', 'pfam_domain']].drop_duplicates() # 3717 Ok. (3581 pfam + 136 pld)
print(f'Generating table pfam_domain.tsv, rows {pfam_domain.shape[0]}')

