# en este script se generan pfam_domain.tsv, protein_has_pfam_domain.tsv y mutation_has_pfam_domain.tsv
import pandas as pd
import numpy as np
import pyranges as pr

pfam = pd.read_csv('../raw_data/pfam_domains.csv').rename(
    columns={'uniprot': 'uniprot_acc', 'pfam_acc': 'id_pfam', 'domain': 'pfam_domain'}
)

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

# ## protein_has_pfam_domain  
# cols: id_protein, id_pfam, start, end, length
protein_has_pfam_domain = pfam.merge(id_protein) # agregar col id_protein
protein_has_pfam_domain['length'] = protein_has_pfam_domain.end - protein_has_pfam_domain.start + 1 # col length
protein_has_pfam_domain.drop(columns=['pfam_domain', 'uniprot_acc'], inplace= True)
protein_has_pfam_domain = protein_has_pfam_domain.sort_values('id_protein')

#protein_has_pfam_domain.duplicated().any() # False, Ok
print(f'Generating table protein_has_pfam_domain.tsv, rows {protein_has_pfam_domain.shape[0]}')
    
# mutation_has_pfam_domain 
df = protein_has_pfam_domain.copy()
df.rename(columns={'id_protein': 'Chromosome', 'start': 'Start', 'end': 'End'}, inplace= True)

# Create the pyranges object of pfam domains
df_py = pr.PyRanges(df)


# Create a PyRanges object of the mutations
mutation = pd.read_csv('../db_tables/mutation.tsv', sep= '\t')

aux_range = mutation[['start_aa', 'end_aa', 'id_mutation', 'id_protein']].copy()
aux_range.rename(columns={'id_protein': 'Chromosome', 'start_aa': 'Start', 'end_aa': 'End'}, inplace= True)
# Pyranges object of mutations
aux_range = pr.PyRanges(aux_range)



# Join both pyranges object: this assings mutations to pfam domains
pfam_py = df_py.join(aux_range, strandedness= False, slack= 1)  # strandedness= False doesnt take count of the chain strand; slack= 1 include bound
#pfam_py.head() # Start and End are from the pfam domain in that protein (a protein may have the same pfam domain repeated at different positions along its sequence).
                # Start_b and End_b are from the mutation in this case

# Pyranges to DataFrame
mutation_has_pfam_domain = pfam_py.df[['id_mutation', 'Chromosome', 'id_pfam', 'Start', 'End']] # cols to keep
mutation_has_pfam_domain.rename(columns={'Chromosome': 'id_protein', 'Start': 'start', 'End': 'end'}, inplace= True)
print(f'Generating table mutation_has_pfam_domain.tsv, rows {mutation_has_pfam_domain.shape[0]}')
# mutation_has_pfam_domain[mutation_has_pfam_domain.id_pfam.str.startswith("PLD")] # 3647 mutations in PLDs
#mutation_has_pfam_domain.to_csv('../db_tables/mutation_has_pfam_domain.tsv', sep='\t', index= False)
# Hasta aca todo ok


######################################################################################################3
def t_pfam_proteinpfam_mutationpfam(id_protein, aux_py):
    #pfam domain table
    pfam = pd.read_csv('../raw_data/tablas_disphase_30-08/pfam_domains.csv').rename(columns={'uniprot': 'uniprot_acc', 'pfam_acc': 'id_pfam', 'domain': 'pfam_domain'})    
    #columns uniprot pfam_acc start end domain
    
    # PLDs table
    plds = pd.read_csv('../raw_data/pld.csv')
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
      
    # Array with unique pfam domains
    pfam_domain = pfam[['id_pfam', 'pfam_domain']].drop_duplicates()
    print(f'Generating table pfam_domain.tsv, rows {pfam_domain.shape[0]}')
    pfam_domain.to_csv('../db_tables/pfam_domain.tsv', sep='\t', index= False)

    # ## protein_has_pfam_domain  
    # cols: id_protein, id_pfam, start, end, length
    protein_has_pfam_domain = pfam.merge(id_protein) # agregar col id_protein
    protein_has_pfam_domain['length'] = protein_has_pfam_domain.end - protein_has_pfam_domain.start + 1 # col length
    protein_has_pfam_domain.drop(columns=['pfam_domain', 'uniprot_acc'], inplace= True)
    protein_has_pfam_domain = protein_has_pfam_domain.sort_values('id_protein')
    
    #protein_has_pfam_domain.duplicated().any()
    print(f'Generating table protein_has_pfam_domain.tsv, rows {protein_has_pfam_domain.shape[0]}')
    protein_has_pfam_domain.to_csv('../db_tables/protein_has_pfam_domain.tsv', sep='\t', index= False)
    
    # mutation_has_pfam_domain 
    df = protein_has_pfam_domain.copy()
    df.rename(columns={'id_protein': 'Chromosome', 'start': 'Start', 'end': 'End'}, inplace= True)
    
    # Create the pyranges object of pfam domains
    df_py = pr.PyRanges(df)
    
    # Join both pyranges object: this assings mutations to pfam domains
    pfam_py = df_py.join(aux_py, strandedness= False, slack= 1)  # strandedness= False doesnt take count of the chain strand; slack= 1 include bound
    #pfam_py.head() # Start and End are from the pfam domain in that protein (a protein may have the same pfam domain repeated at different positions along its sequence).
                    # Start_b and End_b are from the mutation in this case
    
    # Pyranges to DataFrame
    mutation_has_pfam_domain = pfam_py.df[['id_mutation', 'Chromosome', 'id_pfam', 'Start', 'End']] # cols to keep
    mutation_has_pfam_domain.rename(columns={'Chromosome': 'id_protein', 'Start': 'start', 'End': 'end'}, inplace= True)
    print(f'Generating table mutation_has_pfam_domain.tsv, rows {mutation_has_pfam_domain.shape[0]}')
    mutation_has_pfam_domain.to_csv('../db_tables/mutation_has_pfam_domain.tsv', sep='\t', index= False)