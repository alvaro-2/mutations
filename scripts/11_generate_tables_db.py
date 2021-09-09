import pandas as pd
import numpy as np
import pyranges as pr

def t_consequence_and_mutation():   
    #comes from parse_clinvar.py
    mutation = pd.read_csv('../raw_data/mutations_clinvar.tsv.gz', sep='\t', compression='gzip') 
    cols_mut = ['id_protein', 'id_mutation', 'snp_id', 'chromosome', 'start_genomic', 'end_genomic', 'start_aa', 'end_aa', 'notation_aa', 'consequence', 'notation_cds']
    mutation = mutation[cols_mut]
    
    mut2 = pd.read_csv('../raw_data/mutations_cosmic.tsv.gz', sep='\t', compression='gzip') 
    mut2 = mut2[cols_mut]
    
    mutation = mutation.append(mut2, ignore_index = True)
    mutation = mutation.drop_duplicates()
    
    # consequence table
    cf = mutation.consequence.value_counts()
    consequence = pd.DataFrame({'id_consequence': range(1, len(cf)+1), 'consequence': cf.index})
    print(f'Generating table consequence.tsv, rows {consequence.shape[0]}')
    consequence.to_csv('../db_tables/consequence.tsv', sep='\t', index = False)

    mutation = mutation.merge(consequence)

    mutation.drop(columns='consequence', inplace= True)

    mutation = mutation[['id_mutation', 'snp_id', 'chromosome', 'start_genomic', 'end_genomic', 'start_aa','end_aa',
                    'notation_cds', 'notation_aa', 'id_protein', 'id_consequence']].sort_values('id_mutation')

    mutation.chromosome = mutation.chromosome.apply(str)
    print(f'Generating table mutation.tsv, rows {mutation.shape[0]}')
    mutation.to_csv('../db_tables/mutation.tsv', sep='\t', index = False)
    return mutation

def t_pfam_proteinpfam_mutationpfam(id_protein, aux_py):
    #pfam domain table
    pfam = pd.read_csv('../raw_data/tablas_disphase_30-08/pfam.csv').rename(columns={'uniprot': 'uniprot_acc', 'tipo': 'pfam_name'})
    
    # Array with unique pfam domains
    pfam_domain = pfam[['pfam_name', 'pfam_acc']].drop_duplicates()
    pfam_domain.rename(columns={'pfam_acc': 'id_pfam', 'pfam_name': 'pfam_domain'}, inplace= True)
    print(f'Generating table pfam_domain.tsv, rows {pfam_domain.shape[0]}')
    pfam_domain.to_csv('../db_tables/pfam_domain.tsv', sep='\t', index= False)

    # ## protein_has_pfam_domain  
    # cols: id_protein, id_pfam, start, end, length
    protein_has_pfam_domain = pfam.merge(id_protein) # agregar col id_protein
    protein_has_pfam_domain['length'] = protein_has_pfam_domain.end - protein_has_pfam_domain.start + 1 # col length
    protein_has_pfam_domain.drop(columns='pfam_name', inplace= True)
    protein_has_pfam_domain = protein_has_pfam_domain.merge(pfam) # to add the col pfam_id
    protein_has_pfam_domain.rename(columns={'pfam_acc': 'id_pfam'}, inplace= True)
    protein_has_pfam_domain = protein_has_pfam_domain[['id_protein', 'id_pfam', 'start', 'end', 'length']].sort_values('id_protein')
    
    #protein_has_pfam_domain.duplicated().any()
    print(f'Generating table protein_has_pfam_domain.tsv, rows {protein_has_pfam_domain.shape[0]}')
    protein_has_pfam_domain.to_csv('../db_tables/protein_has_pfam_domain.tsv', sep='\t', index= False)
    
    # mutation_has_pfam_domain 
    df = pfam.rename(columns={'pfam_name': 'pfam_domain'}).merge(pfam_domain)
    df = df.merge(id_protein) # mapping uniprot_acc - id_protein
    df.drop(columns='uniprot_acc', inplace= True)
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

def t_lowcomplexity_mutationlc(id_protein, aux_lc_py):
    #low complexity
    low_complexity = pd.read_csv('../raw_data/tablas_disphase_30-08/low_complexity_regions.csv').rename(columns={'uniprot': 'uniprot_acc'})
    low_complexity['id_lc'] = range(1, len(low_complexity)+1)
    # Add length col 
    low_complexity['length'] = low_complexity.end - low_complexity.start + 1 
    # Add id_proteins
    low_complexity = low_complexity.merge(id_protein)
    low_complexity.drop(columns='uniprot_acc', inplace= True)
    # Save
    print(f'Generating table low_complexity.tsv, rows {low_complexity.shape[0]}')
    low_complexity.to_csv('../db_tables/low_complexity.tsv', sep='\t', index= False)
    
    # ## mutation_has_low_complexity  
    # cols: id_mutation, id_lc
    
    # Table for LC data
    lc_has = low_complexity.copy()
    lc_has.rename(columns={'id_protein': 'Chromosome', 'start': 'Start', 'end': 'End'}, inplace= True)
        
    # Create the Pyranges objects
    lc_has_py = pr.PyRanges(lc_has)
    
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
    print(f'Generating table mutation_has_low_complexity.tsv, rows {mutation_has_low_complexity.shape[0]}')
    mutation_has_low_complexity.to_csv('../db_tables/mutation_has_low_complexity.tsv', sep='\t', index= False)

def t_disorder(aux_idr_py):
    # # Disorder Tables
    disorder = pd.read_csv('../raw_data/tablas_disphase_30-08/disordered_regions.csv').rename(columns={'uniprot': 'uniprot_acc'})
    disorder['id_idr'] = range(1, len(disorder)+1)
    # ## disorder_region  
    # cols: id_idr, start, end, length, id_protein
    
    # Add length col 
    disorder['length'] = disorder.end - disorder.start + 1 
    disorder_region = disorder.merge(id_protein).sort_values('id_protein')
    disorder_region.drop(columns='uniprot_acc', inplace= True)
    # Save
    print(f'Generating table disorder_region.tsv, rows {disorder_region.shape[0]}')
    disorder_region.to_csv('../db_tables/disorder_region.tsv', sep='\t', index= False)
    
    # ## mutation_has_disorder_region  
    # cols: id_mutation, id_idr
    
    # Table for IDRs data
    idr_has = disorder_region.copy()
    idr_has.rename(columns={'id_protein': 'Chromosome', 'start': 'Start', 'end': 'End'}, inplace= True)
    
    # Create the Pyranges objects
    idr_has_py = pr.PyRanges(idr_has)
    
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
    print(f'Generating table mutation_has_disorder_region.tsv, rows {mutation_has_disorder_region.shape[0]}')
    mutation_has_disorder_region.to_csv('../db_tables/mutation_has_disorder_region.tsv', sep='\t', index= False)

def t_mlo_db_rol_proteinmlo(id_protein):
    # # Rol table 
    mlo_db_rol = pd.read_csv('../raw_data/mlo_db_rol_cleaned.tsv', sep="\t")
    # cols: id_rol, rol
    
    #database_entrada.rol.unique()
    #database_entrada.rol.value_counts()
    
    rol = pd.DataFrame({'rol': mlo_db_rol.rol.value_counts().index, 'id_rol': range(1, len(mlo_db_rol.rol.value_counts())+1)})
    # Save
    print(f'Generating table rol.tsv, rows {rol.shape[0]}')
    rol.to_csv('../db_tables/rol.tsv', sep='\t', index= False)
    
    # # dataset table  
    # cols: id_dataset, dataset
    #database_entrada.db.value_counts()
    
    dataset = pd.DataFrame({'dataset': mlo_db_rol.db.value_counts().index, 'id_dataset': range(1, len(mlo_db_rol.db.value_counts())+1)})
    # Save
    print(f'Generating table dataset.tsv, rows {dataset.shape[0]}')
    dataset.to_csv('../db_tables/dataset.tsv', sep='\t', index= False)
    
    # # MLOs tables  
    # cols: id_mlo, mlo
    #database_entrada.info()
    #database_entrada.mlo.value_counts()
    #len(database_entrada.uniprot_acc.unique()) # OK
    
    mlo = pd.DataFrame({'mlo': mlo_db_rol.mlo.value_counts().index, 'id_mlo': range(1, len(mlo_db_rol.mlo[mlo_db_rol.mlo.notnull()].unique())+1)})
    # Save
    print(f'Generating table mlo.tsv, rows {mlo.shape[0]}')
    mlo.to_csv('../db_tables/mlo.tsv', sep='\t', index= False)
    
    # ## protein_has_mlo  
    # cols: id_protein, id_mlo, id_rol, id_database
    
    #len(database_entrada.uniprot_acc.unique())
    protein_has_mlo = mlo_db_rol.copy()
    
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
    print(f'Generating table protein_has_mlo.tsv, rows {protein_has_mlo.shape[0]}')
    protein_has_mlo.to_csv('../db_tables/protein_has_mlo.tsv', sep='\t', index= False)

def t_source_mutationsource(mutations_source_clinvar):
    # # source
    source = pd.DataFrame({'id_source': [1, 2, 3], 'source': ['clinvar', 'disgenet', 'uniprot', 'cosmic']})
    # Save
    print(f'Generating table source.tsv, rows {source.shape[0]}')
    source.to_csv('../db_tables/source.tsv', sep='\t', index = False)
    
    mutation_has_sources = pd.read_csv('../raw_data/sources_cosmic.tsv.gz', sep='\t', compression='gzip')
    mutation_has_sources = mutations_source_clinvar.append(mutation_has_sources, ignore_index=True)
    # # mutation_has_source
    mutation_has_source.drop_duplicates(inplace= True)
    mutation_has_source = mutation_has_source.merge(source).drop(columns='source')
    # Save
    print(f'Generating table mutation_has_source.tsv, rows {mutation_has_source.shape[0]}')
    mutation_has_source.to_csv('../db_tables/mutation_has_source.tsv', sep='\t', index= False)


def t_citationsource_mutationcitation(mutations_source_clinvar):
    # # citation_source
    mutation_has_citation = pd.read_csv('../raw_data/var_citations.txt', sep='\t')
    mutation_has_citation.columns = mutation_has_citation.columns.str.lower().str.replace(' ',"_").str.replace("-",'_').str.replace('/','_')
    mutation_has_citation = mutation_has_citation.rename(columns={'variationid': 'id_insource'})
    mutation_has_citation = mutation_has_citation[['id_insource', 'citation_source', 'citation_id']]
    
    mutation_has_citation = mutation_has_citation.merge(mutations_source_clinvar).drop(columns=['id_insource']).drop_duplicates()
    
    pmid = pd.read_csv('../raw_data/pubmeds_cosmic.tsv.gz', sep='\t')
    
    mutation_has_citation = mutation_has_citation.append(pmid, ignore_index=True)
    
    #pmid.variationid.isnull().any() # False. This is the ClinVar ID
    mutation_has_citation.drop_duplicates(inplace= True)
    
    citation_source = pd.DataFrame({'name': mutation_has_citation.citation_source.value_counts().index, 'id_citation_source': range(1, len(mutation_has_citation.citation_source.unique())+1) })
    # Save
    print(f'Generating table citation_source.tsv, rows {citation_source.shape[0]}')
    citation_source.to_csv('../db_tables/citation_source.tsv', sep='\t', index= False)
    
    # # mutation_has_citation
    mutation_has_citation = mutation_has_citation.rename(columns={'citation_id': 'id_citation'})
    mutation_has_citation = mutation_has_citation.merge(citation_source.rename(columns={'name': 'citation_source'})).drop(columns= 'citation_source')
    
    # Save
    print(f'Generating table mutation_has_citation.tsv, rows {mutation_has_citation.shape[0]}')
    mutation_has_citation.to_csv('../db_tables/mutation_has_citation.tsv', sep='\t', index= False) 

if __name__ == "__main__":  

    protein = pd.read_csv('../raw_data/uniprot_all_data_associated_to_mlo_preproc.tsv', sep="\t")
    protein = protein[['id_protein', 'uniprot_acc', 'uniprot_other_accesions']]
    
    id_protein = protein[['id_protein', 'uniprot_other_accesions']].rename(columns={'uniprot_other_accesions': 'uniprot_acc'})
    id_protein['uniprot_acc'] = id_protein['uniprot_acc'].str.split(';')
    id_protein = id_protein.explode('uniprot_acc')
    id_protein['uniprot_acc'] = id_protein['uniprot_acc'].str.strip()
    id_protein = id_protein[~id_protein['uniprot_acc'].isnull()]
    id_protein = id_protein[id_protein['uniprot_acc'] != '']
    
    id_protein = pd.concat([protein[['id_protein', 'uniprot_acc']], id_protein])
      
    #consequence and mutation table
    aux_range = t_consequence_and_mutation()
    
    aux_range = aux_range[['start_aa', 'end_aa', 'id_mutation', 'id_protein']]
    aux_range.rename(columns={'id_protein': 'Chromosome', 'start_aa': 'Start', 'end_aa': 'End'}, inplace= True)
    # Pyranges object of mutations
    aux_range = pr.PyRanges(aux_range)
    
    #tables related to pfam
    t_pfam_proteinpfam_mutationpfam(id_protein, aux_range)
    #tables related to low_complexity
    t_lowcomplexity_mutationlc(id_protein, aux_range)
    #tables related to disorder
    t_disorder(aux_range)
    
    2- agregar metodo para regiones llps en proteina y mutacion
    3- agregar metodo para ptms en proteina
    
    #table related to mlos    
    4- add the field reviewed as varchar 2 LT or HT to 
    t_mlo_db_rol_proteinmlo(id_protein)
    
    mutations_source_clinvar = pd.read_csv('../raw_data/mutations_source_clinvar.tsv.gz', sep='\t', compression='gzip')
    
    #table related to source    
    t_source_mutationsource(mutations_source_clinvar)
    #table related to citation
    t_citationsource_mutationcitation(mutations_source_clinvar)    