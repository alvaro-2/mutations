import pandas as pd
import argparse
import numpy

def parse_args():
    parser = argparse.ArgumentParser(description='Pre-process clinvar file')

    parser.add_argument('--mapping_uniprot',
                        dest='protein',
                        help='Uniprot mapping preproc uniprot_all_data_associated_to_mlo_preproc.tsv')
    
    parser.add_argument('--mobidb_dc',
                        dest='modidb',
                        help='Mobidb disorder content disorder_content.csv')

    opts = parser.parse_args()
    return opts

    parser.parse_args()
    
if __name__ == "__main__":
    
    opts = parse_args()
    #uniprot_acc mappping all data (gene_name, gene_id, hgnc_id, etc) of proteins related to mlo
    protein = pd.read_csv(opts.protein, sep="\t")
    
    # Mobidb disorder content
    mobidb = pd.read_csv(opts.modidb).rename(columns={'uniprot': 'uniprot_acc', 'dc': 'disorder_content'})
    print(f'modibd dc {mobidb.shape[0]}')
    
    #take in count some uniprot can be old accession
    p_new_old = protein[['id_protein', 'uniprot_other_accesions']].rename(columns={'uniprot_other_accesions': 'uniprot_acc'})
    p_new_old['uniprot_acc'] = p_new_old['uniprot_acc'].str.split(';')
    p_new_old = p_new_old.explode('uniprot_acc')
    p_new_old['uniprot_acc'] = p_new_old['uniprot_acc'].str.strip()
    p_new_old = p_new_old[~p_new_old['uniprot_acc'].isnull()]
    p_new_old = p_new_old[p_new_old['uniprot_acc'] != '']
    
    p_new_old = pd.concat([protein[['id_protein', 'uniprot_acc']], p_new_old])
    
    #not human uniprot accessions
    #print(set(mobidb['uniprot_acc'].tolist()).difference(set(p_new_old['uniprot_acc'].tolist())))
    
    mobidb = mobidb.merge(p_new_old)
    mobidb = mobidb[['id_protein', 'disorder_content']].drop_duplicates()
    print(f'mobidb_dc mapped llps proteins, rows {mobidb.shape[0]}, unique id_protein {len(set(mobidb["id_protein"].tolist()))}')
    
    #add disorder content
    protein = protein.merge(mobidb, how= 'left')
    
    protein = protein[['id_protein', 'uniprot_acc', 'uniprot_status', 'length', 'uniprot_name', 'hgnc_id', 'gene_id', 'gene_name', 'disorder_content', 'gene_name_synonyms', 'protein_names', 'sequence']]
    k = "|".join(["ubiquinone", "ATP", "NAD", "NADH", "NADP", "NADPH", "NADP\(\+\)", "NMDA", "ADP\/GDP-forming", "GDP-forming", "ADP-forming", "acylating", "ammonia", "flavin-containing", "isomerizing", "carboxylating", "decarboxylating", "acetyl-transferring"])
    protein['protein_names'] = protein['protein_names'].str.replace(' ?\[('+k+')\]', '', regex=True)
    
    print(f'total of proteins {protein.shape[0]}')
    
    protein.to_csv('../db_tables/protein.tsv', sep='\t', index= False)
    
    '''
    print(numpy.nanmax(protein['uniprot_status'].str.len())) #10
    print(numpy.nanmax(protein['sequence'].str.len())) #34350
    print(numpy.nanmax(protein['hgnc_id'].str.len())) #71
    print(numpy.nanmax(protein['gene_name'].str.len())) #75
    print(numpy.nanmax(protein['gene_id'].str.len())) #73
    print(numpy.nanmax(protein['gene_name_synonyms'].str.len())) #267
    print(numpy.nanmax(protein['protein_names'].str.len())) #631
    print(numpy.nanmax(protein['uniprot_name'].str.len())) #16
    print(numpy.nanmax(protein['uniprot_acc'].str.len())) #10
    '''