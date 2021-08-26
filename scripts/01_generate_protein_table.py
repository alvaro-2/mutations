import pandas as pd
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description='Pre-process clinvar file')

    parser.add_argument('--mapping_uniprot',
                        dest='protein',
                        help='Uniprot mapping preproc uniprot_all_data_associated_to_mlo_preproc.tsv')
    
    parser.add_argument('--mobidb_dc',
                        dest='modidb',
                        help='Mobidb disorder content dc_mobidb_lite.csv')

    opts = parser.parse_args()
    return opts

    parser.parse_args()
    
if __name__ == "__main__":
    
    opts = parse_args()
    #uniprot_acc mappping all data (gene_name, gene_id, hgnc_id, etc) of proteins related to mlo
    protein = pd.read_csv(opts.protein, sep="\t")
    
    # Mobidb disorder content
    mobidb = pd.read_csv(opts.modidb).rename(columns={'uniprot': 'uniprot_acc', 'dc': 'disorder_content'})
    
    #take in count the uniprot O43236 is the current accesion for Q8NEP4 (old)
    mobidb['uniprot_acc'] = mobidb['uniprot_acc'].str.replace("Q8NEP4", "O43236")
    
    #add disorder content
    protein = protein.merge(mobidb, how= 'left')
    
    protein = protein[['id_protein', 'uniprot_acc', 'uniprot_status', 'length', 'uniprot_name', 'hgnc_id', 'gene_id', 'gene_name', 'disorder_content', 'gene_name_synonyms', 'protein_names', 'sequence']]
    k = "|".join(["ubiquinone", "ATP", "NAD", "NADH", "NADP", "NADPH", "NADP\(\+\)", "NMDA", "ADP\/GDP-forming", "GDP-forming", "ADP-forming", "acylating", "ammonia", "flavin-containing", "isomerizing", "carboxylating", "decarboxylating", "acetyl-transferring"])
    protein['protein_names'] = protein['protein_names'].str.replace(' ?\[('+k+')\]', '', regex=True)
    
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