import pandas as pd
import argparse
import re
import os

# python 09_all_mutations.py --inout ../raw_data

def parse_args():
    parser = argparse.ArgumentParser(description='Parse clinvar mutations and merge with target proteins.')

    parser.add_argument('--inout', 
                        dest='inoutfolder',
                        help='Output files: diseases_clinvar.tsv.gz, mutations_with_source_clinvar.tsv.gz, mutations_clinvar.tsv.gz')

    opts = parser.parse_args()
    return opts

if __name__ == "__main__":

    opts = parse_args()
    pd.set_option('display.max_columns', None)
    # mutations clinvar    
    mutations = pd.read_csv(os.path.join(opts.inoutfolder, 'mutations_clinvar.tsv.gz'), sep='\t', dtype="str", compression='gzip')
    mutations = mutations.drop(columns=['refseq', 'type']).drop_duplicates()
    print(f'mutations clinvar {len(set(mutations["id_mutation"].tolist()))}')
    print(f'{mutations[(mutations["notation_cds"] == "c.2573T>G") & (mutations["id_protein"] == "4076")]}')
    print(mutations.info())
    first_id_mut = max(mutations["id_mutation"].astype(int).tolist()) + 1
    #mutation cosmic
    mutations_cosmic_all = pd.read_csv(os.path.join(opts.inoutfolder, 'mutations_parse_cosmic.tsv.gz'), sep='\t', dtype="str", compression='gzip')
    print(f'{mutations_cosmic_all[(mutations_cosmic_all["notation_cds"] == "c.2573T>G") & (mutations_cosmic_all["id_protein"] == "4076")]}')
    
    mutations_cosmic = mutations_cosmic_all.copy()
    mutations_cosmic = mutations_cosmic[['id_protein', 'chromosome', 'start_genomic', 'end_genomic', 'notation_aa', 'notation_cds', 'start_aa', 'end_aa', 'consequence']]
    mutations_cosmic = mutations_cosmic.drop_duplicates()
    print(mutations_cosmic.info())
    print(f'mutations cosmic {mutations_cosmic.shape[0]}')
    print(f'{mutations_cosmic[(mutations_cosmic["notation_cds"] == "c.2573T>G") & (mutations_cosmic["id_protein"] == "4076")]}')
    
    mutations_cosmic = mutations_cosmic.merge(mutations, how='left')
    print(f'{mutations_cosmic[(mutations_cosmic["notation_cds"] == "c.2573T>G") & (mutations_cosmic["id_protein"] == "4076")]}')
    
    cosmic_repeted = mutations_cosmic[~mutations_cosmic['id_mutation'].isnull()]
    mutations_cosmic = mutations_cosmic[mutations_cosmic['id_mutation'].isnull()]
    
    print(f'cosmic and clinvar mutations {cosmic_repeted.shape[0]}')
    print(f'new cosmic mutations {mutations_cosmic.shape[0]}')

    mutations_cosmic['id_mutation'] = list(range(first_id_mut, first_id_mut + len(mutations_cosmic), 1))

    mutations_cosmic = mutations_cosmic.append(cosmic_repeted, ignore_index=True)
    mutations_cosmic.to_csv(os.path.join(opts.inoutfolder, 'mutations_cosmic.tsv.gz'), sep='\t', index= False, compression='gzip')
        
    mutations_cosmic_all = mutations_cosmic_all.merge(mutations_cosmic)    
    
    disease = mutations_cosmic_all[["id_mutation", "Primary_site", "Site_subtype_1", "Site_subtype_2", "Site_subtype_3", "Primary_histology", "Histology_subtype_1", "Histology_subtype_2", "Histology_subtype_3"]]
    disease = disease.drop_duplicates()
    disease.to_csv(os.path.join(opts.inoutfolder, 'diseases_cosmic.tsv.gz'), sep='\t', index= False, compression='gzip')
    
    pubmed = mutations_cosmic_all[["id_mutation", "pubmed"]].rename(columns={"pubmed": 'id_citation'}) 
    pubmed['citation_source'] = 'PubMed'
    pubmed = pubmed.drop_duplicates()
    pubmed = pubmed[pubmed['id_citation'].notnull()]    
    pubmed.to_csv(os.path.join(opts.inoutfolder, 'pubmeds_cosmic.tsv.gz'), sep='\t', index= False, compression='gzip')
    
    sourced = mutations_cosmic_all[["id_mutation", "GENOMIC_MUTATION_ID"]].rename(columns={"GENOMIC_MUTATION_ID": 'id_insource'})
    sourced['source'] = 'cosmic'    
    sourced = sourced.drop_duplicates()
    sourced.to_csv(os.path.join(opts.inoutfolder, 'sources_cosmic.tsv.gz'), sep='\t', index= False, compression='gzip')