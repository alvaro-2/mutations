import pandas as pd
import argparse
import re
import os

# python 03_parse_clinvar.py --clinvar ../raw_data/vs.tsv.gz --refseq ../raw_data/uniprot-refseq-nucleotido.tsv --prot ../db_tables/protein.tsv --out ../raw_data

def parse_args():
    parser = argparse.ArgumentParser(description='Parse clinvar mutations and merge with target proteins.')

    parser.add_argument('--clinvar',
                        dest='clinvar_file',
                        help='pre-processed clinvar variants (vs.tsv.gz)')

    parser.add_argument('--mapping_uniprot', 
                        dest='protein_uniprot',
                        help='File with uniprot mapping (uniprot_all_data_associated_to_mlo_preproc.tsv)')
                        
    parser.add_argument('--out', 
                        dest='outfolder',
                        help='Output files: diseases_clinvar.tsv.gz, mutations_with_source_clinvar.tsv.gz, mutations_clinvar.tsv.gz')

    opts = parser.parse_args()
    return opts

def mapp_by_all_ids(arg_files):
    
    vs_refseq = pd.read_csv(arg_files.clinvar_file, sep='\t', compression='gzip').rename(columns={'start': 'start_genomic', 'stop': 'end_genomic'})
    vs_refseq['index_line'] = list(range(0, len(vs_refseq)))
    #merge at the end to put the data
    vs_ids_other = vs_refseq[['index_line', 'refseq', 'snp_id', 'variationid', 'chromosome', 'start_genomic', 'end_genomic', 'type', 'notation_aa', 'notation_cds', 'origin', 'phenotypeids', 'phenotypelist', 'otherids']]
    #only kept in vs_refseq the 'index_line' and the ids used to mapp
    #by the snp, some refseq can be associated to two or more chromosomes for example to X and Y 
    vs_refseq = vs_refseq[['index_line', 'refseq', 'chromosome']]
    
    print(f'clinvar {vs_refseq.shape[0]}')     
    coln = ['index_line', 'id_protein']
        
    #merge by refseq
    protein_uniprot = pd.read_csv(arg_files.protein_uniprot,  sep="\t").rename(columns={'refseq_canonic': 'refseq'})
    protein_uniprot = protein_uniprot[['id_protein', 'refseq', 'chromosome']].copy()

    #split the field refseq_canonic
    #NP_XXX.x;NM_XXX.x
    #NP_XXX.x;NM_XXX.x.;NP_YYY.y;NM_YYY.y
    #split by ".;"
    protein_uniprot['refseq'] = protein_uniprot['refseq'].str.split('\.;')
    protein_uniprot = protein_uniprot.explode('refseq')
    protein_uniprot['refseq'] = protein_uniprot['refseq'].str.strip()
    protein_uniprot = protein_uniprot[~protein_uniprot['refseq'].isnull()]
    protein_uniprot = protein_uniprot[protein_uniprot['refseq'] != '']

    #split by ";"
    protein_uniprot['refseq'] = protein_uniprot['refseq'].str.split(';')
    protein_uniprot = protein_uniprot.explode('refseq')
    protein_uniprot['refseq'] = protein_uniprot['refseq'].str.strip()
    protein_uniprot = protein_uniprot[~protein_uniprot['refseq'].isnull()]
    protein_uniprot = protein_uniprot[protein_uniprot['refseq'] != '']
    
    #split by "." to delete the version
    protein_uniprot['refseq'] = protein_uniprot['refseq'].map(lambda x: x.split("."))
    protein_uniprot['refseq'] = protein_uniprot['refseq'].str[0]
    
    protein_uniprot_nochr = protein_uniprot[(protein_uniprot['chromosome'].isnull()) | (protein_uniprot['chromosome'] == '')]
    protein_uniprot_nochr = protein_uniprot_nochr[['id_protein', 'refseq']]
    
    #split chromosome
    #split by ";"
    protein_uniprot['chromosome'] = protein_uniprot['chromosome'].str.split(';')
    protein_uniprot = protein_uniprot.explode('chromosome')
    protein_uniprot['chromosome'] = protein_uniprot['chromosome'].str.strip()
    protein_uniprot = protein_uniprot[~protein_uniprot['chromosome'].isnull()]
    protein_uniprot = protein_uniprot[protein_uniprot['chromosome'] != '']
    
    #mapping by refseq and chromosome if exist
    vs_refseq_chr = vs_refseq.merge(protein_uniprot)
    vs_refseq = vs_refseq[~vs_refseq['index_line'].isin(vs_refseq_chr['index_line'].tolist())]
    #mapping refseq for unsertain chromosome
    vs_refseq_nochr = vs_refseq.merge(protein_uniprot_nochr)
    
    vs_refseq = pd.concat([vs_refseq_chr, vs_refseq_nochr])
    print(f'clinvar merge by ref_seq {vs_refseq.shape[0]}, unique {len(vs_refseq["index_line"].unique().tolist())}')
   
    #some id_protein can get more than one refseq
    #kept only the most mutated transcript, ...
    y = vs_refseq[['refseq', 'id_protein']]   
    dd = y.groupby(['id_protein', 'refseq']).size().reset_index(name='counts')
    dd = dd.sort_values(by=['id_protein', 'counts', 'refseq'], ascending=[False, False, True])
    dd = dd.drop_duplicates(subset=['id_protein'], keep='first')
    dd = dd[['id_protein', 'refseq']]
    vs_refseq = vs_refseq.merge(dd)
    print(f'rows mapped by only one refseq by id_protein {vs_refseq.shape[0]}') 
    
    #analizar refseq con disitino id_protein
    y = vs_refseq[['refseq', 'id_protein', 'chromosome']]
    y = y.drop_duplicates()
    dd = y.groupby(['refseq']).size().reset_index(name='counts')
    print('refseq to more than one id_protein should be analized later')
    print('some snps are mapping to two or more chromosomes')
    dd = dd[dd.counts > 1]
    print(dd)
    print(y[y['refseq'].isin(dd['refseq'])])
    
    vs_refseq = vs_refseq[coln]
    print(vs_refseq[vs_refseq['index_line'].isin(vs_refseq[vs_refseq.duplicated(['index_line'])]['index_line'].unique().tolist())])
    
    vs_refseq = vs_refseq.drop_duplicates()
    print(f'clinvar without duplicates {vs_refseq.shape[0]}, unique {len(vs_refseq["index_line"].unique().tolist())}') 
    
    print(f'clinvar proteins {len(vs_refseq["id_protein"].unique())}')
    vs_refseq = vs_refseq.merge(vs_ids_other).drop(columns=['index_line'])
    return vs_refseq

protein_letters_3to1 = {
    'Ala': 'A',
    'Cys': 'C',
    'Asp': 'D',
    'Glu': 'E',
    'Phe': 'F',
    'Gly': 'G',
    'His': 'H',
    'Ile': 'I',
    'Lys': 'K',
    'Leu': 'L',
    'Met': 'M',
    'Asn': 'N',
    'Pro': 'P',
    'Gln': 'Q',
    'Arg': 'R',
    'Ser': 'S',
    'Thr': 'T',
    'Val': 'V',
    'Trp': 'W',
    'Tyr': 'Y',
    'Asx': 'B',
    'Xaa': 'X',
    'Glx': 'Z',
    'Xle': 'J',
    'Sec': 'U',
    'Pyl': 'O',
    'Ter': '*'
}

def separar_en_cols(df, column, conseq, conseq_regex, override=False):
    '''
    recibe un DataFrame, el nombre de una columna auxiliar (column)
    y un string con el tipo de consecuencia (conseq). La col. auxiliar
    es una tupla con los elementos implicados en una mutacion
    como la siguiente (aa1, start_pos, aa2, end_pos, aa/s_nuevos).
    Devuelve el DataFrame df con estas 5 nuevas columnas
    '''
    df_crop = df.copy()
      
    if override:
        df_crop['aux'] = df_crop[column].str.findall(conseq_regex).str[0]
    else:
        df_crop['aux'] = df_crop[column].str.findall('^p\.([\*A-Z])(\d+)_?([\*A-Z])?(\d+)?'+conseq_regex+'(.*)$').str[0]
    
    df_crop = df_crop[~df_crop['aux'].isnull()]
    
    if conseq == "repeted":  #"^p\.(\d+)_?(\d+)?([A-Z]*)\[(\d+)\]$"
        df_crop['start_aa'] = df_crop['aux'].map(lambda x: x[0])
        df_crop['end_aa'] = df_crop['aux'].map(lambda x: x[1] if x[1] != '' else x[0]) 
        
    else:
        # start position
        df_crop['start_aa'] = df_crop['aux'].map(lambda x: x[1])
        # end position
        df_crop['end_aa'] = df_crop['aux'].map(lambda x: x[3] if x[3] != '' else x[1]) # poner en el end el start si no hay end

    # consecuencia de la mutacion
    df_crop['consequence'] = conseq

    df_crop = df_crop.drop(columns=['aux'])
    return df_crop

def put_mutation_consequence(mutations):
    mutations['index_line'] = list(range(0, len(mutations)))
    
    # A DataFrame for each molecular consequence
    synonym = separar_en_cols(mutations, "notation_aa", "synonym", "=$")   
    mutations = mutations[~mutations['index_line'].isin(synonym['index_line'].tolist())]
    print(f"Found {synonym.shape[0]} synonyms")
    
    delins = separar_en_cols(mutations, "notation_aa", "delins", "delins")
    mutations = mutations[~mutations['index_line'].isin(delins['index_line'].tolist())]
    print(f"Found {delins.shape[0]} delins")
    
    deletions = separar_en_cols(mutations, "notation_aa", "deletion", "del") # finish with del
    mutations = mutations[~mutations['index_line'].isin(deletions['index_line'].tolist())]
    print(f"Found {deletions.shape[0]} deletions")
    
    insertions = separar_en_cols(mutations, "notation_aa", "insertion", "ins")
    mutations = mutations[~mutations['index_line'].isin(insertions['index_line'].tolist())]
    print(f"Found {insertions.shape[0]} inserions")
    
    frameshift = separar_en_cols(mutations, "notation_aa", "frameshift", 'fs') #
    mutations = mutations[~mutations['index_line'].isin(frameshift['index_line'].tolist())]
    print(f"Found {frameshift.shape[0]} frameshift")
    
    duplications = separar_en_cols(mutations, "notation_aa", "duplication", "dup")
    mutations = mutations[~mutations['index_line'].isin(duplications['index_line'].tolist())]
    print(f"Found {duplications.shape[0]} duplications")
    
    nostop = separar_en_cols(mutations, "notation_aa", "nostop", "ext") 
    mutations = mutations[~mutations['index_line'].isin(nostop['index_line'].tolist())]
    print(f"Found {nostop.shape[0]} nostop")    
    
    nonsense = separar_en_cols(mutations, "notation_aa", "nonsense", "\*$") # positiv lookbehind search! must have a number before, some delins insert a Ter
    mutations = mutations[~mutations['index_line'].isin(nonsense['index_line'].tolist())]
    print(f"Found {nonsense.shape[0]} nonsense")    
    
    missense = separar_en_cols(mutations, "notation_aa", "missense", '')
    mutations = mutations[~mutations['index_line'].isin(missense['index_line'].tolist())]
    print(f"Found {missense.shape[0]} missense")
    
    repeted = separar_en_cols(mutations, "notation_aa", "repeted", "^p\.(\d+)_?(\d+)?([A-Z]*)\[(\d+)\]$", override=True)   
    mutations = mutations[~mutations['index_line'].isin(repeted['index_line'].tolist())]
    print(f"Found {repeted.shape[0]} repeted")
    
    print(f"No mapped mutations {mutations.shape[0]}")
    print(mutations['notation_aa'].unique().tolist())
    ## Concatenate everything

    tables = [synonym, deletions, delins, duplications, frameshift, insertions, missense, nonsense, nostop, repeted]
    mutations = pd.concat(tables)
    mutations = mutations.drop(columns=['index_line']).drop_duplicates()
    print(f"Total mutations: {mutations.shape[0]}")
   
    return mutations

if __name__ == "__main__":

    opts = parse_args()

    if not os.path.exists(opts.outfolder) or not os.path.isdir(opts.outfolder):
        os.mkdir(opts.outfolder)
    
    mutations = mapp_by_all_ids(opts)    
    print(f'clinvar mutations to mapp AA: {mutations.shape[0]}')
    
    #change the "notation_AA" from 3 letter code to one letter code    
    # Create a regular expression from all of the dictionary keys
    regex = re.compile("|".join(map(re.escape, protein_letters_3to1.keys())))
    # For each match, look up the corresponding value in the dictionary
    mutations['notation_aa'] = [regex.sub(lambda match: protein_letters_3to1[match.group(0)], mut) for mut in mutations['notation_aa'].tolist()]

    mutations = put_mutation_consequence(mutations)
    
    # Assign id_mutation. A mutation is unique by the indexes : 'id_protein', 'notation_aa', 'notation_cds'
    subset = mutations[['id_protein', 'notation_aa', 'notation_cds', 'chromosome', 'start_genomic', 'end_genomic']].drop_duplicates()
    subset['id_mutation'] = range(1, len(subset)+1)
    print(f'Subset (unique mutations) length: {len(subset)}')
    mutations = mutations.merge(subset)
    dd = mutations.groupby(['id_mutation']).size().reset_index(name='counts')
    dd = dd[dd.counts > 1]
    print(dd)
    print(f'All mutations after merge: {mutations.shape[0]}')
    print(mutations.columns.tolist())
    
    
    # Another table for diseases
    diseases = mutations[['id_mutation', 'origin', 'phenotypeids', 'phenotypelist', 'otherids']].copy()
    diseases.to_csv(os.path.join(opts.outfolder, 'diseases_clinvar.tsv.gz'), sep='\t', index= False, compression='gzip')
    
    # Add ClinVar source
    mutations_with_source = mutations[['id_mutation', 'variationid']].drop_duplicates().rename(columns={'variationid': 'id_insource'})
    mutations_with_source['source'] = 'clinvar'
    mutations_with_source.to_csv(os.path.join(opts.outfolder, 'mutations_source_clinvar.tsv.gz'), sep='\t', index= False, compression='gzip')
    
    # drop columns
    mutations.drop(columns=['variationid', 'origin', 'phenotypeids', 'phenotypelist', 'otherids'], inplace= True)
    # one file, one mutation
    mutations.drop_duplicates(ignore_index = True, inplace= True)
    print(f'Length mutations after drop duplicates: {mutations.shape[0]}')
    print(f'clinvar proteins {len(mutations["id_protein"].unique())}')
    mutations.to_csv(os.path.join(opts.outfolder, 'mutations_clinvar.tsv.gz'), sep='\t', index= False, compression='gzip')
