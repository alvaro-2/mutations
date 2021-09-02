import argparse
import pandas as pd
import numpy as np
import os

def parse_args():
    parser = argparse.ArgumentParser(description= 'Pre-process UniProt humsavar.txt dataset')

    parser.add_argument('--cosmic_mutations',
                        dest= 'cosmic',
                        help= 'Cosmic mutations "CosmicMutantExport_sel.tsv"')
    
    parser.add_argument('--mapping_uniprot', 
                        dest='protein_uniprot',
                        help='File with uniprot mapping (uniprot_all_data_associated_to_mlo_preproc.tsv)')
                        
    parser.add_argument('--out', 
                        dest='outfolder',
                        help='Output files: mutations_cosmic.tsv.gz, diseases_cosmic.tsv.gz, pubmeds_cosmic.tsv.gz')
    
    opts = parser.parse_args()
    return opts

def mapping_enst(opts):
    #merge by enst
    protein_uniprot = pd.read_csv(opts.protein_uniprot,  sep="\t").rename(columns={'ensembl_canonic': 'enst'})
    protein_uniprot = protein_uniprot[['id_protein', 'enst', 'chromosome']].copy()

    #split the field enst
    #ENSTXXX;ENSPXXX;ENSGXXX
    #ENSTXXX;ENSPXXX;ENSGXXX.;ENSTYYY;ENSPYYY;ENSGYYY
    #split by ".;"
    protein_uniprot['enst'] = protein_uniprot['enst'].str.split('\.;')
    protein_uniprot = protein_uniprot.explode('enst')
    protein_uniprot['enst'] = protein_uniprot['enst'].str.strip()
    protein_uniprot = protein_uniprot[~protein_uniprot['enst'].isnull()]
    protein_uniprot = protein_uniprot[protein_uniprot['enst'] != '']

    #split by ";", kept the first is the ENST
    protein_uniprot['enst'] = protein_uniprot['enst'].map(lambda x: x.split(";"))
    protein_uniprot['enst'] = protein_uniprot['enst'].str[0]
    
    #split chromosome
    #split by ";"
    protein_uniprot['chromosome'] = protein_uniprot['chromosome'].str.split(';')
    protein_uniprot = protein_uniprot.explode('chromosome')
    protein_uniprot['chromosome'] = protein_uniprot['chromosome'].str.strip()
    protein_uniprot = protein_uniprot[~protein_uniprot['chromosome'].isnull()]
    protein_uniprot = protein_uniprot[protein_uniprot['chromosome'] != '']  
    
    
    #cosmic mutations
    mutation = pd.read_csv(opts.cosmic, sep= '\t', dtype='str').rename(columns={'AA': 'notation_aa', 'CDS': 'notation_cds', 'start': 'start_genomic', 'end': 'end_genomic'})
    
    mutation['chromosome'] = mutation['chromosome'].replace(['23', '24', '25'], ['X', 'Y', 'MT'])
    
    mutation = mutation[['enst', 'notation_aa', 'notation_cds', 'Primary_site', 'Site_subtype_1', 'Site_subtype_2', 'Site_subtype_3', 'Primary_histology', 'Histology_subtype_1', 'Histology_subtype_2', 'Histology_subtype_3', 'GENOMIC_MUTATION_ID', 'chromosome', 'start_genomic', 'end_genomic', 'pubmed']]  
    mutation = mutation.drop_duplicates()
    
    other_data = mutation.copy()  
    
    #to avoid the rows duplication
    mutation = mutation[['enst', 'notation_aa', 'notation_cds', 'chromosome', 'start_genomic', 'end_genomic']].drop_duplicates()
    
    #merge by order 'ensembl_enst'
    mutation = mutation.merge(protein_uniprot)
    print(f'rows mapped by ENST_name {mutation.shape[0]}')
    
    #some id_protein can get more than one transcript
    #kept only the most mutated transcript, ...
    y = mutation[['enst', 'id_protein']]
    dd = y.groupby(['id_protein', 'enst']).size().reset_index(name='counts')
    dd = dd.sort_values(by=['id_protein', 'counts', 'enst'], ascending=[False, False, True])
    dd = dd.drop_duplicates(subset=['id_protein'], keep='first')
    dd = dd[['id_protein', 'enst']]
    mutation = mutation.merge(dd)
    print(f'rows mapped by only one ENST by id_protein {mutation.shape[0]}') 
    
    #analizar enst con disitino id_protein
    y =  mutation[['enst', 'id_protein']]
    y = y.drop_duplicates()
    dd = y.groupby(['enst']).size().reset_index(name='counts')
    print('transcript to more than one id_protein should be change later')
    dd = dd[dd.counts > 1]
    print(dd)
    print(y[y['enst'].isin(dd['enst'])])
    #este es un caso raro ya que la secuencia del transcripto es la union de dos uniprots
    #transcriptos de la misma ENST00000527548 u sus transcriptos equivalentes ENST00000529639, ENST00000531743
    #P35544 + P62861
    #id_protein 123 y 4462    

    print(f'total mutations {mutation.shape[0]}')  
    return mutation, other_data

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
        df_crop['end_aa'] = df_crop['aux'].map(lambda x: x[3] if (len(x) > 3 and x[3] != '') else x[1]) # poner en el end el start si no hay end

    # consecuencia de la mutacion
    df_crop['consequence'] = conseq

    df_crop = df_crop.drop(columns=['aux'])
    return df_crop

def put_mutation_consequence(mutation):
    mutation['index_line'] = list(range(0, len(mutation)))
        
    # A DataFrame for each molecular consequence
    synonym = separar_en_cols(mutation, "notation_aa", "synonym", "=$")   
    mutation = mutation[~mutation['index_line'].isin(synonym['index_line'].tolist())]
    
    aux_other = separar_en_cols(mutation, "notation_aa", "synonym", "^p\.([\*A-Z])+(\d+)=$", override=True)   
    mutation = mutation[~mutation['index_line'].isin(aux_other['index_line'].tolist())]
    synonym = pd.concat([synonym, aux_other])
    print(f"Found {synonym.shape[0]} synonyms")
    
    delins = separar_en_cols(mutation, "notation_aa", "delins", "delins")
    mutation = mutation[~mutation['index_line'].isin(delins['index_line'].tolist())]
    print(f"Found {delins.shape[0]} delins")
    
    deletions = separar_en_cols(mutation, "notation_aa", "deletion", "del$") # finish with 'del' or "del$"
    mutation = mutation[~mutation['index_line'].isin(deletions['index_line'].tolist())]
    print(f"Found {deletions.shape[0]} deletions")
    
    insertions = separar_en_cols(mutation, "notation_aa", "insertion", "ins")
    mutation = mutation[~mutation['index_line'].isin(insertions['index_line'].tolist())]
    print(f"Found {insertions.shape[0]} inserions")
    
    frameshift = separar_en_cols(mutation, "notation_aa", "frameshift", 'fs') #or fs\*
    mutation = mutation[~mutation['index_line'].isin(frameshift['index_line'].tolist())]
    print(f"Found {frameshift.shape[0]} frameshift")
    
    duplications = separar_en_cols(mutation, "notation_aa", "duplication", "dup")
    mutation = mutation[~mutation['index_line'].isin(duplications['index_line'].tolist())]
    print(f"Found {duplications.shape[0]} duplications")
    
    nostop = separar_en_cols(mutation, "notation_aa", "nostop", "(?:del)?ext") #or '(?:del)?ext\*'
    mutation = mutation[~mutation['index_line'].isin(nostop['index_line'].tolist())]
    
    aux_other = separar_en_cols(mutation, "notation_aa", "nostop", "^p\.([\*A-Z]+)(\d+)(?:del)?ext(?:.*)$", override=True)   
    mutation = mutation[~mutation['index_line'].isin(aux_other['index_line'].tolist())]
    nostop = pd.concat([nostop, aux_other])    
    print(f"Found {nostop.shape[0]} nostop")    
    
    nonsense = separar_en_cols(mutation, "notation_aa", "nonsense", "\*$")
    mutation = mutation[~mutation['index_line'].isin(nonsense['index_line'].tolist())]
    print(f"Found {nonsense.shape[0]} nonsense")    
    
    missense = separar_en_cols(mutation, "notation_aa", "missense", "^p\.([\*A-Z])(\d+)([\*A-Z])$", override=True)
    mutation = mutation[~mutation['index_line'].isin(missense['index_line'].tolist())]
    print(f"Found {missense.shape[0]} missense")
    
    repeted = separar_en_cols(mutation, "notation_aa", "repeted", "^p\.(\d+)_?(\d+)?([A-Z]*)\[(\d+)\]$", override=True)   
    mutation = mutation[~mutation['index_line'].isin(repeted['index_line'].tolist())]
    print(f"Found {repeted.shape[0]} repeted")
    
    print(f"No mapped mutations {mutation.shape[0]}")
    print(mutation['notation_aa'].unique().tolist())
    ## Concatenate everything

    tables = [synonym, deletions, delins, duplications, frameshift, insertions, missense, nonsense, nostop, repeted]
    mutation = pd.concat(tables)
        
    mutation['start_aa'] = mutation['start_aa'].astype(int)
    mutation['end_aa'] = mutation['end_aa'].astype(int)
    print(f"Total mutations: {mutation.shape[0]}")
    mutation = mutation.drop(columns=['index_line'])
    return mutation

def remapping_idprotein_ENST00000527548(mutation):
    print('mutation in ENST00000527548 (id protein 123) kept only mutations from 1-74')
    print('mutation in ENST00000527548 (id protein 4462) kept only mutations from 75-133 and change the range to 1-59')
    #change the anotation of mutations for the transcript ENST00000527548,
    #reported of two proteins 4340 and 4342
    #original nuevo	idprotein	uniprot
    #1-74	 1-74 	123   P35544
    #75-133  1-59 	4462  P62861
    y_others = mutation[~mutation['id_protein'].isin([123, 4462])]  
    print(f'mutations not in idprotein 123 or 4462, {y_others.shape[0]}')  
    
    y_123 = mutation[mutation['id_protein'] == 123] 
    y_4462 = mutation[mutation['id_protein'] == 4462] 
    
    print(f'mutations in idprotein 123, {y_123.shape[0]}')
    y_123 = y_123[(y_123['start_aa'] < 75) & (y_123['end_aa'] < 75)]
    print(f'id protein 123, mutations with start and end < 75 are {y_123.shape[0]}')
    
    print(f'mutations in idprotein 4462, {y_4462.shape[0]}')
    y_4462 = y_4462[(y_4462['start_aa'] > 74) & (y_4462['end_aa'] > 74)]
    print(f'id protein 4462, mutations with start and end > 74 are {y_4462.shape[0]}')
        
    #change the start_aa and end_aa for 4342
    #75 should be the 1, and 133 should be 59
    y_4462['start_aa'] = y_4462['start_aa'] - 74
    y_4462['end_aa'] = y_4462['end_aa'] - 74
    #change the notation p. with the new mapping 1-59
    #print(y_4462['notation_aa'].tolist())
    y_4462['notation_aa'] = y_4462['notation_aa'].str.findall('^(p\.[\*A-Z])(\d+)(_)?([\*A-Z])?(\d+)?(.*)$').str[0]
    y_4462['notation_aa'] = y_4462['notation_aa'].map(lambda x: "".join([x[0], str(int(x[1]) - 74), x[2], x[3], str(int(x[4]) - 74) if x[4] != "" else "", x[5]])) 
    #print(y_4462['notation_aa'].tolist())    
    
    mutation = pd.concat([y_others, y_123, y_4462], ignore_index=True)
    mutation = mutation.drop_duplicates()
    print(f'total of mutations {mutation.shape[0]}')
    return mutation

if __name__ == '__main__':

    opts = parse_args()
    
    if not os.path.exists(opts.outfolder) or not os.path.isdir(opts.outfolder):
        os.mkdir(opts.outfolder)
    pd.set_option('display.max_columns', None)
    #return mutation columns: 'index_line', 'id_protein'   
    mutation, other_data = mapping_enst(opts)
    
    #return mutation columns: 'enst', 'id_protein', 'notation_aa', notation_cds, 'start_aa', 'end_aa', 'consequence'  
    mutation = put_mutation_consequence(mutation)
        
    #select the other_data
    mutation = mutation.merge(other_data)
    
    #re-mapping position of mutation for  ENST00000527548 proteins 
    mutation = remapping_idprotein_ENST00000527548(mutation)
    print(f'total sample-mutations {mutation.shape[0]}')
    y = mutation[["id_protein", "notation_aa", "notation_cds", 'chromosome', 'start_genomic', 'end_genomic']]
    y = y.drop_duplicates() 
    print(f'unique idprot-notation_aa_cds {y.shape[0]}, proteins {len(set(y["id_protein"].tolist()))}')
     
    mutation.to_csv(os.path.join(opts.outfolder, 'mutations_parse_cosmic.tsv.gz'), sep='\t', index= False, compression='gzip')

    
    
