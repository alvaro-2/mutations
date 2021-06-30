import pandas as pd
import argparse
import re
import numpy as np
import os

# python 03_parse_clinvar.py --clinvar ../raw_data/vs.tsv.gz --refseq ../raw_data/uniprot-refseq-nucleotido.tsv --prot ../db_tables/protein.tsv --out ../raw_data

def parse_args():
    parser = argparse.ArgumentParser(description='Parse clinvar mutations and merge with target proteins.')

    parser.add_argument('--clinvar',
                        dest='clinvar_file',
                        help='pre-processed clinvar variants (vs.tsv.gz)')

    parser.add_argument('--refseq', 
                        dest='refseq_file',
                        help='File with uniprot vs refseq (uniprot-refseq-nucleotido.tsv)')
    
    parser.add_argument('--prot', 
                        dest='prot_file',
                        help='File with target proteins (protein.tsv)')
                        
    parser.add_argument('--out', 
                        dest='outfolder',
                        help='Output file')

    opts = parser.parse_args()
    return opts

    parser.parse_args()

def mapp_by_all_ids(arg_files):
    all_protein = pd.read_csv(arg_files.prot_file, sep="\t")
    all_protein = all_protein[['id_protein', 'hgnc_id', 'uniprot_acc', 'gene_name', 'gene_id']]
    
    protein_gene = all_protein[['id_protein', 'gene_name']].copy()
    protein_gene['gene_name'] = protein_gene['gene_name'].str.split(';')
    protein_gene = protein_gene.explode('gene_name')
    protein_gene['gene_name'] = protein_gene['gene_name'].str.strip()
    protein_gene = protein_gene[~protein_gene['gene_name'].isnull()]
    protein_gene = protein_gene[protein_gene['gene_name'] != '']
    
    vs_all = pd.read_csv(arg_files.clinvar_file, sep='\t', compression='gzip')  
    vs_all['index_line'] = list(range(0, len(vs_all)))
    #merge at the end to put the data
    vs_ids_other = vs_all[['index_line', 'snp_id', 'variationid', 'chromosome', 'start', 'stop', 'type', 'cambio', 'cambio_nt', 'origin', 'phenotypeids', 'phenotypelist', 'otherids']]
    #only kept in vs_all the 'index_line' and the ids used to mapp
    vs_all = vs_all[['index_line', 'gene_name', 'gene_by_name', 'hgnc_id', 'gene_id', 'refseq']]
    
    print(f'clinvar {vs_all.shape[0]}') 
    
    coln = ['index_line', 'id_protein']
    
    vs1 = vs_all.copy()
    #merge by gene_name
    vs1['gene_name'] = vs1['gene_name'].str.split(';')
    vs1 = vs1.explode('gene_name')  
    vs1['gene_name'] = vs1['gene_name'].str.strip()
    print(f'clinvar explode gene_name {vs1.shape[0]}') 
    vs1 = vs1.merge(protein_gene)
    vs1 = vs1[coln]
    print(f'clinvar merge by gene_name {vs1.shape[0]}')
    
    #merge by gene_by_name
    protein_gene = protein_gene.rename(columns={'gene_name': 'gene_by_name'})
    vs2 = vs_all.copy()
    vs2 = vs2.merge(protein_gene)
    vs2 = vs2[coln]
    print(f'clinvar merge by gene_by_name {vs2.shape[0]}')
    
    #merge by hgnc_id
    protein_hgnc = all_protein[['id_protein', 'hgnc_id']].copy()
    protein_hgnc['hgnc_id'] = protein_hgnc['hgnc_id'].str.split(';')
    protein_hgnc = protein_hgnc.explode('hgnc_id')
    protein_hgnc['hgnc_id'] = protein_hgnc['hgnc_id'].str.strip()
    protein_hgnc = protein_hgnc[~protein_hgnc['hgnc_id'].isnull()]
    protein_hgnc = protein_hgnc[protein_hgnc['hgnc_id'] != '']
    vs3 = vs_all.copy()
    
    vs3['hgnc_id'] = vs3['hgnc_id'].str.strip()
    vs3['hgnc_id'] = vs3['hgnc_id'].str.lstrip("HGNC:")
    vs3 = vs3[~vs3['hgnc_id'].isnull()]
    vs3 = vs3[vs3['hgnc_id'] != '']   
    vs3 = vs3.merge(protein_hgnc)
    vs3 = vs3[coln]
    print(f'clinvar merge by hgnc_id {vs3.shape[0]}')
    
    #juntar por gene_id
    protein_geneid = all_protein[['id_protein', 'gene_id']].copy()
    protein_geneid['gene_id'] = protein_geneid['gene_id'].str.split(';')
    protein_geneid = protein_geneid.explode('gene_id')
    protein_geneid['gene_id'] = protein_geneid['gene_id'].str.strip()
    protein_geneid = protein_geneid[~protein_geneid['gene_id'].isnull()]
    protein_geneid = protein_geneid[protein_geneid['gene_id'] != '']
    vs4 = vs_all.copy()
    
    vs4 = vs4[~vs4['gene_id'].isnull()]
    vs4['gene_id'] = [str(int(x)) for x in vs4['gene_id'].tolist()]
    vs4['gene_id'] = vs4['gene_id'].str.strip()
    vs4 = vs4[vs4['gene_id'] != '']    
    vs4 = vs4.merge(protein_geneid)
    vs4 = vs4[coln]
    print(f'clinvar merge by gene_id {vs4.shape[0]}')
    
    #juntar por refseq
    refseq_nt = pd.read_csv(arg_files.refseq_file,  sep="\t")
    protein_uniprot = all_protein[['id_protein', 'uniprot_acc']].copy()
    protein_uniprot = protein_uniprot.merge(refseq_nt).drop(columns=['uniprot_acc', 'version'])
    vs5 = vs_all.copy()
    vs5 = vs5.merge(protein_uniprot)
    vs5 = vs5[coln]
    print(f'clinvar merge by ref_seq {vs5.shape[0]}')
       
    #juntar vs, rest_vs, rest_vs2, rest_vs3
    mut_cv = pd.concat([vs1, vs2, vs3, vs4, vs5], ignore_index=True)
    mut_cv = mut_cv.drop_duplicates()
    print(f'clinvar without duplicates {mut_cv.shape[0]}') 
    
    print(f'clinvar proteins {len(mut_cv["id_protein"].unique())}')
    mut_cv = mut_cv.merge(vs_ids_other).drop(columns=['index_line'])
    return mut_cv

def seq3to1(seq):
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
    seq2 = re.findall('[A-Z][a-z]{2}', seq)  
    return "".join(protein_letters_3to1.get(aa, "X") for aa in seq2)

def delete_brackets(x):
    if re.search('\[\d+\]', x):
        z = re.search('(\[\d+\])', x)
        mlen = len(z.groups()[0])
        return x[:-mlen]
    else:
        return x

def separar_en_cols(df, column, conseq, conseq_regex, override=False):
    '''
    recibe un DataFrame, el nombre de una columna auxiliar (column)
    y un string con el tipo de consecuencia (conseq). La col. auxiliar
    es una tupla con los elementos implicados en una mutacion
    como la siguiente (aa1, start_pos, aa2, end_pos, aa/s_nuevos).
    Devuelve el DataFrame df con estas 5 nuevas columnas
    '''
    df_crop = df[df[column].str.contains(conseq_regex)].copy()
      
    if override:
        df_crop['aux'] = df_crop[column].str.findall(conseq_regex).str[0]
    else:
        df_crop['aux'] = df_crop[column].str.findall('^([A-Z][a-z]{2})(\d+)_?([A-Z][a-z]{2})?(\d+)?'+conseq_regex+'(.*)$').str[0]
        
    if conseq == "missense":
        df_crop['start_aa'] = df_crop['aux'].map(lambda x: x[1])
        df_crop['end_aa'] = df_crop['aux'].map(lambda x: x[1])  # mismo start y end para los que no tienen un end
        df_crop['from'] = df_crop['aux'].map(lambda x: seq3to1(x[0]))       
        df_crop['to'] = df_crop['aux'].map(lambda x: seq3to1(x[2]))
    else:
    
        # start position
        df_crop['start_aa'] = df_crop['aux'].map(lambda x: x[1])

        # end position
        df_crop['end_aa'] = df_crop['aux'].map(lambda x: int(x[3]) if x[3] != '' else int(x[1])) # poner en el end el start si no hay end

        # from: es el/los aa que cambian
        df_crop['from'] = df_crop['aux'].map(lambda x: seq3to1(x[0]) + seq3to1(x[2])) # concateno si existe mas de un aa que cambia (o sea, si es un rango)

        # to: aa/s nuevos
        if conseq == "nonsense":
            df_crop['to'] = "*"
        else:
            df_crop['to'] = df_crop['aux'].map(lambda x: seq3to1(x[4]) if x[4] != '' else np.nan)

    # consecuencia de la mutacion
    df_crop['consequence'] = conseq

    df_crop = df_crop.drop(columns=['aux'])

    return df_crop

def separar_en_cols_raros(df, column, conseq, conseq_regex, override=False):
    df_crop = df[df[column].str.contains(conseq_regex)].copy()
      
    if override:
        df_crop['aux'] = df_crop[column].str.findall(conseq_regex).str[0]
    else:
        df_crop['aux'] = df_crop[column].str.findall('^([A-Z][a-z]{2})(\d+)_?([A-Z][a-z]{2})?(\d+)?'+conseq_regex+'(.*)$').str[0]
            
    # start position
    df_crop['start_aa'] = df_crop['aux'].map(lambda x: x[0])
    df_crop.start_aa = df_crop.start_aa.apply(int)
    
    # end position
    df_crop['end_aa'] = df_crop['aux'].map(lambda x: int(x[1]) if x[1] != '' else int(x[0]))
    
    if conseq != 'insertion':
        # from: es el/los aa que cambian
        df_crop['from'] = df_crop['aux'].map(lambda x: x[2]) # concateno si existe mas de un aa que cambia (o sea, si es un rango)
        #df_crop['from'] = df_crop['from'].map(lambda x: seq3to1(x)) # ya viene en codigo de 1 letra
        df_crop['from'] = df_crop['from'].apply(str)
        
        df_crop['to'] = np.nan
    else:
        # from: es el/los aa que se insertan
        df_crop['to'] = df_crop['aux'].map(lambda x: x[2]) 
        df_crop['to'] = df_crop['to'].apply(str)
        
        df_crop['from'] = np.nan

    # consecuencia de la mutacion
    df_crop['consequence'] = conseq

    df_crop = df_crop.drop(columns=['aux'])

    return df_crop


if __name__ == "__main__":

    opts = parse_args()

    if not os.path.exists(opts.outfolder) or not os.path.isdir(opts.outfolder):
        os.mkdir(opts.outfolder)
    
    proteins_clinvar_total = mapp_by_all_ids(opts)
    
    # separate synonyms mutations. We don't take them into account
    syn = proteins_clinvar_total[proteins_clinvar_total.cambio.str.endswith('=')]
    # delete them
    cond = proteins_clinvar_total.index.isin(syn.index) # bool array
    proteins_clinvar_total = proteins_clinvar_total.drop(proteins_clinvar_total[cond].index) # drop those rows
    print(f'Found {syn.shape[0]} synonym mutations')
    print(f'clinvar mutations to mapp: {proteins_clinvar_total.shape[0]}')
    
    # A DataFrame for each molecular consequence
    delins = separar_en_cols(proteins_clinvar_total, "cambio", "delins", "delins")
    deletions = separar_en_cols(proteins_clinvar_total, "cambio", "deletion", "del$") # finish with 'del'
    insertions = separar_en_cols(proteins_clinvar_total, "cambio", "insertion", "(?<!del)ins") # Negative lookbehind search!
    frameshift = separar_en_cols(proteins_clinvar_total, "cambio", "frameshit", '^([A-Z][a-z]{2})(\d+)_?(?!Ter)([A-Z][a-z]{2})?(\d+)?fs$(.*)$', override= True) # expressions as 'Lys1254Terfs' are nonsense, not frameshift
    nonsense = separar_en_cols(proteins_clinvar_total, "cambio", "nonsense", "(?<=\d)Ter") # positiv lookbehind search! must have a number before, some delins insert a Ter
    missense = separar_en_cols(proteins_clinvar_total, "cambio", "missense", '^([A-Z][a-z]{2})(\d+)(?!Ter)([A-Z][a-z]{2})$', override=True)
    duplications = separar_en_cols(proteins_clinvar_total, "cambio", "duplication", "dup")

    # Get indexes
    ix_targets_lists = [list(x) for x in [delins.index, deletions.index, insertions.index, frameshift.index, nonsense.index, missense.index, duplications.index]]
    ix_targets = [y for x in ix_targets_lists for y in x]
    cond = proteins_clinvar_total.index.isin(ix_targets)                                    # bool array
    box1_clinvar_left = proteins_clinvar_total.drop(proteins_clinvar_total[cond].index) # drop those rows

    # Drop the '?' and Nan
    box1_clinvar_left = box1_clinvar_left[(box1_clinvar_left.cambio != '?') & (box1_clinvar_left.cambio != '(?') & (box1_clinvar_left.cambio.notnull())]
    box1_clinvar_left.cambio = box1_clinvar_left.cambio.apply(delete_brackets)

    dup = box1_clinvar_left[box1_clinvar_left.type == 'Duplication']
    dup_raros = separar_en_cols_raros(dup, "cambio", "duplication", '^(\d+)_?(\d+)?([A-Za-z]*)', override=True)

    duplications2 = pd.concat([duplications, dup_raros])

    delet = box1_clinvar_left[box1_clinvar_left.type == 'Deletion']
    delet2 = separar_en_cols_raros(delet, "cambio", "deletion", '^(\d+)_?(\d+)?([A-Za-z]*)', override=True)

    # agrego la posicion de fin faltante (o más, busca Nan en el end_aa)
    ix = delet2.end_aa.isna()
    for i in np.where(ix)[0]:
        delet2['end_aa'].iloc[i] = delet2.start_aa.iloc[i] + len(delet2['from'].iloc[i]) - 1
    deletions2 = pd.concat([deletions, delet2])

    inser = box1_clinvar_left[box1_clinvar_left.type == 'Insertion']
    insert2 = separar_en_cols_raros(inser, "cambio", "insertion", '^(\d+)_?(\d+)?([A-Za-z]*)', override=True)
    insertions2 = pd.concat([insertions, insert2])

    # Final round of cleaning up
    ix_targets_lists = [list(x) for x in [delins.index, deletions2.index, insertions2.index, frameshift.index, nonsense.index, missense.index, duplications2.index]]
    ix_targets = [y for x in ix_targets_lists for y in x]
    cond = proteins_clinvar_total.index.isin(ix_targets)                       # es un array de bool
    box1_clinvar_leftovers = proteins_clinvar_total.drop(proteins_clinvar_total[cond].index)   # drop esas filas
    box1_clinvar_leftovers = box1_clinvar_leftovers[(box1_clinvar_leftovers.cambio != '?') & (box1_clinvar_leftovers.cambio != '(?') & (box1_clinvar_leftovers.cambio.notnull())]
    box1_clinvar_leftovers.cambio = box1_clinvar_leftovers.cambio.apply(delete_brackets)
    box1_clinvar_leftovers.shape

    leftovers2check = separar_en_cols_raros(box1_clinvar_leftovers, "cambio", "to_check", '^(\d+)_?(\d+)?([A-Za-z]*)', override=True )
    # si el rango de de las posiciones coincide con el nro de letras: es delecion
    l = []
    for i in leftovers2check.index:
        start = int(leftovers2check.loc[i]["start_aa"])
        try:
            end = int(leftovers2check.loc[i]["end_aa"])
        except:
            pass
        length = end - start + 1
        aa = len(leftovers2check.loc[i]["from"])
        l.append(length == aa)


    leftovers2check['is_del'] = l
    delet3 = leftovers2check[leftovers2check.is_del == True]
    delet3["consequence"] = "deletion"
    delet3 = delet3.drop(columns=["is_del"])
    deletions3 = pd.concat([deletions2, delet3])
    
    delins2 = leftovers2check[leftovers2check.is_del == False]
    delins2["consequence"] = "delins"
    delins2 = delins2.drop(columns=["is_del"])
    delins3 = pd.concat([delins, delins2])

    ix_targets_lists = [list(x) for x in [delins3.index, deletions3.index, insertions2.index, frameshift.index, nonsense.index, missense.index, duplications2.index]]
    ix_targets = [y for x in ix_targets_lists for y in x]
    cond = proteins_clinvar_total.index.isin(ix_targets)                       # es un array de bool
    box1_clinvar_leftovers_notmapped = proteins_clinvar_total.drop(proteins_clinvar_total[cond].index)   # drop esas filas
    print(f"Found {box1_clinvar_leftovers_notmapped.shape[0]} not mapped")
    print(box1_clinvar_leftovers_notmapped['cambio'])

    ## Concatenate everything
    print(f"Found {deletions3.shape[0]} deletions")
    print(f"Found {delins3.shape[0]} delins")
    print(f"Found {duplications2.shape[0]} duplications")
    print(f"Found {frameshift.shape[0]} frameshift")
    print(f"Found {insertions2.shape[0]} inserions")
    print(f"Found {missense.shape[0]} missense")
    print(f"Found {nonsense.shape[0]} nonsense")

    tables = [deletions3, delins3, duplications2, frameshift, insertions2, missense, nonsense]
    mutations = pd.concat(tables)
    print(f"Total mutations: {mutations.shape[0]}")

    # Assign id_mutation. A mutation is unique by the indexes : 'uniprot_acc', 'chromosome', 'start', 'stop', 'cambio', 'cambio_nt'
    subset = mutations[['id_protein', 'chromosome', 'start', 'stop', 'cambio', 'cambio_nt']].drop_duplicates()
    subset['id_mutation'] = range(1, len(subset)+1)
    print(f'Subset (unique mutations) length: {len(subset)}')
    print(f'All mutations (vs) before: {mutations.shape[0]}')
    mutations = mutations.merge(subset)
    print(f'All mutations after merge: {mutations.shape[0]}')

    # Another table for diseases
    diseases = mutations[['id_mutation', 'origin', 'phenotypeids', 'phenotypelist', 'otherids']].copy()
    diseases.to_csv(os.path.join(opts.outfolder, 'diseases.tsv.gz'), sep='\t', index= False, compression='gzip')
    
    # Add ClinVar source
    mutations_with_source = mutations[['id_mutation', 'variationid']].drop_duplicates()
    mutations_with_source['source'] = 'clinvar'
    mutations_with_source.to_csv(os.path.join(opts.outfolder, 'mutations_with_source.tsv.gz'), sep='\t', index= False, compression='gzip')
    
    # drop columns
    mutations.drop(columns=['variationid', 'origin', 'phenotypeids', 'phenotypelist', 'otherids'], inplace= True)
    # one file, one mutation
    mutations.drop_duplicates(ignore_index = True, inplace= True)
    print(f'Length mutations after drop duplicates: {mutations.shape[0]}')
    print(f'clinvar proteins {len(mutations["id_protein"].unique())}')
    mutations.to_csv(os.path.join(opts.outfolder, 'mutations.tsv.gz'), sep='\t', index= False, compression='gzip')
