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
    
    vs_all = pd.read_csv(arg_files.clinvar_file, sep='\t', compression='gzip')  
    vs_all['index_line'] = list(range(0, len(vs_all)))
    #merge at the end to put the data
    vs_ids_other = vs_all[['index_line', 'snp_id', 'variationid', 'chromosome', 'start', 'stop', 'type', 'cambio', 'cambio_nt', 'origin', 'phenotypeids', 'phenotypelist', 'otherids']]
    #only kept in vs_all the 'index_line' and the ids used to mapp
    vs_all = vs_all[['index_line', 'gene_name', 'gene_by_name', 'hgnc_id', 'gene_id', 'refseq']]
    
    print(f'clinvar {vs_all.shape[0]}')     
    coln = ['index_line', 'id_protein']
        
    all_protein = pd.read_csv(arg_files.prot_file, sep="\t")
    all_protein = all_protein[['id_protein', 'hgnc_id', 'uniprot_acc', 'gene_name', 'gene_id']]
    
    #merge by refseq
    refseq_nt = pd.read_csv(arg_files.refseq_file,  sep="\t")
    protein_uniprot = all_protein[['id_protein', 'uniprot_acc']].copy()
    protein_uniprot = protein_uniprot.merge(refseq_nt).drop(columns=['uniprot_acc', 'version'])
    vs_refseq = vs_all.copy()
    vs_refseq = vs_refseq.merge(protein_uniprot)
    vs_refseq = vs_refseq[coln]
    print(f'clinvar merge by ref_seq {vs_refseq.shape[0]}, unique {len(vs_refseq["index_line"].unique().tolist())}')
    print(vs_refseq[vs_refseq['index_line'].isin(vs_refseq[vs_refseq.duplicated(['index_line'])]['index_line'].unique().tolist())])
    
    vs_all = vs_all[~vs_all['index_line'].isin(vs_refseq['index_line'].tolist())]
    
    #merge by gene_id
    protein_geneid = all_protein[['id_protein', 'gene_id']].copy()
    protein_geneid['gene_id'] = protein_geneid['gene_id'].str.split(';')
    protein_geneid = protein_geneid.explode('gene_id')
    protein_geneid['gene_id'] = protein_geneid['gene_id'].str.strip()
    protein_geneid = protein_geneid[~protein_geneid['gene_id'].isnull()]
    protein_geneid = protein_geneid[protein_geneid['gene_id'] != '']
    vs_geneid = vs_all.copy()
    
    vs_geneid = vs_geneid[~vs_geneid['gene_id'].isnull()]
    vs_geneid['gene_id'] = [str(int(x)) for x in vs_geneid['gene_id'].tolist()]
    vs_geneid['gene_id'] = vs_geneid['gene_id'].str.strip()
    vs_geneid = vs_geneid[vs_geneid['gene_id'] != '']    
    vs_geneid = vs_geneid.merge(protein_geneid)
    vs_geneid = vs_geneid[coln]
    print(f'clinvar merge by gene_id {vs_geneid.shape[0]}, unique {len(vs_geneid["index_line"].unique().tolist())}')
    print(vs_geneid[vs_geneid['index_line'].isin(vs_geneid[vs_geneid.duplicated(['index_line'])]['index_line'].unique().tolist())])
    
    vs_all = vs_all[~vs_all['index_line'].isin(vs_geneid['index_line'].tolist())]
    
    #merge by hgnc_id
    protein_hgnc = all_protein[['id_protein', 'hgnc_id']].copy()
    protein_hgnc['hgnc_id'] = protein_hgnc['hgnc_id'].str.split(';')
    protein_hgnc = protein_hgnc.explode('hgnc_id')
    protein_hgnc['hgnc_id'] = protein_hgnc['hgnc_id'].str.strip()
    protein_hgnc = protein_hgnc[~protein_hgnc['hgnc_id'].isnull()]
    protein_hgnc = protein_hgnc[protein_hgnc['hgnc_id'] != '']
    vs_hgnc = vs_all.copy()
    
    vs_hgnc['hgnc_id'] = vs_hgnc['hgnc_id'].str.strip()
    vs_hgnc['hgnc_id'] = vs_hgnc['hgnc_id'].str.lstrip("HGNC:")
    vs_hgnc = vs_hgnc[~vs_hgnc['hgnc_id'].isnull()]
    vs_hgnc = vs_hgnc[vs_hgnc['hgnc_id'] != '']   
    vs_hgnc = vs_hgnc.merge(protein_hgnc)
    vs_hgnc = vs_hgnc[coln]
    print(f'clinvar merge by hgnc_id {vs_hgnc.shape[0]}, unique {len(vs_hgnc["index_line"].unique().tolist())}')
    
    vs_all = vs_all[~vs_all['index_line'].isin(vs_hgnc['index_line'].tolist())]
    
        #merge by gene_name
    protein_gene = all_protein[['id_protein', 'gene_name']].copy()
    protein_gene['gene_name'] = protein_gene['gene_name'].str.split(';')
    protein_gene = protein_gene.explode('gene_name')
    protein_gene['gene_name'] = protein_gene['gene_name'].str.strip()
    protein_gene = protein_gene[~protein_gene['gene_name'].isnull()]
    protein_gene = protein_gene[protein_gene['gene_name'] != '']
    
    #merge by gene_by_name
    protein_gene = protein_gene.rename(columns={'gene_name': 'gene_by_name'})
    vs_byname = vs_all.copy()
    vs_byname = vs_byname.merge(protein_gene)
    vs_byname = vs_byname[coln]
    print(f'clinvar merge by gene_by_name {vs_byname.shape[0]}, unique {len(vs_byname["index_line"].unique().tolist())}')
    
    vs_all = vs_all[~vs_all['index_line'].isin(vs_byname['index_line'].tolist())]
    
    
    protein_gene = protein_gene.rename(columns={'gene_by_name': 'gene_name'})
    vs_genename = vs_all.copy()
    vs_genename['gene_name'] = vs_genename['gene_name'].str.split(';')
    vs_genename = vs_genename.explode('gene_name')  
    vs_genename['gene_name'] = vs_genename['gene_name'].str.strip()
    print(f'\tclinvar explode gene_name {vs_genename.shape[0]}') 
    vs_genename = vs_genename.merge(protein_gene)
    vs_genename = vs_genename[coln]
    print(f'clinvar merge by gene_name {vs_genename.shape[0]}, unique {len(vs_genename["index_line"].unique().tolist())}')
    #one gene_name more than one idprotein
    print(vs_genename[vs_genename['index_line'].isin(vs_genename[vs_genename.duplicated(['index_line'])]['index_line'].unique().tolist())])
    
    vs_all = vs_all[~vs_all['index_line'].isin(vs_genename['index_line'].tolist())]   
    
    #concat vs_refseq, vs_genename, vs_byname, vs_hgnc, vs_geneid
    mut_cv = pd.concat([vs_refseq, vs_genename, vs_byname, vs_hgnc, vs_geneid], ignore_index=True)
    mut_cv = mut_cv.drop_duplicates()
    print(f'clinvar without duplicates {mut_cv.shape[0]}, unique {len(mut_cv["index_line"].unique().tolist())}') 
    
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
        df_crop['aux'] = df_crop[column].str.findall('^([A-Z][a-z]{2})(\d+)_?([A-Z][a-z]{2})?(\d+)?'+conseq_regex+'(.*)$').str[0]
    
    df_crop = df_crop[~df_crop['aux'].isnull()]
    
    if conseq == "repeted": 
        df_crop['start_aa'] = df_crop['aux'].map(lambda x: x[0])
        df_crop['end_aa'] = df_crop['aux'].map(lambda x: x[1] if x[1] != '' else x[0]) 
        df_crop['from_aa'] = df_crop['aux'].map(lambda x: x[2]) #one aa letter
        
    else:
        # start position
        df_crop['start_aa'] = df_crop['aux'].map(lambda x: x[1])
        # end position
        df_crop['end_aa'] = df_crop['aux'].map(lambda x: x[3] if x[3] != '' else x[1]) # poner en el end el start si no hay end
        # from: es el/los aa que cambian # concateno si existe mas de un aa que cambia (o sea, si es un rango)
        df_crop['from_aa'] = df_crop['aux'].map(lambda x: seq3to1(str(x[0])) + ('-' + seq3to1(str(x[2])) if x[2] != '' and x[3] != '' else ''))   
    
    if conseq == "missense":     
        df_crop['to_aa'] = df_crop['aux'].map(lambda x: seq3to1(x[2]))
    elif conseq == "nonsense":
        df_crop['to_aa'] = "*"
    elif conseq in ["frameshift", 'duplications', 'nostop', 'deletion', 'repeted', 'synonym']:
        df_crop['to_aa'] = "" #ver si poner las repeticiones para repeted
    elif conseq in ['delins', 'insertion']: # to: aa/s nuevos
        df_crop['to_aa'] = df_crop['aux'].map(lambda x: seq3to1(x[4]))
    
    # consecuencia de la mutacion
    df_crop['consequence'] = conseq

    df_crop = df_crop.drop(columns=['aux'])

    return df_crop


def separar_en_cols_cds(df, column, conseq, conseq_regex, override=False):
    '''
    recibe un DataFrame, el nombre de una columna auxiliar (column)
    y un string con el tipo de consecuencia (conseq). La col. auxiliar
    es una tupla con los elementos implicados en una mutacion
    como la siguiente (nucleatides1, nucleatides2).
    Devuelve el DataFrame df con estas 2 nuevas columnas
    '''
    
    df_crop = df.copy()
    if override:
        df_crop['aux'] = df_crop[column].str.findall(conseq_regex).str[0]
    else:
        df_crop['aux'] = df_crop[column].str.findall('^[-\+\*\d]+_?(?:[-\+\*\d]+)?([ACTG]+)?'+conseq_regex+'([ACTG]+)?$').str[0]
  
    df_crop = df_crop[~df_crop['aux'].isnull()]
    
    print(df_crop[['aux', 'cambio_nt']])
    df_crop['from_genomic'] = ''
    df_crop['to_genomic'] = ''
    
    if conseq in ["delins", 'insertion', 'missense']: 
        df_crop['to_genomic'] = df_crop['aux'].map(lambda x: x[1])
    if conseq == 'missense': 
        df_crop['from_genomic'] = df_crop['aux'].map(lambda x: x[0])

    # consecuencia de la mutacion
    df_crop['consequence_genomic'] = conseq
    #print(df_crop[df_crop['badformat'] == True].sort_values(by='AA', key=lambda col: col.str.len(), ascending=False, na_position='last'))
    
    df_crop = df_crop.drop(columns=['aux'])

    return df_crop

if __name__ == "__main__":

    opts = parse_args()

    if not os.path.exists(opts.outfolder) or not os.path.isdir(opts.outfolder):
        os.mkdir(opts.outfolder)
    
    proteins_clinvar_total = mapp_by_all_ids(opts)    
    print(f'clinvar mutations to mapp AA: {proteins_clinvar_total.shape[0]}')
    
    proteins_clinvar_total['index_line'] = list(range(0, len(proteins_clinvar_total)))
    
    # A DataFrame for each molecular consequence
    synonym = separar_en_cols(proteins_clinvar_total, "cambio", "synonym", "=$")   
    proteins_clinvar_total = proteins_clinvar_total[~proteins_clinvar_total['index_line'].isin(synonym['index_line'].tolist())]
    print(f"Found {synonym.shape[0]} synonyms")
    
    delins = separar_en_cols(proteins_clinvar_total, "cambio", "delins", "delins")
    proteins_clinvar_total = proteins_clinvar_total[~proteins_clinvar_total['index_line'].isin(delins['index_line'].tolist())]
    print(f"Found {delins.shape[0]} delins")
    
    deletions = separar_en_cols(proteins_clinvar_total, "cambio", "deletion", "del") # finish with 'del'
    proteins_clinvar_total = proteins_clinvar_total[~proteins_clinvar_total['index_line'].isin(deletions['index_line'].tolist())]
    print(f"Found {deletions.shape[0]} deletions")
    
    insertions = separar_en_cols(proteins_clinvar_total, "cambio", "insertion", "ins")
    proteins_clinvar_total = proteins_clinvar_total[~proteins_clinvar_total['index_line'].isin(insertions['index_line'].tolist())]
    print(f"Found {insertions.shape[0]} inserions")
    
    frameshift = separar_en_cols(proteins_clinvar_total, "cambio", "frameshift", 'fs') # expressions as 'Lys1254Terfs' es un frameshift donde el primer aa queda como *
    proteins_clinvar_total = proteins_clinvar_total[~proteins_clinvar_total['index_line'].isin(frameshift['index_line'].tolist())]
    print(f"Found {frameshift.shape[0]} frameshift")
    
    duplications = separar_en_cols(proteins_clinvar_total, "cambio", "duplication", "dup")
    proteins_clinvar_total = proteins_clinvar_total[~proteins_clinvar_total['index_line'].isin(duplications['index_line'].tolist())]
    print(f"Found {duplications.shape[0]} duplications")
    
    nonsense = separar_en_cols(proteins_clinvar_total, "cambio", "nonsense", "Ter$") # positiv lookbehind search! must have a number before, some delins insert a Ter
    proteins_clinvar_total = proteins_clinvar_total[~proteins_clinvar_total['index_line'].isin(nonsense['index_line'].tolist())]
    print(f"Found {nonsense.shape[0]} nonsense")
    
    nostop = separar_en_cols(proteins_clinvar_total, "cambio", "nostop", "ext")
    proteins_clinvar_total = proteins_clinvar_total[~proteins_clinvar_total['index_line'].isin(nostop['index_line'].tolist())]
    print(f"Found {nostop.shape[0]} nostop")
    
    missense = separar_en_cols(proteins_clinvar_total, "cambio", "missense", '')
    proteins_clinvar_total = proteins_clinvar_total[~proteins_clinvar_total['index_line'].isin(missense['index_line'].tolist())]
    print(f"Found {missense.shape[0]} missense")
    
    repeted = separar_en_cols(proteins_clinvar_total, "cambio", "repeted", "^(\d+)_?(\d+)?([A-Z]*)\[(\d+)\]$", override=True)   
    proteins_clinvar_total = proteins_clinvar_total[~proteins_clinvar_total['index_line'].isin(repeted['index_line'].tolist())]
    print(f"Found {repeted.shape[0]} repeted")
    
    print(f"No mapped mutations {proteins_clinvar_total.shape[0]}")
    print(proteins_clinvar_total['cambio'].unique().tolist())
    ## Concatenate everything

    tables = [synonym, deletions, delins, duplications, frameshift, insertions, missense, nonsense, nostop, repeted]
    mutations = pd.concat(tables)
    print(f"Total mutations: {mutations.shape[0]}")

    #split cds part
    proteins_clinvar_total = mutations
    print(f'clinvar mutations to mapp CDS: {proteins_clinvar_total.shape[0]}')
      
    delins = separar_en_cols_cds(proteins_clinvar_total, "cambio_nt", "delins", "delins")
    proteins_clinvar_total = proteins_clinvar_total[~proteins_clinvar_total['index_line'].isin(delins['index_line'].tolist())]
    print(f"Found {delins.shape[0]} delins")
    
    deletions = separar_en_cols_cds(proteins_clinvar_total, "cambio_nt", "deletion", "del") # finish with 'del'
    proteins_clinvar_total = proteins_clinvar_total[~proteins_clinvar_total['index_line'].isin(deletions['index_line'].tolist())]
    print(f"Found {deletions.shape[0]} deletions")
    
    insertions = separar_en_cols_cds(proteins_clinvar_total, "cambio_nt", "insertion", "ins")
    proteins_clinvar_total = proteins_clinvar_total[~proteins_clinvar_total['index_line'].isin(insertions['index_line'].tolist())]
    print(f"Found {insertions.shape[0]} inserions")
    
    duplications = separar_en_cols_cds(proteins_clinvar_total, "cambio_nt", "duplication", "dup")
    proteins_clinvar_total = proteins_clinvar_total[~proteins_clinvar_total['index_line'].isin(duplications['index_line'].tolist())]
    print(f"Found {duplications.shape[0]} duplications")
    
    inversions = separar_en_cols_cds(proteins_clinvar_total, "cambio_nt", "inversion", "inv")
    proteins_clinvar_total = proteins_clinvar_total[~proteins_clinvar_total['index_line'].isin(inversions['index_line'].tolist())]
    print(f"Found {inversions.shape[0]} inversions")
    
    missense = separar_en_cols_cds(proteins_clinvar_total, "cambio_nt", "missense", '>')
    proteins_clinvar_total = proteins_clinvar_total[~proteins_clinvar_total['index_line'].isin(missense['index_line'].tolist())]
    print(f"Found {missense.shape[0]} missense")
    
    print(f"No mapped mutations {proteins_clinvar_total.shape[0]}")
    print(proteins_clinvar_total['cambio_nt'].unique().tolist())
    ## Concatenate everything

    tables = [deletions, delins, duplications, inversions, insertions, missense]
    mutations = pd.concat(tables)
    print(f"Total mutations: {mutations.shape[0]}")
    
    # Assign id_mutation. A mutation is unique by the indexes : 'uniprot_acc', 'chromosome', 'start', 'stop', 'cambio', 'cambio_nt'
    subset = mutations[['id_protein', 'chromosome', 'start', 'stop', 'start_aa', 'end_aa', 'from_aa', 'to_aa', 'from_genomic', 'to_genomic',]].drop_duplicates()
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
