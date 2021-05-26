import pandas as pd
import argparse
import re
import numpy as np

def parse_args():
    parser = argparse.ArgumentParser(description='Parse clinvar mutations and merge with target proteins.')

    parser.add_argument('--clinvar',
                        dest='clinvar_file',
                        help='pre-processed clinvar variants vs.csv.gz')

    parser.add_argument('--prot', 
                        dest='prot_file',
                        help='File with target proteins (box1_proteins.csv)')
                        
    parser.add_argument('--out', 
                        dest='outfile',
                        help='Output file')

    opts = parser.parse_args()
    return opts

    parser.parse_args()


def seq3(seq):
    
    protein_letters_1to3 = {
        "A": "Ala",
        "C": "Cys",
        "D": "Asp",
        "E": "Glu",
        "F": "Phe",
        "G": "Gly",
        "H": "His",
        "I": "Ile",
        "K": "Lys",
        "L": "Leu",
        "M": "Met",
        "N": "Asn",
        "P": "Pro",
        "Q": "Gln",
        "R": "Arg",
        "S": "Ser",
        "T": "Thr",
        "V": "Val",
        "W": "Trp",
        "Y": "Tyr",
        "B": "Asx",
        "X": "Xaa",
        "Z": "Glx",
        "J": "Xle",
        "U": "Sec",
        "O": "Pyl"
    }
    
    return "".join(protein_letters_1to3.get(aa, "Xaa") for aa in seq)


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
        df_crop['end_aa'] = np.nan
        df_crop['from'] = df_crop['aux'].map(lambda x: x[0])        
        df_crop['to'] = df_crop['aux'].map(lambda x: x[2])
    else:
    
        # start position
        df_crop['start_aa'] = df_crop['aux'].map(lambda x: x[1])

        # end position
        df_crop['end_aa'] = df_crop['aux'].map(lambda x: int(x[3]) if x[3] != '' else np.nan)

        # from: es el/los aa que cambian
        df_crop['from'] = df_crop['aux'].map(lambda x: x[0] + x[2]) # concateno si existe mas de un aa que cambia (o sea, si es un rango)

        # to: aa/s nuevos
        if conseq == "nonsense":
            df_crop['to'] = "Ter"
        else:
            df_crop['to'] = df_crop['aux'].map(lambda x: x[4] if x[4] != '' else np.nan)

    # consecuencia de la mutacion
    df_crop['consequence'] = conseq

    df_crop = df_crop.drop(columns=['aux'])

    return df_crop[['cambio', 'start_aa', 'end_aa', 'from', 'to', 'consequence']]

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
    df_crop['end_aa'] = df_crop['aux'].map(lambda x: int(x[1]) if x[1] != '' else np.nan)
    
    # from: es el/los aa que cambian
    df_crop['from'] = df_crop['aux'].map(lambda x: x[2]) # concateno si existe mas de un aa que cambia (o sea, si es un rango)
    df_crop['from'] = df_crop['from'].map(lambda x: seq3(x))
    df_crop['from'] = df_crop['from'].apply(str)
    
    df_crop['to'] = np.nan
        
    # consecuencia de la mutacion
    df_crop['consequence'] = conseq

    df_crop = df_crop.drop(columns=['aux'])

    return df_crop[['cambio', 'start_aa', 'end_aa', 'from', 'to', 'consequence']]




if __name__ == "__main__":

    opts = parse_args()

    vs = pd.read_csv(opts.clinvar_file)
    box = pd.read_csv(opts.prot_file)

    # Mergeo con el dataset completo
    box1_clinvar_total = box.merge(vs)

    # Para generar una col con los codigos NM... estos son los id de los transcriptos
    box1_clinvar_total['nuccore_id'] = box1_clinvar_total.name.map(lambda x: re.findall('[A-Z]{2}\_[0-9]+\.[0-9]*', x))
    box1_clinvar_total['nuccore_id'] = box1_clinvar_total.nuccore_id.str[0]

    # Subset mutations with "p." only
    box1_clinvar_total['cambio'] = box1_clinvar_total.name.map(lambda x: re.findall('\(p\..*\)$', x))
    box1_clinvar_total['cambio'] = box1_clinvar_total.cambio.str[0]
    box1_clinvar_total.cambio = box1_clinvar_total.cambio.str.strip('()')  # para quitar los parentesis
    box1_clinvar_total.cambio = box1_clinvar_total.cambio.str.lstrip('p.')                             # se usa lstrip xq strip tambien saca las p del final 

    # saco los nans -- > Qué son los nan? 9433 NaNs probably intronic mutations?
    box1_clinvar_total = box1_clinvar_total[box1_clinvar_total.cambio.notnull()]

    # separo las mutaciones sinonimas, no las tenemos en cuenta xq no surgen cambio en la proteina
    syn = box1_clinvar_total[box1_clinvar_total.cambio.str.endswith('=')]
    # y las borro
    cond = box1_clinvar_total.index.isin(syn.index) # es un array de bool
    box1_clinvar_total = box1_clinvar_total.drop(box1_clinvar_total[cond].index) # drop esas filas

    delins = separar_en_cols(box1_clinvar_total, "cambio", "delins", "delins")
    deletions = separar_en_cols(box1_clinvar_total, "cambio", "deletion", "del$") # finish with 'del'
    insertions = separar_en_cols(box1_clinvar_total, "cambio", "insertion", "(?<!del)ins") # Negative lookbehind search!
    frameshift = separar_en_cols(box1_clinvar_total, "cambio", "frameshit", "fs$")
    nonsense = separar_en_cols(box1_clinvar_total, "cambio", "nonsense", "(?<=\d)Ter") # positiv lookbehind search! must have a number before, some delins insert a Ter
    missense = separar_en_cols(box1_clinvar_total, "cambio", "missense", '^([A-Z][a-z]{2})(\d+)(?!Ter)([A-Z][a-z]{2})$', override=True)
    duplications = separar_en_cols(box1_clinvar_total, "cambio", "duplication", "dup")

    ix_targets_lists = [list(x) for x in [delins.index, deletions.index, insertions.index, frameshift.index, nonsense.index, missense.index, duplications.index]]
    ix_targets = [y for x in ix_targets_lists for y in x]
    cond = box1_clinvar_total.index.isin(ix_targets)                       # es un array de bool
    box1_clinvar_left = box1_clinvar_total.drop(box1_clinvar_total[cond].index)   # drop esas filas

    # elimino los '?' y los Nan
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
        delet2['end_aa'].iloc[i] = delet2.start_aa.iloc[i] + len(delet2['from'].iloc[i]) / 3 - 1
    deletions2 = pd.concat([deletions, delet2])

    inser = box1_clinvar_left[box1_clinvar_left.type == 'Insertion']
    insert2 = separar_en_cols_raros(inser, "cambio", "insertion", '^(\d+)_?(\d+)?([A-Za-z]*)', override=True)
    insertions2 = pd.concat([insertions, insert2])

    # Final round of cleaning up
    ix_targets_lists = [list(x) for x in [delins.index, deletions2.index, insertions2.index, frameshift.index, nonsense.index, missense.index, duplications2.index]]
    ix_targets = [y for x in ix_targets_lists for y in x]
    cond = box1_clinvar_total.index.isin(ix_targets)                       # es un array de bool
    box1_clinvar_leftovers = box1_clinvar_total.drop(box1_clinvar_total[cond].index)   # drop esas filas
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
        aa = len(leftovers2check.loc[i]["from"]) / 3
        l.append(length == aa)


    leftovers2check['is_del'] = l
    delet3 = leftovers2check[leftovers2check.is_del == True]
    delet3["consequence"] = "deletion"
    delet3 = delet3.drop(columns=["is_del"])
    deletions3 = pd.concat([deletions2, delet3])
    
    delins2 = leftovers2check[leftovers2check.is_del == False]
    delins2.consequence = "delins"
    delins2 = delins2.drop(columns=["is_del"])
    delins3 = pd.concat([delins, delins2])


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
    mutations['source'] = 'clinvar'
    mutations.to_csv(opts.outfile, index= False, compression='gzip')