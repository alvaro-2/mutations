# preproc_uniprot.py

import pandas as pd
import numpy as np
import re
import argparse
import os
#os.chdir('G:\My Drive\FIL\project')


def parse_args():
    parser = argparse.ArgumentParser(description= 'Pre-process UniProt humsavar.txt dataset')

    parser.add_argument('--in',
                        dest= 'input',
                        help= 'UniProt humsavar.txt')
    
    parser.add_argument('--out',
                        dest= 'output',
                        help= 'Processed output file humsavar.tsv')

    opts = parser.parse_args()
    return opts

    parser.parse_args()

if __name__ == '__main__':

    opts = parse_args()
    basedir = os.path.dirname(opts.input)
    if not os.path.exists(basedir):
        os.makedirs(basedir)

    '''humsavar.txt (https://www.uniprot.org/docs/humsavar.txt):
    Index of manually curated Human polymorphisms and disease mutations from UniProtKB/Swiss-Prot.
    This file lists all missense variants annotated in UniProtKB/Swiss-Prot human
    entries. It provides a variant classification which is intended for research
    purposes only, not for clinical and diagnostic use.
    - The column 'Variant category' shows the classification of the variant using
    the American College of Medical Genetics and Genomics/Association for
    Molecular Pathology (ACMG/AMP) terminology (Richards et al. PubMed:25741868)
    into the following categories:
    
    LP/P = likely pathogenic or pathogenic
    LB/B = likely benign or benign
    US   = uncertain significance

    These categories are assigned based on the variant annotation in the
    corresponding UniProtKB/Swiss-Prot entries that is curated from literature
    reports. The classification may change over time and must not be considered
    as a definitive statement about the pathogenic role of a variant.

    - The column 'Disease name' shows the name of the disease or the disease sample
    in which variants have been found. Names are only provided for diseases
    catalogued in OMIM and for cancer samples.
    '''
    with open(opts.input) as f:
        gene_name=[]
        uniprot=[]
        ft_id=[]
        change=[]
        category=[]
        snp_id=[]
        disease_name=[]
        for line in f:
            stripped_line = line.strip()
            gene_name.append(stripped_line[0:10])
            uniprot.append(stripped_line[10:21])
            ft_id.append(stripped_line[21:33])
            change.append(stripped_line[33:48])
            category.append(stripped_line[48:57])
            snp_id.append(stripped_line[57:72])
            disease_name.append(stripped_line[72:])

    # Create the dataframe: humsavar
    humsavar = pd.DataFrame(list(zip(gene_name, uniprot, ft_id, change, category, snp_id, disease_name)),
                            columns=['gene_name', 'uniprot', 'ft_id', 'change', 'category', 'snp_id', 'disease_name'])
    humsavar = humsavar.drop([0]).reset_index(drop=True)

    # Drop blank spaces
    humsavar = humsavar.applymap(lambda x: x.strip())

    # Add omim accession col
    humsavar['mim'] = humsavar.disease_name.map(lambda x: re.findall('\[(.*?)\]', x))
    humsavar['mim'] = humsavar.mim.str[0]

    #l= re.findall('(.*?) ?(\[.*?\])?$', 'Parietal foramina 2 (PFM2) [MIM: 05255256]')
    #l
    #l[0][0]

    # disease col: get everithing except the omim accession
    humsavar['disease'] = humsavar.disease_name.map(lambda x: re.findall('(.*?) ?(\[.*?\])?$', x)) # el primer grupo atrapa cualquier cosa, luego viene un espacio, y el segundo atrapa el codigo [MIM], este o no
    humsavar['disease'] = humsavar.disease.str[0].str[0]

    # replace "-" by NaNs
    humsavar.disease.replace('-', np.nan, inplace= True)
    humsavar.snp_id.replace('-', np.nan, inplace= True)
    humsavar.gene_name.replace('-', np.nan, inplace= True)
    # Format change col
    humsavar.change.str.startswith('p.').value_counts() # 79376, all entries starts with 'p.'
    humsavar.change = humsavar.change.str.lstrip('p.')

    # humsavar[humsavar.disease_name != "-"] now I can drop disease_name col
    humsavar.drop(columns='disease_name', inplace= True)
    humsavar.rename(columns={'uniprot': 'uniprot_acc'}, inplace= True)

    #humsavar.change[humsavar.change.notnull()] # 79376, all entries have a protein change
    print(f'Uniprot total non-null entries: {humsavar.shape[0]}')

    #humsavar.info()
    # Save
    #humsavar.to_csv('humsavar.csv', index=False)
    humsavar.to_csv(opts.output, index=False, sep= '\t')
    print(f'Saved in raw_data as {opts.output}')