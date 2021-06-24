import pandas as pd
import numpy as np
import re
import argparse
import os

# python preproc_clinvar.py --in ./variant_summary.txt.gz --out vs.csv.gz

def parse_args():
    parser = argparse.ArgumentParser(description='Pre-process clinvar file')

    parser.add_argument('--in',
                        dest='input',
                        help='Clinvar variant_summary.txt.gz')

    parser.add_argument('--out', 
                        dest='output',
                        help='Parsed output file')

    opts = parser.parse_args()
    return opts

    parser.parse_args()

if __name__ == "__main__":
    
    opts = parse_args()
    basedir = os.path.dirname(opts.input)
    if not os.path.exists(basedir):
        os.makedirs(basedir)
    clinvar = pd.read_csv(opts.input, sep= '\t')
    clinvar.columns = clinvar.columns.str.lower().str.replace(' ',"_").str.replace("-",'_').str.replace('/','_')
    clinvar = clinvar.rename(columns={'#alleleid' : 'alleleid', 'rs#_(dbsnp)': 'snpid', 'origin':'allelic_origins', 'originsimple':'origin', 'geneid': 'gene_id', 'genesymbol': 'gene_name'})

    # Subset GRCh38 genome assembly
    clinvar = clinvar[clinvar.assembly == 'GRCh38']
    vs = clinvar[['hgnc_id', 'snpid', 'variationid', 'chromosome', 'start', 'stop', 'type', 'name', 'origin', 'phenotypeids', 'phenotypelist', 'otherids']]
    vs.to_csv(opts.output, index= False, compression= 'gzip')

    # Subset mutations with "p." only
    vs['cambio'] = vs.name.map(lambda x: re.findall('\(p\..*\)$', x)).str[0]
    vs.cambio = vs.cambio.str.strip('()')
    vs.cambio = vs.cambio.str.lstrip('p.') # se usa lstrip xq strip tambien saca las p del final 
    # drop the nans
    vs = vs[vs.cambio.notnull()]

    # Cambio en nucleotido
    vs['cambio_nt'] = vs.name.map(lambda x: re.findall(': ?c\.(.*) ?\(p\.', x)).str[0]
    nt_null = vs.name[vs.cambio_nt.isnull()]
    print(nt_null.head(10))

    vs.to_csv(opts.output, index= False, compression= 'gzip')