import pandas as pd
import numpy as np
import re
import argparse
import os

# python 02_preproc_clinvar.py --in ../raw_data/variant_summary.txt.gz --out ../raw_data/vs.tsv.gz

def parse_args():
    parser = argparse.ArgumentParser(description='Pre-process clinvar file')

    parser.add_argument('--in',
                        dest='input',
                        help='Clinvar variant_summary.txt.gz')

    parser.add_argument('--out', 
                        dest='output',
                        help='Parsed output file vs.tsv.gz')

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

    clinvar = clinvar.rename(columns={'#alleleid' : 'alleleid', 'rs#_(dbsnp)': 'snp_id', 'origin':'allelic_origins', 'originsimple':'origin', 'geneid': 'gene_id', 'genesymbol': 'gene_name'})
    # Subset GRCh38 genome assembly
    clinvar = clinvar[clinvar.assembly == 'GRCh38']
    vs = clinvar[['hgnc_id', 'gene_name', 'gene_id', 'snp_id', 'variationid', 'chromosome', 'start', 'stop', 'type', 'name', 'origin', 'phenotypeids', 'phenotypelist', 'otherids']]

    # Subset mutations with "p." only: change in protein
    vs['cambio'] = vs.name.map(lambda x: re.findall('\(p\.(.*)\)$', x))
    vs['cambio'] = vs.cambio.str[0]
    #vs.cambio = vs.cambio.str.strip('()')
    #vs.cambio = vs.cambio.str.lstrip('p.')
    # drop the nans
    vs = vs[vs.cambio.notnull()]
    # drop the proteins as '(?)' or '?', '=' because don't have the protein notation    
    vs = vs[~vs.cambio.isin(['(?)', '?', '='])]

    # Change in nucleotide
    vs['cambio_nt'] = vs.name.map(lambda x: re.findall(': *c\.(.*) ?\(p\.', x))
    vs['cambio_nt'] = vs.cambio_nt.str[0]
    
    # Create a nuccore id col (transcripts accession)
    #para las proteinas comienza con NM_
    vs['nuccore_id'] = vs.name.map(lambda x: re.findall('[A-Z]{2}\_[0-9]+\.[0-9]*', x))
    vs['nuccore_id'] = vs.nuccore_id.str[0]
    vs['refseq'] = vs['nuccore_id'].str.split(".").str[0]
    vs['refseq_version'] = vs['nuccore_id'].str.split(".").str[1]
    
    # Create a gene_by_name col (gene obtained from field name)
    vs['gene_by_name'] = vs.name.map(lambda x: re.findall('\(?([-_A-Za-z0-9]+)\)?: *c\.', x))
    vs['gene_by_name'] = vs.gene_by_name.str[0]
    
    # Replace missing values by NaNs
    vs.hgnc_id.replace(['-', ''], np.nan, inplace= True)      
    vs.gene_name.replace(['-', ''], np.nan, inplace= True)  
    vs.snp_id.replace(-1, np.nan, inplace= True)
    vs.gene_id.replace(-1, np.nan, inplace= True)    
    
    #put the rs before in snpid number
    vs['snp_id'] = vs.snp_id.map(lambda x: 'rs' + str(int(x)) if not np.isnan(x) else x)

    vs.to_csv(opts.output, index= False, compression= 'gzip', sep="\t")
