import pandas as pd
import argparse
import os


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
    clinvar = clinvar.rename(columns={'#alleleid' : 'alleleid', 'rs#_(dbsnp)': 'snpid', 'origin':'allelic_origins', 'originsimple':'origin'})

    # Subset GRCh38 genome assembly
    clinvar = clinvar[clinvar.assembly == 'GRCh38'] # me quedo solo con las entradas del GRCh38
    vs = clinvar[['geneid', 'genesymbol', 'hgnc_id', 'snpid', 'alleleid', 'chromosomeaccession', 'chromosome', 'start', 'stop', 'type', 'name', 'origin', 'phenotypeids', 'phenotypelist', 'otherids']]
    vs.to_csv(opts.output, index= False, compression= 'gzip')