import pandas as pd
clinvar = pd.read_csv('variant_summary.txt.gz', sep= '\t')
clinvar.columns
pmid = pd.read_csv('var_citations.txt', sep='\t')
# AlleleID                       integer value as stored in the AlleleID field in ClinVar  (//Measure/@ID in the XML)
# VariationID                    The identifier ClinVar uses to anchor its default display. (in the XML,  //MeasureSet/@ID)
# rs                             rs identifier from dbSNP, null if missing
# nsv                            nsv identifier from dbVar, null if missing
# citation_source                The source of the citation, either PubMed, PubMedCentral, or the NCBI Bookshelf
# citation_id                    The identifier used by that source
pmid.columns
pmid.info()
pmid.VariationID.isnull().any() # False



# Subset GRCh38 genome assembly
clinvar = clinvar[clinvar.Assembly == 'GRCh38'] # me quedo solo con las entradas del GRCh38
clinvar.VariationID.isnull().any()


clinvar_pmid = clinvar.merge(pmid)



clinvar.columns = clinvar.columns.str.lower().str.replace(' ',"_").str.replace("-",'_').str.replace('/','_')
clinvar = clinvar.rename(columns={'#alleleid' : 'alleleid', 'rs#_(dbsnp)': 'snpid', 'origin':'allelic_origins', 'originsimple':'origin'})

# Subset GRCh38 genome assembly
clinvar = clinvar[clinvar.assembly == 'GRCh38'] # me quedo solo con las entradas del GRCh38
vs = clinvar[['geneid', 'genesymbol', 'hgnc_id', 'snpid', 'alleleid', 'chromosomeaccession', 'chromosome', 'start', 'stop', 'type', 'name', 'origin', 'phenotypeids', 'phenotypelist', 'otherids']]
vs.to_csv(opts.output, index= False, compression= 'gzip')