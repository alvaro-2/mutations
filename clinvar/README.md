00 -> clinvar_box1.ipynb; no mirar, la version mejorada y completa es el de abajo

01 -> clinvar_box1_mutaciones_en_proteinas_v2.ipynb; este script procesa el dataset de ClinVar, vs (Variant Summary).
Requiere descargar el dataset vs.csv.gz de aqui: https://drive.google.com/file/d/1SXbzWKRhBeCvp8-FFrsi-mhnvtNRjyL6/view?usp=sharing

o se puede generar a partir de la database completa asi:

Clinvar database completa -> variant_summary.txt.gz, descargar de aqui:
https://drive.google.com/file/d/1voGrP3aeA5JUdJK4xdUJvlFIGdNxDpu8/view?usp=sharing
```
clinvar = pd.read_csv('variant_summary.txt', sep= '\t')
clinvar.columns = clinvar.columns.str.lower().str.replace(' ',"_").str.replace("-",'_').str.replace('/','_')
clinvar = clinvar.rename(columns={'#alleleid' : 'alleleid', 'rs#_(dbsnp)': 'snpid', 'origin':'allelic_origins', 'originsimple':'origin'})
clinvar = clinvar[clinvar.assembly == 'GRCh38'] # me quedo solo con las entradas del GRCh38
#Subset del df de ClinVar: vs
vs = clinvar[['geneid', 'genesymbol', 'hgnc_id', 'snpid', 'alleleid', 'chromosomeaccession', 'chromosome', 'start', 'stop', 'type', 'name', 'origin', 'phenotypeids', 'phenotypelist', 'otherids']]
vs.to_csv('datasets/vs.csv.gz', index= False, compression= 'gzip')
```

# dataset falses: son los aa que no coinciden con la seq canonica de uniprot, al googlear algunos de esos snps encontras que son en otras isoformas.
  con el codigo de acceso a nuccore (ej. NM_001002295.2, esta en la tabla) se podria traer la seq de esa isoforma (creo)
