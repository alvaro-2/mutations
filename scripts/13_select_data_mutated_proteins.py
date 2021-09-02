import pandas as pd
#generar los files
files_to_sel = ["protein", 
                "disorder_region", 
                "low_complexity", 
                "protein_has_pfam_domain", 
                "protein_has_mlo"]
#protein_llps, 
files_mutations = ["mutation"]

#cargar cada file de files_mutations, quedarnos solo con la columna id_protein
protein_mutated = set()
for x in files_mutations:
    y = pd.read_csv("../db_tables/" + x + ".tsv", sep="\t")
    protein_mutated = protein_mutated.union(set(y['id_protein'].unique()))

protein_mutated = list(protein_mutated)
print(f'Mutated proteins {len(protein_mutated)}')
#en cada fichero de files_to_sel, eliminar las filas donde id_protein no sean protein_mutated
for x in files_to_sel:
    y = pd.read_csv("../db_tables/" + x + ".tsv", sep="\t")
    print(f'Table {x}.tsv rows {y.shape[0]}')
    y = y[y['id_protein'].isin(protein_mutated)]
    print(f'Table {x}_sel.tsv rows {y.shape[0]}')
    #guardar cada fichero como el nombre concatenado con _sel.tsv, sep="\t"
    y.to_csv("../db_tables/" + x + "_sel.tsv", sep="\t", index=False)
 
