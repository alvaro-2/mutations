vda_gene.drop_duplicates()
import pandas as pd

disnames = pd.read_csv('../raw_data/disease_names_clinvar.txt', sep= '\t')
disdata = pd.read_csv("../raw_data/diseases_clinvar.tsv.gz", sep="\t", compression='gzip')

y = disdata[['phenotypeids', 'phenotypelist']]
y = y.drop_duplicates()
z = y.groupby(['phenotypelist']).size().reset_index(name="counts")
z = z[z['counts'] > 1]

#hay 21 phenotypelist con disitinto texto en phenotypeids 
#debido a que uno es un poco mas corto pues le falta algun identificador
#ordenaramenos la tabla por la longitud del campo phenotypeids,
#asi cuando generemos los datos para la disease tendra todos los ids
disdata['size_phenotypeids'] = disdata['phenotypeids'].str.len()
disdata = disdata.sort_values(by='size_phenotypeids', ascending=False)

#buscar el duplicado "Colon cancer, susceptibility to"  y COLON CANCER, SUSCEPTIBILITY TO
# y remplazar el texto
disdata['phenotypelist'] = disdata['phenotypelist'].str.replace('COLON CANCER, SUSCEPTIBILITY TO', 'Colon cancer, susceptibility to', regex=False, case=True)

#split by | for each disease
disdata['phenotypeids'] = disdata['phenotypeids'].str.split("|")
disdata['phenotypelist'] = disdata['phenotypelist'].str.split("|")
disdata.head()

#make the relationship of disease name with disease id
#and mutation with disease name
mutation_has_disease = [] 
ph_all = {}

for index, row in disdata.iterrows():  
    ph_id = [a.split(';') for a in row['phenotypeids']]
    ph_name = [a.split(';') if len(ph_id[index]) > 1 else [a] for index, a in enumerate(row['phenotypelist'])]
    #aplanar la lista
    ph_id = [y for x in ph_id for y in x]
    ph_name = [y for x in ph_name for y in x]
    found_d = [ph_all.get(x, None) for x in ph_name]   
    add_not = {ph_name[k]: [a.split(":", 1) for a in ph_id[k].split(",") ]
               for k in range(0, len(found_d)) if found_d[k] == None}
    #kept only the list with two values [db, id] and when db not is 'Gene'
    add_not = {k: [s for s in v if len(s) == 2 and s[0] != 'Gene'] 
               for k, v in add_not.items()}
    #add not exist phenotype
    ph_all.update(add_not)
      
    #add diseases to the mutation
    mutation_has_disease.extend([[row['id_mutation'], x] for x in ph_name])

# lineas con ; en list y no en ids.
# bad = [(4977, 1), (5366, 2), (35852, 1), (35853, 0), (35854, 0), (35855, 0), (35856, 0), (35857, 0), (35991, 0), (42615, 2), (42800, 3), (43119, 1), (138229, 0), (166250, 2), (171075, 5), (171307, 4), (171489, 2), (171885, 3), (202204, 0), (244515, 1), (246025, 0), (246026, 0), (246027, 0), (246028, 0), (246029, 0), (246030, 0), (246031, 0), (246032, 0), (246033, 0), (246034, 0), (261693, 1), (264349, 1), (264350, 0), (264351, 0), (264352, 0), (264441, 1), (264442, 0)]

#falsas enfermedades 'not provided', 'not specified', '10 conditions', etc...
#eliminar las falsas enfermedades de mut_dis y ph_all
y = list(set(list(ph_all.keys())))
y.sort() 
print(y)

#manually selected as no valid disease
no_valid_disease = ['not provided', 'not specified', '', '-', '10 conditions', '11 conditions', '12 conditions', 
                    '13 conditions', '14 conditions', '15 conditions', '16 conditions', '17 conditions', '18 conditions', 
                    '19 conditions', '20 conditions', '21 conditions', '22 conditions', '24 conditions', '27 conditions', 
                    '36 conditions', '6 conditions', '7 conditions', '8 conditions', '9 conditions']

ph_data = {key:ph_all[key] for key in no_valid_disease}
#only two mutations (mutation id 8897 and 249100) has the disease '' and that have the {'OMIM': '114500'}

#removing not valid disease
ph_all = {k:v for k, v in ph_all.items() if k not in no_valid_disease}
#total of 4375 disease

#disease without other ids
ph_no_ids = {k:v for k, v in ph_all.items() if len(v) == 0}
# 501 diseases whithout other ids

#generar tabla cross_reference
#tiene el id_cross, name, ...
cross_reference = pd.DataFrame({'id_cross': list(range(1,8)),
                                'name': ['MedGen', 'OMIM', 'MeSH', 'MONDO', 'Human Phenotype Ontology', 'Orphanet', 'EFO'],
                                'url': ['https://www.ncbi.nlm.nih.gov/medgen/', 'https://www.omim.org/', 'https://www.ncbi.nlm.nih.gov/mesh/', 'https://www.ebi.ac.uk/ols/ontologies/mondo/', 'https://hpo.jax.org/', 'https://www.orpha.net/','https://www.ebi.ac.uk/ols/ontologies/efo']
                               })
#save the disease 
print(f'Generating table cross_reference.tsv, rows {cross_reference.shape[0]}')
cross_reference.to_csv("../db_tables/cross_reference.tsv", sep="\t", index=False)

#generar tabla disease
#tiene el id_disease, name, ...
disease = list(ph_all.keys())
disease = pd.DataFrame({'id_disease': list(range(1, len(disease) + 1)),
                       'name': disease})
#save the disease
print(f'Generating table disease.tsv, rows {disease.shape[0]}') 
disease.to_csv("../db_tables/disease.tsv", sep="\t", index=False)

#generar tabla disease_has_cross
disease = disease.rename(columns={'name': 'disease'})
#aplanarlo para ponerle el la disease y hacer el data frame
ph_plane = [[k] + e for k, l in ph_all.items() for e in l]

disease_has_cross_reference = pd.DataFrame(ph_plane, columns =['disease', 'name', 'id_incross'])

#add the disease and cross_reference id to disease_has_cross
disease_has_cross_reference = disease_has_cross_reference.merge(disease)
disease_has_cross_reference = disease_has_cross_reference.merge(cross_reference)
disease_has_cross_reference = disease_has_cross_reference.drop(['disease', 'name', 'url'], axis=1)

#save the disease_has_hpo
print(f'Generating table disease_has_cross_reference.tsv, rows {disease_has_cross_reference.shape[0]}') 
disease_has_cross_reference.to_csv("../db_tables/disease_has_cross_reference.tsv", sep="\t", index=False)

#generar tabla mutation_has_disease
mutation_has_disease = pd.DataFrame(mutation_has_disease, columns =['id_mutation', 'disease'])
#eliminar los nombres duplicados
mutation_has_disease = mutation_has_disease.drop_duplicates()

#add the disease id to mutation_has_disease
mutation_has_disease = mutation_has_disease.merge(disease)
mutation_has_disease = mutation_has_disease.drop(['disease'], axis=1)

#save the mutation_has_disease
print(f'Generating table mutation_has_disease.tsv, rows {mutation_has_disease.shape[0]}') 
mutation_has_disease.to_csv("../db_tables/mutation_has_disease.tsv", sep="\t", index=False)


# de la tabla de deases names clinvar quedarnos con: 
#Finding
#Disease
#Infectious disease

disnames.columns = disnames.columns.str.lower().str.replace('#', '')
disnames.info()

# List with diseases names in clinvar
disnames.diseasename.isin(disnames)