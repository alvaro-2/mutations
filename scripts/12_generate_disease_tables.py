import pandas as pd

disease = pd.read_csv("../raw_data/MRCONSO.RRF", dtype="str", sep="|", names=['cui', 'lat', 'ts', 'lui', 'stt', 'sui', 'ispref', 'aui', 'saui', 'scui', 'sdui', 'sab', 'tty', 'code', 'str', 'srl', 'suppress', 'cvf', 'other'])

disease = disease[(disease['lat'] == 'ENG')]
disease_synonyms = disease[['cui', 'str', 'sab', 'code']].copy()
disease_synonyms['str'] = disease_synonyms['str'].str.replace('"', "")
disease_synonyms = disease_synonyms.drop_duplicates(ignore_index=True)

#keept only one name per term
disease = disease[(disease['ts'] == 'P') & (disease['stt'] == 'PF') & (disease['ispref'] == 'Y')]
disease = disease[['cui', 'str']].drop_duplicates(ignore_index=True).rename(columns={'str': 'name'})

#disease clinvar
disdata = pd.read_csv("../raw_data/diseases_clinvar.tsv.gz", sep="\t", dtype="str", compression='gzip')
disdata = disdata[['id_mutation', 'phenotypeids']].drop_duplicates(ignore_index=True)

#split by | for each disease
disdata['phenotypeids'] = disdata['phenotypeids'].str.split("|")
disdata = disdata.explode('phenotypeids')
disdata['phenotypeids'] = disdata['phenotypeids'].str.strip()
disdata['phenotypeids'] = disdata['phenotypeids'].str.split(",")
disdata = disdata.explode('phenotypeids')
disdata['phenotypeids'] = disdata['phenotypeids'].str.strip()

disdata['cui'] = disdata['phenotypeids'].str.findall('^(.+):(.+)$').str[0]
disdata = disdata[~disdata['cui'].isnull()]

disdata['cui'] = disdata['cui'].map(lambda x: x[1] if x[0] == 'MedGen' else '') 

disdata = disdata[disdata['cui'] != '']
disdata = disdata[['id_mutation', 'cui']].drop_duplicates(ignore_index=True)
disdata = disdata[disdata['cui'].isin(disease['cui'].unique().tolist())]


#disease cosmic
dis_cosmic = pd.read_csv("../raw_data/cosmic_v94/classification.csv", dtype="str")
dis_cosmic = dis_cosmic.rename(columns={'SITE_PRIMARY_COSMIC': 'Primary_site', 
                                        'SITE_SUBTYPE1_COSMIC': 'Site_subtype_1', 
                                        'SITE_SUBTYPE2_COSMIC': 'Site_subtype_2', 
                                        'SITE_SUBTYPE3_COSMIC': 'Site_subtype_3', 
                                        'HISTOLOGY_COSMIC': 'Primary_histology', 
                                        'HIST_SUBTYPE1_COSMIC': 'Histology_subtype_1', 
                                        'HIST_SUBTYPE2_COSMIC': 'Histology_subtype_2', 
                                        'HIST_SUBTYPE3_COSMIC': 'Histology_subtype_3', 
                                        'NCI_CODE': 'NCI'})
dis_cosmic = dis_cosmic[['NCI', 'Primary_site', 'Site_subtype_1', 'Site_subtype_2', 'Site_subtype_3', 'Primary_histology', 'Histology_subtype_1', 'Histology_subtype_2', 'Histology_subtype_3']]
dis_cosmic = dis_cosmic.drop_duplicates(ignore_index=True)

mutation_has_disease = pd.read_csv("../raw_data/diseases_cosmic.tsv.gz", sep="\t", dtype="str", compression='gzip')
mutation_has_disease = mutation_has_disease.merge(dis_cosmic, on=['Primary_site', 'Site_subtype_1', 'Site_subtype_2', 'Site_subtype_3', 'Primary_histology', 'Histology_subtype_1', 'Histology_subtype_2', 'Histology_subtype_3'])

nci_to_cui = disease_synonyms[disease_synonyms['sab'] == 'NCI']
nci_to_cui = nci_to_cui[['cui', 'code']].drop_duplicates(ignore_index=True).rename(columns={'code': 'NCI'})
mutation_has_disease = mutation_has_disease.merge(nci_to_cui, on=['NCI'])
mutation_has_disease = mutation_has_disease[['id_mutation', 'cui']]

#all mutation and disease
mutation_has_disease = disdata.append(mutation_has_disease)
mutation_has_disease = mutation_has_disease.drop_duplicates(ignore_index=True)

#save the mutation_has_disease
print(f'Generating table mutation_has_disease.tsv, rows {mutation_has_disease.shape[0]}') 
mutation_has_disease.to_csv("../db_tables/mutation_has_disease.tsv", sep="\t", index=False)

#save the disease
disease = disease[disease['cui'].isin(mutation_has_disease['cui'].unique().tolist())]
disease = disease.rename(columns={'name':'disease_name'})
print(f'Generating table disease.tsv, rows {disease.shape[0]}') 
disease.to_csv("../db_tables/disease.tsv", sep="\t", index=False)

#save the disease_has_synonyms
disease_synonyms = disease_synonyms[disease_synonyms['cui'].isin(disease['cui'].unique().tolist())]

cross_reference = disease_synonyms[['sab']].copy().drop_duplicates(ignore_index=True)
cross_reference['id_cross'] = range(1, len(cross_reference) + 1)

disease_synonyms = disease_synonyms.merge(cross_reference, on="sab").drop(columns=['sab']).rename(columns={'str':'name'})

disease_synonyms = disease_synonyms.drop_duplicates(ignore_index=True)

cross_reference = cross_reference.rename(columns={'sab':'cross_name'})
print(f'Generating table cross_references.tsv, rows {cross_reference.shape[0]}') 
cross_reference.to_csv("../db_tables/cross_reference.tsv", sep="\t", index=False)

print(f'Generating table disease_has_synonyms.tsv, rows {disease_synonyms.shape[0]}') 
disease_synonyms.to_csv("../db_tables/disease_has_synonyms.tsv", sep="\t", index=False)

'''
de la tabla de deases names clinvar quedarnos con: 
Finding
Disease
Infectious disease
'''