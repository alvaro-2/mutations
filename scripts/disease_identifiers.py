import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

protein = pd.read_csv('../db_tables/protein.tsv', sep='\t')
mutation = pd.read_csv('../db_tables/mutation.tsv', sep='\t')
disease = pd.read_csv('../db_tables/disease.tsv', sep='\t')
disease_has_synonyms = pd.read_csv('../db_tables/disease_has_synonyms.tsv', sep='\t')
mutation_has_disease = pd.read_csv('../db_tables/mutation_has_disease.tsv', sep='\t')
cross_reference = pd.read_csv('../db_tables/cross_reference.tsv', sep='\t')

mutation_has_disease.merge(disease).merge(disease_has_synonyms)


# Busco las mutaciones para FUS (id 1)
fus_mut = mutation[mutation.id_protein == 1]
# solo aquellas con snpid
fus_mut = fus_mut[fus_mut.snp_id.notnull()]
# tomo los CUI para esa mutacion (p.R521C; rs121909668; 168326)
cuis = mutation_has_disease[mutation_has_disease.id_mutation == 168326]

# busco los sinonimos
terminos = disease_has_synonyms[disease_has_synonyms.cui.isin(cuis.cui)].merge(cross_reference).drop_duplicates()