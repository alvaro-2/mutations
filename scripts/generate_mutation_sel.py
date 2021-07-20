import os
from numpy.core.fromnumeric import sort
import pandas as pd
import numpy as np

os.getcwd()
os.chdir('C:\\Users\\User\\Documents\\mutations')

protein_sel = pd.read_csv('db_tables/protein_sel.tsv', sep='\t')
#protein = pd.read_csv('db_tables/protein.tsv', sep='\t')
consequence= pd.read_csv('db_tables/consequence.tsv', sep='\t')
mutation = pd.read_csv('db_tables/mutation.tsv', sep='\t').merge(consequence)
#mutation = mutation.merge(protein)
#mutation_insertion = mutation[mutation.consequence == 'insertion'].copy()
#mutation = mutation[['id_mutation', 'start_aa',	'end_aa', 'from_aa', 'to_aa', 'consequence']].copy()


# Missense (substitution): Format: “prefix”“amino_acid”“position”“new_amino_acid”, e.g. p.(Arg54Ser)
missense = mutation[mutation.consequence == 'missense']

missense['change'] = np.nan
for i in missense.index:
    missense['change'][i] = missense.from_aa[i] + str(missense.start_aa[i]) + missense.to_aa[i]

# frameshit: Format: “prefix”“amino_acid”position”new_amino_acid”“fs”“Ter”“position_termination_site”, e.g. p.(Arg123LysfsTer34)
# frame shifts can also be described using a short format; p.Arg123fs,
# i.e. indicating the first amino acid changed, its position and “fs” without further detail.
frameshit = mutation[mutation.consequence == 'frameshit']

frameshit['change'] = np.nan
for i in frameshit.index:
    frameshit['change'][i] = frameshit.from_aa[i][0] + str(frameshit.start_aa[i]) + "fs"


# Nonsense (substitution): Format: “prefix”“amino_acid”“position”“new_amino_acid”, e.g. p.Trp24Ter(p.Trp24*)
nonsense = mutation[mutation.consequence == 'nonsense']

nonsense['change'] = np.nan
for i in nonsense.index:
    nonsense['change'][i] = nonsense.from_aa[i] + str(nonsense.start_aa[i]) + nonsense.to_aa[i]


# Deletions: Format: “prefix”“amino_acid(s)+position(s)_deleted”“del”, e.g. p.(Cys76_Glu79del)
deletion = mutation[mutation.consequence == 'deletion']

deletion['change'] = np.nan
for i in deletion.index:
    # several aminoacids
    if len(deletion.from_aa[i]) >= 2:
        deletion['change'][i] = deletion.from_aa[i][0] + str(deletion.start_aa[i]) + "_" + deletion.from_aa[i][-1] + str(deletion.end_aa[i]) + "del"
    # one aminoacid
    else:
        deletion['change'][i] = deletion.from_aa[i] + str(deletion.start_aa[i]) + "del"


# Insertion: “prefix”“amino_acids+positions_flanking”“ins”“inserted_sequence”, e.g. p.(Lys23_Leu24insArgSerGln)
insertion = mutation[mutation.consequence == 'insertion']

# Hay 11 inserciones que tienen NaN en from_aa. Descartar por ahora, retomar luego...
#insertion.sequence[10749][1172:1174]
#insertion.sequence[40463][78]
#insertion.sequence[46625][230:232]
#insertion.sequence[61625][1944]
#insertion.sequence[65679][1302:1304]

insertion = insertion[insertion.from_aa.notnull()]

insertion['change'] = np.nan
for i in insertion.index:
    # several aminoacids
    if len(insertion.from_aa[i]) >= 2:
        insertion['change'][i] = insertion.from_aa[i][0] + str(insertion.start_aa[i]) + "_" + insertion.from_aa[i][-1] + str(insertion.end_aa[i]) + "ins" + insertion.to_aa[i]


# Delins: “prefix”“amino_acid(s)+position(s)_deleted”“delins”“inserted_sequence”, e.g. p.(Arg123_Lys127delinsSerAsp)
delins = mutation[mutation.consequence == 'delins']
# Hay 21 delins que tienen NaN en to_aa. Descartar por ahora, retomar luego...
delins = delins[delins.to_aa.notnull()]

delins['change'] = np.nan
for i in delins.index:
    # several aminoacids
    if len(delins.from_aa[i]) >= 2:
        delins['change'][i] = delins.from_aa[i][0] + str(delins.start_aa[i]) + "_" + delins.from_aa[i][-1] + str(delins.end_aa[i]) + "delins" + delins.to_aa[i] 
    # one aminoacid
    else:
        delins['change'][i] = delins.from_aa[i] + str(delins.start_aa[i]) + "delins" + delins.to_aa[i]


# Duplication: “prefix”“amino_acid(s)+position(s)_duplicated”“dup”, e.g. p.(Cys76_Glu79dup)
duplication = mutation[mutation.consequence == 'duplication']

# VER esta entrada salen nans en from y to: NM_001429.4(EP300):c.7070_7096dup (p.Asn2357_Ala2365dup)
duplication['from_aa'][3589] = "NA"

duplication['change'] = np.nan
for i in duplication.index:
    # several aminoacids
    if len(duplication.from_aa[i]) >= 2:
        duplication['change'][i] = duplication.from_aa[i][0] + str(duplication.start_aa[i]) + "_" + duplication.from_aa[i][-1] + str(duplication.end_aa[i]) + "dup" 
    # one aminoacid
    else:
        duplication['change'][i] = duplication.from_aa[i] + str(duplication.start_aa[i]) + "dup"

tables = [missense, frameshit, nonsense, deletion, insertion, delins, duplication]

mutations = pd.concat(tables).sort_values('id_mutation')
#mutations.to_csv('db_tables\mutation_sel.tsv', sep='\t', index= False)
