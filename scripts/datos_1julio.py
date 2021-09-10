import os
from numpy.core.fromnumeric import sort
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pandas.core import groupby
import seaborn as sns
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

ax = mutations.consequence.value_counts().plot(kind='pie', autopct='%.0f%%')
ax.set_title("Mutations by consequence")
plt.legend()
plt.show()

sns.set()
#cq = mutations.groupby('consequence').size()
cq.plot(kind='pie', title= "Mutations by consequence", autopct=lambda p: '{:.0f}'.format((p/100)*cq.sum()))
plt.show()


################## Mutations by consequence barra ########################
cq = mutations.consequence.value_counts()
# Horizontal barplot
sns.set_style("whitegrid")
bar,ax = plt.subplots(figsize=(10,6))
ax = sns.barplot(x=cq, y= cq.index, palette="muted", orient='h')
ax.set_title("Mutations by consequence", fontsize=20, weight='bold')
ax.set_xlabel ("Number of mutations", fontsize=15, weight='bold')
ax.set_ylabel ("Consequence", fontsize=15, weight='bold')
for i, v in enumerate(cq):
    ax.text(v, i, str(v), weight='bold')

####################### Mutations by consequence donut #######################################
labels = cq.index
sizes = cq
pcts = [f'{s} {l}\n({s*100/sum(sizes):.1f}%)' for s,l in zip(sizes, labels)]
#width = 0.35

fig, ax = plt.subplots(figsize=(10, 6), subplot_kw=dict(aspect="equal"))

wedges, texts = ax.pie(sizes, wedgeprops=dict(width=0.5), startangle=-40, labels= pcts)

bbox_props = dict(boxstyle="square,pad=0.3", fc="w", ec="k", lw=0.72)
kw = dict(arrowprops=dict(arrowstyle="-"), bbox=bbox_props, zorder=0, va="center")

for i, p in enumerate(wedges):
    ang = (p.theta2 - p.theta1)/2. + p.theta1
    y = np.sin(np.deg2rad(ang))
    x = np.cos(np.deg2rad(ang))
    horizontalalignment = {-1: "right", 1: "left"}[int(np.sign(x))]
    connectionstyle = "angle,angleA=0,angleB={}".format(ang)
    kw["arrowprops"].update({"connectionstyle": connectionstyle})
    ax.annotate(cq.index[i], xy=(x, y), xytext=(1.35*np.sign(x), 1.4*y),
                horizontalalignment=horizontalalignment, **kw)

ax.set_title("Mutations by consequence", fontsize=20, weight= 'bold')

plt.show()


#################### Diseases ###########################
disease = pd.read_csv('db_tables/disease.tsv', sep='\t')
mutation_has_disease = pd.read_csv('db_tables/mutation_has_disease.tsv', sep='\t').merge(disease)

# Verificar si estan las mutaciones que se eliminaron
len(# estas son las que no van (32)
eliminar = mutation[~mutation.id_mutation.isin(mutations.id_mutation)]

mutation_has_disease = mutation_has_disease[~mutation_has_disease.id_mutation.isin(eliminar.id_mutation)]
#mutation_has_disease.to_csv('db_tables\mutation_has_disease_sel.tsv', sep='\t', index= False)


# groupby the column to group on"name". Then, ["id_mutation"] to specify the col to perform the actual aggregation.
disease_mutations = mutation_has_disease.groupby('name')['id_mutation'].count().sort_values(ascending= False)
topten = disease_mutations[:10]
topten.plot(kind='pie', ylabel='')
plt.title("Top ten diseases in DisPhaseDB", fontsize=20, weight= 'bold')
plt.show()

############### Top ten diseases donut #######################

fig, ax = plt.subplots(figsize=(10, 6), subplot_kw=dict(aspect="equal"))


wedges, texts = ax.pie(topten.index, wedgeprops=dict(width=0.5), startangle=-40)

bbox_props = dict(boxstyle="square,pad=0.3", fc="w", ec="k", lw=0.72)
kw = dict(arrowprops=dict(arrowstyle="-"), bbox=bbox_props, zorder=0, va="center")

for i, p in enumerate(wedges):
    ang = (p.theta2 - p.theta1)/2. + p.theta1
    y = np.sin(np.deg2rad(ang))
    x = np.cos(np.deg2rad(ang))
    horizontalalignment = {-1: "right", 1: "left"}[int(np.sign(x))]
    connectionstyle = "angle,angleA=0,angleB={}".format(ang)
    kw["arrowprops"].update({"connectionstyle": connectionstyle})
    ax.annotate(topten[i], xy=(x, y), xytext=(1.35*np.sign(x), 1.4*y),
                horizontalalignment=horizontalalignment, **kw)

ax.set_title("Top ten diseases in DisPhaseDB", fontsize=20, weight= 'bold')

plt.show()

################ Top ten diseases barra ############################
# Horizontal barplot
sns.set_style("whitegrid")
bar,ax = plt.subplots(figsize=(10,6))
ax = sns.barplot(x= topten, y= topten.index, palette="muted", orient='h')
ax.set_title("Top ten diseases in DisPhaseDB", fontsize=20, weight='bold')
ax.set_xlabel ("Number of mutations", fontsize=15, weight='bold')
ax.set_ylabel ("Consequence", fontsize=15, weight='bold')
for i, v in enumerate(topten):
    ax.text(v, i, str(v), weight='bold')



################# Ejemplo ###################
fig, ax = plt.subplots(figsize=(6, 3), subplot_kw=dict(aspect="equal"))

recipe = ["225 g flour",
          "90 g sugar",
          "1 egg",
          "60 g butter",
          "100 ml milk",
          "1/2 package of yeast"]

data = [225, 90, 50, 60, 100, 5]

wedges, texts = ax.pie(data, wedgeprops=dict(width=0.5), startangle=-40)

bbox_props = dict(boxstyle="square,pad=0.3", fc="w", ec="k", lw=0.72)
kw = dict(arrowprops=dict(arrowstyle="-"),
          bbox=bbox_props, zorder=0, va="center")

for i, p in enumerate(wedges):
    ang = (p.theta2 - p.theta1)/2. + p.theta1
    y = np.sin(np.deg2rad(ang))
    x = np.cos(np.deg2rad(ang))
    horizontalalignment = {-1: "right", 1: "left"}[int(np.sign(x))]
    connectionstyle = "angle,angleA=0,angleB={}".format(ang)
    kw["arrowprops"].update({"connectionstyle": connectionstyle})
    ax.annotate(recipe[i], xy=(x, y), xytext=(1.35*np.sign(x), 1.4*y),
                horizontalalignment=horizontalalignment, **kw)

ax.set_title("Matplotlib bakery: A donut")

plt.show()
