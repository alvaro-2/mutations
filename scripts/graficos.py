import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
os.getcwd()
os.chdir('C:\\Users\\User\\Documents\\mutations')

consequence= pd.read_csv('db_tables/consequence.tsv', sep='\t')
mutations = pd.read_csv('db_tables\mutation.tsv', sep='\t').merge(consequence)

# Diseases - tabla sin las 32 mutaciones eliminadas
mutation_has_disease = pd.read_csv('db_tables\mutation_has_disease.tsv', sep='\t')

disease_mutations = mutation_has_disease.groupby('name')['id_mutation'].count().sort_values(ascending= False)
others = pd.Series({'others': disease_mutations[10:].sum()})
topten = disease_mutations[:10]
topten_others= pd.concat([others, topten])

# %% Mutations by consequence Donut
plt.style.use('seaborn-pastel')
cq = mutations.consequence.value_counts()
labels = cq.index
sizes = cq
pcts = [f'{s} \n({s*100/sum(sizes):.1f}%)' for s in sizes]
#width = 0.35

fig, ax = plt.subplots(figsize=(10, 6), subplot_kw=dict(aspect="equal"))

wedges, texts = ax.pie(sizes, wedgeprops=dict(width=0.5), startangle= 50, labels= pcts) #, rotatelabels=True

bbox_props = dict(boxstyle="square,pad=0.3", fc="w", ec="k", lw=0.72)
kw = dict(arrowprops=dict(arrowstyle="-"), bbox=bbox_props, zorder=0, va="center")

for i, p in enumerate(wedges):
    ang = (p.theta2 - p.theta1)/2. + p.theta1
    y = np.sin(np.deg2rad(ang))
    x = np.cos(np.deg2rad(ang))
    horizontalalignment = {-1: "right", 1: "left"}[int(np.sign(x))]
    connectionstyle = "angle,angleA=0,angleB={}".format(ang)
    kw["arrowprops"].update({"connectionstyle": connectionstyle})
    ax.annotate(cq.index[i], xy=(x, y), xytext=(2*np.sign(x), 1.5*y),
                horizontalalignment=horizontalalignment, **kw)

ax.set_title("Mutations by consequence\n", fontsize=20, weight= 'bold')

plt.show()

#%% Top Ten diseases - Pie chart
# groupby the column to group on"name". Then, ["id_mutation"] to specify the col to perform the actual aggregation.
topten_others.plot(kind='pie', ylabel='')
plt.title("Top ten diseases in DisPhaseDB", fontsize=20, weight= 'bold')
plt.show()

#%% Top ten diseases donut 
sns.set_palette('turbo', 10)
fig, ax = plt.subplots(figsize=(12, 8), subplot_kw=dict(aspect="equal"))

pcts = [f'{s} {l}\n({s*100/sum(topten):.1f}%)' for s,l in zip(topten, topten.index)]

wedges, texts = ax.pie(topten, wedgeprops=dict(width=0.5), startangle=-35)

bbox_props = dict(boxstyle="square,pad=0.3", fc="w", ec="k", lw=0.72)
kw = dict(arrowprops=dict(arrowstyle="-"), bbox=bbox_props, zorder=0, va="center")

for i, p in enumerate(wedges):
    ang = (p.theta2 - p.theta1)/2. + p.theta1
    y = np.sin(np.deg2rad(ang))
    x = np.cos(np.deg2rad(ang))
    horizontalalignment = {-1: "right", 1: "left"}[int(np.sign(x))]
    connectionstyle = "angle,angleA=0,angleB={}".format(ang)
    kw["arrowprops"].update({"connectionstyle": connectionstyle})
    ax.annotate(topten.index[i], xy=(x, y), xytext=(1.35*np.sign(x), 1.4*y),
                horizontalalignment=horizontalalignment, **kw, fontsize=15)

ax.set_title("Top ten diseases in DisPhaseDB", fontsize=25, weight= 'bold')
plt.show()

#%% Top ten diseases barra (incluyendo others)
# Horizontal barplot
sns.set_style("whitegrid")
sns.set_palette('turbo', 10)

bar,ax = plt.subplots(figsize=(14,8))
ax = sns.barplot(x= topten_others, y= topten_others.index, palette="turbo", orient='h')
ax.set_title("Top ten diseases in DisPhaseDB\n", fontsize=25, weight='bold')
ax.set_xlabel ("Number of mutations", fontsize=20, weight='bold')
ax.set_ylabel ("Consequence", fontsize=20, weight='bold')
ax.tick_params(labelsize=15)
for i, v in enumerate(topten_others):
    ax.text(v, i, str(v), weight='bold', fontsize=14)


# %% DOMINIOS
mutation_pfam = pd.read_csv('db_tables\mutation_has_pfam_domain.tsv', sep='\t')

# Verificar si estan las mutaciones que se eliminaron
# Este es un subset son las que ya no van (en este caso 9)
eliminar = mutation_pfam[~mutation_pfam.id_mutation.isin(mutations.id_mutation)]
# Con esto las elimino:
mutation_pfam = mutation_pfam[~mutation_pfam.id_mutation.isin(eliminar.id_mutation)]

# %% Traer el nombre de los PFAM acc
pfam_domain = pd.read_csv('db_tables\pfam_domain.tsv', sep='\t')
mutation_pfam = mutation_pfam.merge(pfam_domain)
# %% Tengo que agrupar 
# groupby the column to group on "name". Then, ["id_mutation"] to specify the col to perform the actual aggregation.
pfam = mutation_pfam.groupby('pfam_domain')['id_mutation'].count().sort_values(ascending= False)

# %%
#mutation_pfam.groupby(['id_mutation', 'id_protein']).size().sort_values(ascending=False).reset_index(name='count').drop_duplicates(subset='id_protein')
# %%
pfam_others = pd.Series({'others': pfam[10:].sum()})
pfam_topten = pfam[:10]
pfam_topten_others= pd.concat([pfam_others, pfam_topten])
# %% PLOT Top Ten mutations in PFAM domains
# Horizontal barplot
sns.set_style("whitegrid")
sns.set_palette("husl", 11)

bar,ax = plt.subplots(figsize=(14,8))
ax = sns.barplot(x= pfam_topten, y= pfam_topten.index, orient='h')
ax.set_title("PFAM domains with most mutations in DisPhaseDB\n", fontsize=25, weight='bold')
ax.set_xlabel ("Number of mutations", fontsize=20, weight='bold')
ax.set_ylabel ("PFAM domain name", fontsize=20, weight='bold')
ax.tick_params(labelsize=15)
for i, v in enumerate(pfam_topten):
    ax.text(v, i, str(v), weight='bold', fontsize=14)

plt.show()
# %% DISORDER REGIONS: no tienen nombres. Como identificarlas mas alla del id numerico?

#%% PROTEINS with most mutations
len(mutations.id_protein.unique()) # 2473

# Traer los nombres de las prote
proteins = pd.read_csv('db_tables\protein_sel.tsv', sep='\t', usecols=['id_protein', 'uniprot_acc', 'uniprot_name'])

# %%
prot_mutation = mutations.merge(proteins)[['id_mutation',	'snp_id', 'start_aa',	'end_aa',	'consequence',	'uniprot_acc',	'uniprot_name']]
# %% Group
prot = prot_mutation.groupby('uniprot_name')['id_mutation'].count().sort_values(ascending= False)
prot_others = pd.Series({'others': prot[10:].sum()})
prot_topten = prot[:1000]
prot_topten_others= pd.concat([prot_others, prot_topten])
# %% PLOT top ten mutations in LLPS proteins
#sns.set_style("whitegrid")
sns.set_palette("husl", 10)

bar,ax = plt.subplots(figsize=(14,8))
ax = sns.scatterplot(x= prot_topten.index, y= prot_topten)
ax.set_title("Phase Separation proteins with most mutations in DisPhaseDB\n", fontsize=25, weight='bold')
ax.set_ylabel ("Number of mutations", fontsize=20, weight='bold')
#ax.set_xlabel ("Protein name", fontsize=10, weight='bold')
#ax.tick_params(labelsize=10, rotation=90)
# for i, v in enumerate(prot_topten):
#     ax.text(v, i, str(v), weight='bold', fontsize=14)

plt.show()
# %%
