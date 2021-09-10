#%% Import dependecies
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from brokenaxes import brokenaxes

#%% Read data
protein = pd.read_csv('../db_tables/protein.tsv', sep='\t')

consequence = pd.read_csv('../db_tables/consequence.tsv', sep='\t')
mutation = pd.read_csv('../db_tables/mutation.tsv', sep='\t').merge(consequence).drop(columns='id_consequence')
rol = pd.read_csv('../db_tables/rol.tsv', sep='\t')
mlo = pd.read_csv('../db_tables/mlo.tsv', sep='\t')
protein_has_mlo = pd.read_csv('../db_tables/protein_has_mlo.tsv', sep='\t').merge(rol).merge(mlo).drop(columns=['id_rol', 'id_mlo'])
disease = pd.read_csv('../db_tables/disease.tsv', sep='\t')
mutation_has_disease = pd.read_csv('../db_tables/mutation_has_disease.tsv', sep='\t').merge(disease)

#%%
# Mutaciones unicas
len(mutation.id_mutation.unique()) # 169793
# En 2473 proteinas de LLPS
len(mutation.id_protein.unique())

# En cuantas disease?
len(mutation_has_disease.id_disease.unique()) # 4242

# Enfermedades con mas de 10 mutaciones
mutation_has_disease.name.value_counts()[mutation_has_disease.name.value_counts() >= 10] # 1231


disease_mutations = mutation_has_disease.groupby('name')['id_mutation'].count().sort_values(ascending= False)
others = pd.Series({'others': disease_mutations[10:].sum()})
topten = disease_mutations[:10]
topten_others= pd.concat([others, topten])

# consequences
cq = mutation.consequence.value_counts()
cq.rename(index={'frameshit': 'frameshift'}, inplace= True)
#[f'{s} ({s*100/sum(cq):.1f}%)' for s in cq]
pcts = [f'({s*100/sum(cq):.1f}%)' for s in cq]
labels = cq.index + " " + pcts
legend = cq.index + [f" {s}" for s in cq]
# Mutaciones por MLOs
#mut_mlo = protein.merge(mutation, on ='id_protein').merge(protein_has_mlo, on='id_protein')

#%% MUTATION BY CONSEQUENCE (no queda bien)
fig, ax = plt.subplots(figsize=(12, 8), subplot_kw=dict(aspect="equal"))

recipe = ["225 g flour",
          "90 g sugar",
          "1 egg",
          "60 g butter",
          "100 ml milk",
          "1/2 package of yeast"]

data = [225, 90, 50, 60, 100, 5]

wedges, texts = ax.pie(cq, wedgeprops=dict(width= 0.4), startangle= 40)

bbox_props = dict(boxstyle="square,pad=0.5", fc="w", ec="k", lw=2)
kw = dict(arrowprops=dict(arrowstyle="-"),
          bbox=bbox_props, zorder=1, va="center")

for i, p in enumerate(wedges):
    ang = (p.theta2 - p.theta1)/2. + p.theta1
    y = np.sin(np.deg2rad(ang))
    x = np.cos(np.deg2rad(ang))
    horizontalalignment = {-1: "right", 1: "left"}[int(np.sign(x))]
    connectionstyle = "angle,angleA=0,angleB={}".format(ang)
    kw["arrowprops"].update({"connectionstyle": connectionstyle})
    ax.annotate(labels[i], xy=(x, y), xytext=(1.4*np.sign(x), 1.35*y),
                horizontalalignment=horizontalalignment, **kw)

ax.set_title("Mutations by consequence\n", fontsize=25, weight= 'bold')

plt.show()

#%% MUTATION BY CONSEQUENCE (OK!)
sns.set_palette('terrain', 10)
fig, ax = plt.subplots(figsize=(12, 8), subplot_kw=dict(aspect="equal"))

# Create a circle at the center of
# the plot
my_circle = plt.Circle((0, 0), 0.7, color='white')

# Explode tuple
explode = (0.0, 0.0, 0.0, 0.15, 0.3, 0.44, 0.6) 

# Give color names
plt.pie(cq,  startangle= 40,explode= explode, labels= labels, textprops={'fontsize': 18})


p = plt.gcf()
p.gca().add_artist(my_circle)

# Add Legends
#plt.legend(legend, bbox_to_anchor=(0,1))

ax.set_title("Mutations by consequence", fontsize=25, weight= 'bold')

plt.savefig("mutations_consequence.png", format='png', dpi=500)
# Show the graph
plt.show()


#%% TOP 10 DISEASES (donut)
#sns.set_palette('turbo', 10)
fig, ax = plt.subplots(figsize=(12, 20), subplot_kw=dict(aspect="equal"))

pcts = [f'{s} {l}\n({s*100/sum(topten):.1f}%)' for s,l in zip(topten, topten.index)]

wedges, texts = ax.pie(topten, wedgeprops=dict(width=0.5), startangle=-20)

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

ax.set_title("Top ten diseases in DisPhaseDB\n", fontsize=25, weight='bold')
plt.show()

# %% Mutations by consequence Donut
#plt.style.use('seaborn-pastel')

pcts = [f'{s} \n({s*100/sum(cq):.1f}%)' for s in cq]
#width = 0.35

fig, ax = plt.subplots(figsize=(10, 6), subplot_kw=dict(aspect="equal"))

wedges, texts = ax.pie(cq, wedgeprops=dict(width=0.5), startangle= 25, labels= pcts) #, rotatelabels=True

bbox_props = dict(boxstyle="square,pad=0.3", fc="w", ec="k", lw=0.72)
kw = dict(arrowprops=dict(arrowstyle="-"), bbox=bbox_props, zorder=0, va="center")

for i, p in enumerate(wedges):
    ang = (p.theta2 - p.theta1)/2. + p.theta1
    y = np.sin(np.deg2rad(ang))
    x = np.cos(np.deg2rad(ang))
    horizontalalignment = {-1: "right", 1: "left"}[int(np.sign(x))]
    connectionstyle = "angle,angleA=0,angleB={}".format(ang)
    kw["arrowprops"].update({"connectionstyle": connectionstyle})
    ax.annotate(cq.index[i], xy=(x, y), xytext=(1.5*np.sign(x), 1.5*y),
                horizontalalignment=horizontalalignment, **kw)

ax.set_title("Mutations by consequence\n", fontsize=20, weight= 'bold')

plt.show()

# %% Mutations by consequence Bar (Sale bueno, pero usare donut)
# Horizontal barplot
sns.set_style("whitegrid")
bar,ax = plt.subplots(figsize=(14,8))
ax = sns.barplot(x= cq, y= cq.index, palette="turbo", orient='h')
ax.set_title("Mutations by consequence\n", fontsize=25, weight='bold')
ax.set_xlabel ("Number of mutations", fontsize=20, weight='bold')
ax.set_ylabel ("Consequence", fontsize=20, weight='bold')
ax.tick_params(labelsize=15)
for i, v in enumerate(cq):
    ax.text(v, i, str(v), weight='bold', fontsize=14)
plt.show()

#%% Top ten diseases barra (OK)
disease_mutations = mutation_has_disease.groupby('name')['id_mutation'].count().sort_values(ascending= False)
#others = pd.Series({'others': disease_mutations[10:].sum()})
topten = disease_mutations[:10]
# Horizontal barplot
sns.set_style("whitegrid")
sns.set_palette('turbo', 10)

bar,ax = plt.subplots(figsize=(14,10))
ax = sns.barplot(x= topten, y= topten.index, palette= 'rainbow',  orient='h')
ax.set_title("Top ten diseases by number of mutations\n", fontsize=30, weight='bold')
ax.set_xlabel ("Number of mutations", fontsize=25, weight='bold')
ax.set_ylabel ("Disease", fontsize=25, weight='bold')
ax.tick_params(labelsize=20)
for i, v in enumerate(topten):
    ax.text(v, i, str(v), weight='bold', fontsize=18)

plt.savefig("diseases_mutations.png", format='png', dpi=500, bbox_inches='tight')
plt.show()

#%% Top ten diseases barra (OK) - AGRUPAR
top = disease_mutations[:30].copy()
top.rename({
    'Hereditary cancer-predisposing syndrome':'Cancer',
    'Hereditary breast and ovarian cancer syndrome': 'Cancer',
    'Breast-ovarian cancer, familial 1': 'Cancer',
    'Breast-ovarian cancer, familial 2': 'Cancer',
    'Familial cancer of breast': 'Cancer',
    'Rhabdoid tumor predisposition syndrome 2': 'Cancer',
    'Gastrointestinal stromal tumor': 'Cancer',
    'DICER1-related pleuropulmonary blastoma cancer predisposition syndrome': 'Cancer',
    'Hereditary nonpolyposis colorectal neoplasms': 'Cancer',
    'Epidermolysis bullosa simplex with muscular dystrophy': 'Epidermolysis bullosa simplex',
    'Epidermolysis bullosa simplex with pyloric atresia': 'Epidermolysis bullosa simplex',
    'Epidermolysis bullosa simplex, Ogna type': 'Epidermolysis bullosa simplex',
    'Epidermolysis bullosa simplex with nail dystrophy': 'Epidermolysis bullosa simplex'
    }, inplace= True)
top = top.sum(level=0)[:10]
#%% Top ten diseases barra GROUPED (OK!!)
# Horizontal barplot
sns.set_style("whitegrid")
sns.set_palette('turbo')

f, (ax1, ax2) = plt.subplots(1, 2, sharey=True, figsize=(14,10))

sns.barplot(x= top, y= top.index, palette= 'rainbow',  orient='h', ax=ax1)
sns.barplot(x= top, y= top.index, palette= 'rainbow',  orient='h', ax=ax2)

# zoom-in / limit the view to different portions of the data
ax1.set_xlim(0, 6500) # most of the data
ax2.set_xlim(58000,59500) # outliers only

# Make the spacing between the two axes a bit smaller
plt.subplots_adjust(wspace=0.1)

ax = ax1
ax.set_title("                                           Top ten diseases by number of mutations\n", fontsize=30, weight='bold')
ax.set_xlabel ("                                                  Number of mutations", fontsize=25, weight='bold')
ax.set_ylabel ("Disease", fontsize=25, weight='bold')
ax.tick_params(labelsize=20)
ax2.set(ylabel= None, xlabel= None)
ax2.set(xticklabels=[58000, '', '', '', '', '', 59500])

ax2.text(top[0], 0, str(top[0]), weight='bold', fontsize=20)
for i, v in enumerate(top[1:]):
    ax1.text(v, i+1, str(v), weight='bold', fontsize=20)
    #ax2.text(v, i, str(v), weight='bold', fontsize=20)

#plt.savefig("diseases_mutations_1.png", format='png', dpi=500, bbox_inches='tight')
plt.show()

#%% PROBANDO 2... (no sirve)
# sns.set_style("whitegrid")
sns.set_palette('turbo')

fig = plt.figure(figsize=(14,10))
baxes = brokenaxes(xlims=((0, 7500), (58000,60000)), hspace= .15)
ax = sns.barplot(x= top, y= top.index, palette= 'rainbow',  orient='h')
ax.yaxis.grid(False) # Show the vertical gridlines



ax.set_title("Top ten diseases by number of mutations\n", fontsize=30, weight='bold')
ax.set_xlabel ("Number of mutations", fontsize=25, weight='bold')
ax.set_ylabel ("Disease", fontsize=25, weight='bold')
ax.tick_params(labelsize=20)
for i, v in enumerate(top):
    ax.text(v, i, str(v), weight='bold', fontsize=18)

#plt.savefig("diseases_mutations.png", format='png', dpi=500, bbox_inches='tight')
plt.show()
#%% PROBANDO



ax.set_title("Top ten diseases by number of mutations\n", fontsize=30, weight='bold')
ax.set_xlabel ("Number of mutations", fontsize=25, weight='bold')
ax.set_ylabel ("Disease", fontsize=25, weight='bold')
ax.tick_params(labelsize=20)
#for i, v in enumerate(top):
    #ax2.text(v, i, str(v), weight='bold', fontsize=18)

#plt.savefig("diseases_mutations.png", format='png', dpi=500, bbox_inches='tight')
plt.show()

#%% PROTEINS with most mutations
len(mutation.id_protein.unique()) # 2473
prot_mutation = mutation.merge(protein)[['id_mutation',	'snp_id', 'start_aa',	'end_aa',	'from_aa',	'to_aa',	'consequence',	'uniprot_acc',	'uniprot_name']]
# %% Group
prot = prot_mutation.groupby('uniprot_name')['id_mutation'].count().sort_values(ascending= False)
prot_others = pd.Series({'others': prot[10:].sum()})
prot_topten = prot[:10]
prot_topten_others= pd.concat([prot_others, prot_topten])
# %% PLOT top ten mutations in LLPS proteins
sns.set_style("whitegrid")
sns.set_palette("husl", 10)

bar,ax = plt.subplots(figsize=(7,14))
ax = sns.barplot(x= prot_topten, y= prot_topten.index, orient='h')
ax.set_title("Phase separation proteins with most mutations\nin DisPhaseDB\n", fontsize=30, weight='bold')
ax.set_xlabel ("Number of mutations", fontsize=25, weight='bold')
ax.set_ylabel ("Protein name", fontsize=25, weight='bold')
ax.tick_params(labelsize=20)
for i, v in enumerate(prot_topten):
    ax.text(v, i, str(v), weight='bold', fontsize=19)

plt.savefig("proteins_mutations.png", format='png', dpi=500, bbox_inches='tight')
plt.show()

# %% Cuantas mutaciones caen en regiones IDR?
mutation_has_disorder_region = pd.read_csv('../db_tables/mutation_has_disorder_region.tsv', sep='\t')
len(mutation_has_disorder_region.id_mutation.unique()) # 20550 mutaciones en IDRs
len(mutation_has_disorder_region.id_idr.unique()) # 2073 IDRs unicos con mutaciones

# %% Cuantas en LOW-COMPLEXITY?
mutation_has_low_complexity = pd.read_csv('../db_tables/mutation_has_low_complexity.tsv', sep='\t')
len(mutation_has_low_complexity.id_mutation.unique()) # 17369 mutaciones en LC
len(mutation_has_low_complexity.id_lc.unique()) # 3477 LC con mutaciones
# %% Cuantas en PFAM
mutation_has_pfam_domain = pd.read_csv('../db_tables/mutation_has_pfam_domain.tsv', sep='\t')
len(mutation_has_pfam_domain.id_mutation.unique()) # 80837 mutaciones en PFAM domains
len(mutation_has_pfam_domain.id_pfam.unique()) # 1639 PFAM domains con mutaciones

# %% Mutaciones en DOMINIOS PFAM
mutation_pfam = pd.read_csv('../db_tables/mutation_has_pfam_domain.tsv', sep='\t')

# Verificar si estan las mutaciones que se eliminaron
# Este es un subset son las que ya no van (en este caso 9)
#eliminar = mutation_pfam[~mutation_pfam.id_mutation.isin(mutations.id_mutation)]
# Con esto las elimino:
#mutation_pfam = mutation_pfam[~mutation_pfam.id_mutation.isin(eliminar.id_mutation)]

# %% Traer el nombre de los PFAM acc
pfam_domain = pd.read_csv('../db_tables/pfam_domain.tsv', sep='\t')
mutation_pfam = mutation_pfam.merge(pfam_domain)
# %% Tengo que agrupar 
# groupby the column to group on "name". Then, ["id_mutation"] to specify the col to perform the actual aggregation.
pfam = mutation_pfam.groupby('pfam_domain')['id_mutation'].count().sort_values(ascending= False)
# %%
pfam_others = pd.Series({'others': pfam[10:].sum()})
pfam_topten = pfam[:10]
pfam_topten_others= pd.concat([pfam_others, pfam_topten])
# %% PLOT Top Ten mutations in PFAM domains
# Horizontal barplot
sns.set_style("whitegrid")
sns.set_palette("cool", 10)

bar,ax = plt.subplots(figsize=(7,14))
ax = sns.barplot(x= pfam_topten, y= pfam_topten.index, orient='h')
ax.set_title("PFAM domains with most mutations\nin DisPhaseDB\n", fontsize=30, weight='bold')
ax.set_xlabel ("Number of mutations", fontsize=30, weight='bold')
ax.set_ylabel ("domain name", fontsize=25, weight='bold')
ax.tick_params(labelsize=22)
for i, v in enumerate(pfam_topten):
    ax.text(v, i, str(v), weight='bold', fontsize=18)

plt.savefig("pfam_mutations.png", format='png', dpi=500, bbox_inches='tight')
plt.show()
# %%
