# Este script es para un Venn diagram de 4 sets de las mutaciones

setwd("C:/Users/User/Documents/mutations/disease_plots/")
#getwd()

#install.packages('VennDiagram')
#install.packages('RColorBrewer')
library(VennDiagram)
library(RColorBrewer)

# Read data
mutation_has_source <- read.csv("../db_tables/mutation_has_source_new.tsv", header = TRUE, sep= '\t')
source = read.csv("../db_tables/source.tsv", header = TRUE, sep='\t')

mutation_has_source <- merge(mutation_has_source, source, by="id_source", all.mutation_has_source = TRUE)

# No. unique mutations
#length(unique(mutation_has_source$id_mutation)) # 1_660_059 OK

# Debo generar 4 sets con las mutaciones en cada db
clinvar <- subset(mutation_has_source, source == "clinvar")$id_mutation
cosmic <- subset(mutation_has_source, source == "cosmic")$id_mutation
disgenet <- subset(mutation_has_source, source == "disgenet")$id_mutation
uniprot <- subset(mutation_has_source, source == "uniprot")$id_mutation

# Paleta
myCol <- brewer.pal(4, "Pastel2")

# Plot
plt <- venn.diagram(
  x = list(clinvar, cosmic, disgenet, uniprot),
  category.names = c("ClinVar", "COSMIC", "DisGeNET", "UniProt"),
  filename = "venn_mutation.png",
  output = TRUE,
  
  # Output features
  imagetype = 'png',
  #height = 500,
  #width= 500,
  #resolution= 300,
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  
  # Numbers
  cex = 1.2,
  #fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 1.2,
  cat.fontface = "bold",
  #cat.default.pos = "outer",
  #cat.pos = c(-27, 27, 135),
  #cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  #rotation = 1
)
