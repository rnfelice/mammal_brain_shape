# this script imports the time-scaled mammal tree from Alvarez-Carretero et al 2021
# edits tip labels to match species names in dataset (for example adding species names or removing subspecies names)
# The edited tree is saved as a nexus file for use in subsequent analyses.

# load packages
library(ape)
library(geiger)
library(paleomorph)
library(tidyverse)
library(phytools)
library(diversitree)
library(scales)
library(RRphylo)


# import time-scaled mammal tree (Alvarez-Carretero et al 2021)
full_mamTree <- read.nexus("./phylogeny_construction/input_trees/4705sp_mammal-time.txt")

# ultrametricise tree to account for rounding errors during pruning
full_mamTree <- force.ultrametric(full_mamTree, method = "extend")

# Scale tree to units of 1 Ma
full_mamTree$edge.length <- (full_mamTree$edge.length) * 100

# capitalise generic name of each species
full_mamTree$tip.label <- gsub("\\b([A-Za-z])", "\\U\\1", full_mamTree$tip.label, perl = TRUE)


# rename tree taxa, where needed
full_mamTree$tip.label <- replace(full_mamTree$tip.label, full_mamTree$tip.label == "Abrothrix_olivaceus", "Abrothrix_olivacea")
full_mamTree$tip.label <- replace(full_mamTree$tip.label, full_mamTree$tip.label == "Anomalurus", "Anomalurus_beecrofti")
full_mamTree$tip.label <- replace(full_mamTree$tip.label, full_mamTree$tip.label == "Bos_bison", "Bison_bison")
full_mamTree$tip.label <- replace(full_mamTree$tip.label, full_mamTree$tip.label == "Myodes_glareolus", "Clethrionomys_glareolus")
full_mamTree$tip.label <- replace(full_mamTree$tip.label, full_mamTree$tip.label == "Ctenodactylus_gundi", "Ctenodactylus_sp")
full_mamTree$tip.label <- replace(full_mamTree$tip.label, full_mamTree$tip.label == "Ellobius_talpinus", "Ellobius_sp")
full_mamTree$tip.label <- replace(full_mamTree$tip.label, full_mamTree$tip.label == "Equus_grevyi", "Equus_quagga")
full_mamTree$tip.label <- replace(full_mamTree$tip.label, full_mamTree$tip.label == "Cercopithecus_ascanius_ascanius", "Cercopithecus_ascanius")
full_mamTree$tip.label <- replace(full_mamTree$tip.label, full_mamTree$tip.label == "Rheomys_raptor", "Chibchanomys_trichotis")
full_mamTree$tip.label <- replace(full_mamTree$tip.label, full_mamTree$tip.label == "Eptesicus_bottae", "Eptesicus_andinus")
full_mamTree$tip.label <- replace(full_mamTree$tip.label, full_mamTree$tip.label == "Eulemur_macaco_macaco", "Eulemur_macaco")
full_mamTree$tip.label <- replace(full_mamTree$tip.label, full_mamTree$tip.label == "Graphiurus_murinus", "Graphiurus_nagtglasii")
full_mamTree$tip.label <- replace(full_mamTree$tip.label, full_mamTree$tip.label == "Hipposideros_turpis_turpis", "Hipposideros_turpis")
full_mamTree$tip.label <- replace(full_mamTree$tip.label, full_mamTree$tip.label == "Lophocebus_aterrimus", "Lophocebus_albigena")
full_mamTree$tip.label <- replace(full_mamTree$tip.label, full_mamTree$tip.label == "Macrotarsomys", "Macrotarsomys_bastardi")
full_mamTree$tip.label <- replace(full_mamTree$tip.label, full_mamTree$tip.label == "Megaladapis_edwardsi", "Megaladapis_madagascariensis")
full_mamTree$tip.label <- replace(full_mamTree$tip.label, full_mamTree$tip.label == "Microgale_majori", "Microgale_drouhardi")
full_mamTree$tip.label <- replace(full_mamTree$tip.label, full_mamTree$tip.label == "Myospalax_aspalax", "Myospalax_myospalax")
full_mamTree$tip.label <- replace(full_mamTree$tip.label, full_mamTree$tip.label == "Nyctimene_robinsoni", "Nyctimene_rabori")
full_mamTree$tip.label <- replace(full_mamTree$tip.label, full_mamTree$tip.label == "Petaurus_breviceps", "Petaurus_australis")
full_mamTree$tip.label <- replace(full_mamTree$tip.label, full_mamTree$tip.label == "Pipistrellus_pipistrellus", "Pipistrellus_hesperus")
full_mamTree$tip.label <- replace(full_mamTree$tip.label, full_mamTree$tip.label == "Typhlomys_cinereus", "Platacanthomys_lasiurus")
full_mamTree$tip.label <- replace(full_mamTree$tip.label, full_mamTree$tip.label == "Rhipidomys_nitela", "Rhipidomys_fulviventer")
full_mamTree$tip.label <- replace(full_mamTree$tip.label, full_mamTree$tip.label == "Manis_tricuspis", "Smutsia_gigantea")
full_mamTree$tip.label <- replace(full_mamTree$tip.label, full_mamTree$tip.label == "Sylvilagus_aquaticus", "Sylvilagus_bachmani")
full_mamTree$tip.label <- replace(full_mamTree$tip.label, full_mamTree$tip.label == "Trichechus_manatus", "Trichechus_senegalensis")
full_mamTree$tip.label <- replace(full_mamTree$tip.label, full_mamTree$tip.label == "Varecia_variegata_variegata", "Varecia_variegata")
full_mamTree$tip.label <- replace(full_mamTree$tip.label, full_mamTree$tip.label == "Zaglossus_bruijni", "Zaglossus_bartoni")
full_mamTree$tip.label <- replace(full_mamTree$tip.label, full_mamTree$tip.label == "Felis_silvestris", "Felis_lybica")
full_mamTree$tip.label <- replace(full_mamTree$tip.label, full_mamTree$tip.label == "Hystrix_cristata", "Hystrix_indica")
full_mamTree$tip.label <- replace(full_mamTree$tip.label, full_mamTree$tip.label == "Blanfordimys_afghanus", "Microtus_yuldaschi")
full_mamTree$tip.label <- replace(full_mamTree$tip.label, full_mamTree$tip.label == "Mops_bakarii", "Mops_thersites")
full_mamTree$tip.label <- replace(full_mamTree$tip.label, full_mamTree$tip.label == "Monachus_monachus", "Neomonachus_tropicalis")
full_mamTree$tip.label <- replace(full_mamTree$tip.label, full_mamTree$tip.label == "Speothos_venaticus", "Speothos_pacivorus")
full_mamTree$tip.label <- replace(full_mamTree$tip.label, full_mamTree$tip.label == "Steatomys_parvus", "Steatomys_pratensis")
full_mamTree$tip.label <- replace(full_mamTree$tip.label, full_mamTree$tip.label == "Tamandua_tetradactyla", "Tamandua_sp")
full_mamTree$tip.label <- replace(full_mamTree$tip.label, full_mamTree$tip.label == "Urocyon_littoralis_santarosae", "Urocyon_littoralis")
full_mamTree$tip.label <- replace(full_mamTree$tip.label, full_mamTree$tip.label == "Vandeleuria_sp_kcr_2008", "Vandeleuria_oleracea")

# import species data
classifier <- read.csv("./input_data/mammal_data.csv", fileEncoding = "latin1")

# get list of species in dataset
species <- classifier$species

# create list of species missing from tree
missing_taxa_A <- setdiff(species, full_mamTree$tip.label)

# create list of taxa on tree minus missing taxa
species_scaled <- species[!species %in% missing_taxa_A]

write.nexus(full_mamTree, file = "./input_data/mammaltree4705_editednames.nex")
