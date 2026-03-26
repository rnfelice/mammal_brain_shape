#this script grafts the marsupial tree from Weisbecker and Yalkut 2020 onto the full tree created in the previous step

# load packages
library(ape)
library(geiger)
library(paleomorph)
library(tidyverse)
library(RRphylo)

# load marsupial tree
mamTree_mars <-read.nexus("./phylogeny_construction/input_trees/Weisbecker_Yalk_resolved.tre")
# find missing tips from original tree
mamTree_mars <- mamTree_mars$Tree1

which(is.na(match(species, mamTree_mars$tip.label)))

# list specimens in dataset missing from marsupial tree
missing_taxa_C <- setdiff(species, mamTree_mars$tip.label)

# list taxa in marsupial tree missing from original tree
missing_taxa_AB <- setdiff(species, merged_tree_B$tip.label)

# set difference
diff_AC <-setdiff(missing_taxa_AB, missing_taxa_C)


# make matrix of taxa to add to full tree
sisters_C <-lapply(c(1:length(diff_AC)), function(x) phytools::getSisters(mamTree_mars, node = diff_AC[x], mode="label"))
sisters_C <-as.matrix(unlist(sisters_C))
pairs_C<-cbind(diff_AC, sisters_C[,1])
rownames(pairs_C) = NULL

colnames(pairs_C)<-c("bind","reference")
pairs_C <-as_tibble(pairs_C)
pairs_C$poly <- "FALSE" # add column specifying if graft is poltomy (always FALSE for our purposes)
View(pairs_C)

# tools for finding sister clades
plot(extract.clade(merged_tree_B,4738), cex = 0.5)
findMRCA(merged_tree_B,c("Caenolestes_fuliginosus","Thylacinus_cynocephalus"))

# Add refence clade names
pairs_C[1,]$reference  <- "Wallabia_bicolor-Potorous_tridactylus"
pairs_C[2,]$reference  <- "Antechinus_flavipes-Planigale_ingrami"
pairs_C[3,]$reference  <- "Caenolestes_fuliginosus-Wallabia_bicolor" #Borhyaena_tuberata
pairs_C[5,]$reference  <- "Macrotis_lagotis-Isoodon_macrourus"
pairs_C[7,]$reference  <- "Neohelos_stirtoni" #Nimbadon_lavarackorum
pairs_C[8,]$reference  <- "Zygomaturus_trilobus-Nimbadon_lavarackorum"
pairs_C[9,]$reference  <- "Dromiciops_gliroides-Planigale_ingrami"

# convert to data frame
pairs_C2 <- data.frame(pairs_C)


# create list of tip ages of taxa to be grafted
tip.ages_C <- c(classifier$tip_age[(which(classifier$species %in% pairs_C2$bind))])
names(tip.ages_C) = pairs_C2$bind # add taxon names to tip age list


# merge taxa with full tree
merged_tree_C <- tree.merger(backbone=merged_tree_B,data=pairs_C2,
                           tip.ages=tip.ages_C, 
                           plot=F)

