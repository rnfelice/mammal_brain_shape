# this script grafts together the Alvarez-Carretero et al 2021 tree and the Faurby and Svenning 2015 tree,
# to create a more complete tree that includes fossils for use in subsequent analyses. The grafting is done using the tree.merger function in the RRphylo package,
# which grafts taxa onto a backbone tree based on their sister relationships.

# load packages
library(ape)
library(geiger)
library(paleomorph)
library(tidyverse)
library(phytools)
library(diversitree)
library(scales)
library(RRphylo)

# The grafted tree is then saved as a nexus file for use in subsequent analyses.
# first load in the Faurby and Svenning tree
mamTree_FS <- read.tree("./phylogeny_construction/input_trees/Complete_phylogeny_Faurby_Svenning.txt")


# rename tree taxa, where needed
mamTree_FS$tip.label <- replace(mamTree_FS$tip.label, mamTree_FS$tip.label == "Abrothrix_olivaceus", "Abrothrix_olivacea")
mamTree_FS$tip.label <- replace(mamTree_FS$tip.label, mamTree_FS$tip.label == "Archaeolemur_edwardsi", "Archaeolemur_sp")
mamTree_FS$tip.label <- replace(mamTree_FS$tip.label, mamTree_FS$tip.label == "Canis_dirus", "Aenocyon_dirus")
mamTree_FS$tip.label <- replace(mamTree_FS$tip.label, mamTree_FS$tip.label == "Cavia_tschudii", "Cavia_porcellus")
mamTree_FS$tip.label <- replace(mamTree_FS$tip.label, mamTree_FS$tip.label == "Catonyx_cuvieri", "Catonyx_tarijensis")
mamTree_FS$tip.label <- replace(mamTree_FS$tip.label, mamTree_FS$tip.label == "Myodes_glareolus", "Clethrionomys_glareolus")
mamTree_FS$tip.label <- replace(mamTree_FS$tip.label, mamTree_FS$tip.label == "Ctenodactylus_gundi", "Ctenodactylus_sp")
mamTree_FS$tip.label <- replace(mamTree_FS$tip.label, mamTree_FS$tip.label == "Ellobius_alaicus", "Ellobius_sp")
mamTree_FS$tip.label <- replace(mamTree_FS$tip.label, mamTree_FS$tip.label == "Felis_silvestris", "Felis_lybica")
mamTree_FS$tip.label <- replace(mamTree_FS$tip.label, mamTree_FS$tip.label == "Glyptodon_clavipes", "Glyptodon_sp")
mamTree_FS$tip.label <- replace(mamTree_FS$tip.label, mamTree_FS$tip.label == "Hystrix_cristata", "Hystrix_indica")
mamTree_FS$tip.label <- replace(mamTree_FS$tip.label, mamTree_FS$tip.label == "Lonchophylla_thomasi", "Hsunycteris_thomasi")
mamTree_FS$tip.label <- replace(mamTree_FS$tip.label, mamTree_FS$tip.label == "Nothrotheriops_shastense", "Nothrotheriops_shastensis")
mamTree_FS$tip.label <- replace(mamTree_FS$tip.label, mamTree_FS$tip.label == "Blanfordimys_afghanus", "Microtus_yuldaschi")
mamTree_FS$tip.label <- replace(mamTree_FS$tip.label, mamTree_FS$tip.label == "Tadarida_thersites", "Mops_thersites")
mamTree_FS$tip.label <- replace(mamTree_FS$tip.label, mamTree_FS$tip.label == "Monachus_tropicalis", "Neomonachus_tropicalis")
mamTree_FS$tip.label <- replace(mamTree_FS$tip.label, mamTree_FS$tip.label == "Mylodon_darwinii", "Mylodon_sp")
mamTree_FS$tip.label <- replace(mamTree_FS$tip.label, mamTree_FS$tip.label == "Palaeolama_weddeli", "Palaeolama_sp")
mamTree_FS$tip.label <- replace(mamTree_FS$tip.label, mamTree_FS$tip.label == "Speothos_venaticus", "Speothos_pacivorus")
mamTree_FS$tip.label <- replace(mamTree_FS$tip.label, mamTree_FS$tip.label == "Martes_pennanti", "Pekania_pennanti")
mamTree_FS$tip.label <- replace(mamTree_FS$tip.label, mamTree_FS$tip.label == "Sicista_armenica", "Sicista_concolor")
mamTree_FS$tip.label <- replace(mamTree_FS$tip.label, mamTree_FS$tip.label == "Tamandua_tetradactyla", "Tamandua_sp")
mamTree_FS$tip.label <- replace(mamTree_FS$tip.label, mamTree_FS$tip.label == "Tarsius_bancanus", "Tarsius_tarsier")
# mamTree_FS$tip.label <-replace(mamTree_FS$tip.label, mamTree_FS$tip.label == "Toxodon_platensis", "Nesodon_imbricatus")
mamTree_FS$tip.label <- replace(mamTree_FS$tip.label, mamTree_FS$tip.label == "Trigonodops_lopesi", "Adinotherium_ovinum")
mamTree_FS$tip.label <- replace(mamTree_FS$tip.label, mamTree_FS$tip.label == "Alopex_lagopus", "Vulpes_lagopus")
mamTree_FS$tip.label <- replace(mamTree_FS$tip.label, mamTree_FS$tip.label == "Maokopia_ronaldi", "Neohelos_stirtoni")
mamTree_FS$tip.label <- replace(mamTree_FS$tip.label, mamTree_FS$tip.label == "Otolemur_garnetti", "Otolemur_garnettii")

# create list of species missing from tree
missing_taxa_FS <- setdiff(species, mamTree_FS$tip.label)

# create list of taxa on tree minus missing taxa
species_FS_tree <- species[!species %in% missing_taxa_FS]

# modify tree to match dataset
mamTree_B <- keep.tip(mamTree_FS, species_FS_tree)

# ultrametricise tree to account for rounding errors during pruning
mamTree_B <- force.ultrametric(mamTree_B, method = "extend")

# plot(mamTree_B, type = "fan", cex = 0.2)

# Find missing taxa in Tree B
which(is.na(match(species, mamTree_B$tip.label)))
missing_taxa_B <- setdiff(species, mamTree_B$tip.label)

# Find taxa in Tree B that are not in Tree A
diff_AB <- setdiff(missing_taxa_A, missing_taxa_B)

length(diff_AB)


sisters_B <- lapply(c(1:length(diff_AB)), function(x) phytools::getSisters(mamTree_B, node = diff_AB[x], mode = "label"))
sisters_B <- as.matrix(unlist(sisters_B))
pairs_B <- cbind(diff_AB, sisters_B[, 1])
rownames(pairs_B) <- NULL

colnames(pairs_B) <- c("bind", "reference")
pairs_B <- as_tibble(pairs_B)
pairs_B$poly <- "FALSE" # add column specifying if graft is poltomy (always FALSE for our purposes)

View(pairs_B)
#### tools for finding position of taxa on tree
# find sister clade to missing taxa
plot(extract.clade(mamTree_B, 561), cex = 0.5)
# find most recent common ancestor
findMRCA(full_mamTree, c("Megaladapis_madagascariensis", "Homo_sapiens"))

# Add reference clade names
pairs_B[1, ]$reference <- "Equus_quagga-Tapirus_terrestris" # Adinotherium
pairs_B[2, ]$reference <- "Lupulella_mesomelas-Canis_lupus" # Aenocyon
pairs_B[3, ]$reference <- "Abrothrix_olivacea-Aepeomys_lugens" # Anotomys
pairs_B[4, ]$reference <- "Babakotia_radofilai-Propithecus_verreauxi" # Archaeolemur_sp
pairs_B[5, ]$reference <- "Tremarctos_ornatus" # 	Arctodus_simus
pairs_B[6, ]$reference <- "Propithecus_verreauxi-Indri_indri" # Babakotia
pairs_B[7, ]$reference <- "Choloepus_hoffmanni" # Catonyx
pairs_B[8, ]$reference <- "Isoodon_obesulus-Peroryctes_broadbenti" # Chaeropus
pairs_B[9, ]$reference <- "Lasiorhinus_latifrons-Vombatus_ursinus" # Diprotodon
pairs_B[10, ]$reference <- "Catonyx_tarijensis" # Glossotherium_robustum
pairs_B[11, ]$reference <- "Chlamyphorus_truncatus-Cabassous_unicinctus" # Glyptodon_sp
pairs_B[12, ]$reference <- "Panthera_uncia-Felis_lybica" # Homotherium
pairs_B[13, ]$reference <- "Bradypus_tridactylus" # Megalonyx
pairs_B[14, ]$reference <- "Bradypus_tridactylus-Megalonyx_jeffersonii" # Megatherium
pairs_B[15, ]$reference <- "Glossotherium_robustum" # Mylodon
pairs_B[16, ]$reference <- "Diprotodon_optatum-Zygomaturus_trilobus" # Neohelos
pairs_B[17, ]$reference <- "Megalonyx_jeffersonii" # Nothrotheriops
pairs_B[18, ]$reference <- "Lama_glama-Lama_guanicoe" # Palaeolama_sp
pairs_B[19, ]$reference <- "Mylodon_sp-Glossotherium_robustum" # Paramylodon
pairs_B[20, ]$reference <- "Tayassu_pecari-Pecari_tajacu" # Platygonus_compressus
pairs_B[21, ]$reference <- "Dorcopsis_veterum-Setonix_brachyurus" # Simosthenurus
pairs_B[23, ]$reference <- "Myrmecobius_fasciatus-Dasyurus_hallucatus" # Thylacinus
pairs_B[24, ]$reference <- "Phascolarctos_cinereus-Lasiorhinus_latifrons" # Thylacoleo
pairs_B[25, ]$reference <- "Molossus_molossus-Tadarida_brasiliensis" # Tomopeas
pairs_B[26, ]$reference <- "Anomalurus_beecrofti-Idiurus_macrotis" # Zenkerella_insignis
pairs_B[27, ]$reference <- "Diprotodon_optatum" # Zygomaturus_trilobus

# convert to data frame
pairs_B2 <- data.frame(pairs_B)

# check list
setdiff(species, classifier$species)


# create list of tip ages of taxa to be grafted
tip.ages_B <- c(classifier$tip_age[(which(classifier$species %in% pairs_B2$bind))])
names(tip.ages_B) <- pairs_B2$bind # add taxon names to tip age list


# Merge new taxa with full tree
merged_tree_B <- tree.merger(
    backbone = full_mamTree, data = pairs_B2,
    tip.ages = tip.ages_B,
    plot = F
)

plot(merged_tree_B, type = "fan", cex = 0.35)
