#this script grafts individual fossil taxa onto the full tree created in the previous step
#using the tree.merger function in the RRphylo()

# load packages
library(ape)
library(geiger)
library(paleomorph)
library(tidyverse)
library(RRphylo)

# create list of species missing from tree
missing_taxa_D <- setdiff(species, merged_tree_C$tip.label)


# create list of taxa on tree minus missing taxa
species_D <- species[ ! species %in% missing_taxa_D ]



# find sister clade to missing taxa
plot(extract.clade(mamTree_scaled_prefossil,437), cex = 0.6, type = "fan")

findMRCA(mamTree_scaled_prefossil,c("Chlamyphorus_truncatus","Mus_musculus"))

# Redo fossil species graft
# Manually attach fossil taxa to tree
# pairs_D <- data.frame(bind=c(setdiff(species,merged_tree_C$tip.label)),
pairs_D <- data.frame(bind=c(setdiff(species,merged_tree_C$tip.label)),
                      reference=c("Loris_tardigradus-Lemur_catta", # Adapis_parisiensis
                                  "Cephalophus_leucogaster-Orcaella_heinsohni", # Agriochoerus_sp
                                  "Caperea_marginata", # Albacetus_salvifactus
                                  "Ceratotherium_simum-Rhinoceros_unicornis", # Amynodon_advenus
                                  "Agriochoerus_sp-Cephalophus_leucogaster", # Anoplotherium_sp
                                  "Quercygale_angustidens-Lynx_rufus", # Apterodon_macrognathus ####
                                  "Hippopotamus_amphibius-Orcaella_heinsohni", # Archaeotherium_mortoni
                                  "Cebochoerus_lacustris-Cephalophus_leucogaster", # Arctocyon_primaevus
                                  "Delphinus_delphis-Physeter_catodon", # Argyrocetus_joaquinensis
                                  "Physeter_catodon", # Aulophyseter_morricei
                                  "Homo_sapiens", # Australopithecus_africanus
                                  "Cephalophus_leucogaster-Equus_quagga", # Barylambda_schmidti
                                  "Anoplotherium_sp", # Caenomeryx_filholi
                                  "Leptocyon_sp-Lycaon_pictus", # Carpocyon_webbi
                                  "Mouillacitherium_elegans-Cephalophus_leucogaster", # Cebochoerus_lacustris
                                  "Ratufa_affinis-Sciurus_vulgaris", # Cedromus_wilsoni
                                  "Manis_javanica-Patriomanis_americana", # Chadronia_margaretae ###########
                                  "Equus_quagga-Equus_caballus", # Cormohipparion_occidentale
                                  "Lycaon_pictus-Hesperocyon_gregarius", # Cynodictis_cayluxi
                                  "Apterodon_macrognathus", # Cynohyaenodon_cayluxi
                                  "Platygonus_compressus-Tayassu_pecari", # Desmathyus_sp
                                  "Lama_glama-Cephalophus_leucogaster", # Dichobune_leporina
                                  "Odobenus_rosmarus-Neomonachus_tropicalis", # Enaliarctos_sp
                                  "Nandinia_binotata-Procyon_lotor", # Eusmilus_bidentatus
                                  "Castor_canadensis-Castor_fiber", # Eutypomys_thompsoni
                                  "Cynodictis_cayluxi", # Gustafsonia_sp
                                  "Dugong_dugon", # Halitherium_schinzi
                                  "Megalonyx_jeffersonii", # Hapalops_sp
                                  "Tapirus_terrestris-Tapirus_indicus", # Heptodon_sp
                                  "Lycaon_pictus-Carpocyon_webbi", # Hesperocyon_gregarius
                                  "Leontinia_gaudryi-Nesodon_imbricatus", # Homalodotherium_cunninghami
                                  "Amynodon_advenus-Ceratotherium_simum", # Hyrachyus_modestus
                                  "Mesohippus_bairdi-Equus_quagga", # Hyracotherium_sp
                                  "Loris_tardigradus-Homo_sapiens", # Ignacius_graybullianus
                                  "Aplodontia_rufa-Sciurus_vulgaris", # Ischyromys_typus
                                  "Orcaella_heinsohni-Protocetus_atavus", # Khirtharia_inflata
                                  "Acrobates_pygmaeus-Mus_musculus", # Kryptobaatar_dashzevegi
                                  "Adinotherium_ovinum-Nesodon_imbricatus", # Leontinia_gaudryi
                                  "Agriochoerus_sp", # Leptauchenia_sp
                                  "Chlamyphorus_truncatus-Mus_musculus", # Leptictis_sp
                                  "Lycaon_pictus-Urocyon_littoralis", # Leptocyon_sp
                                  "Ceratotherium_simum-Tapirus_terrestris", # Megacerops_coloradensis
                                  "Lepus_europaeus-Oryctolagus_cuniculus", # Megalagus_turgidus
                                  "Equus_quagga-Cormohipparion_occidentale", # Merychippus_severus
                                  "Leptauchenia_sp", # Merycochoerus_proprius
                                  "Megacerops_coloradensis", # Mesatirhinus_junius
                                  "Aplodontia_rufa", # Mesogaulus_paniensis
                                  "Equus_quagga-Merychippus_severus", # Mesohippus_bairdi
                                  "Tarsius_tarsier", # Microchoerus_erinaceus
                                  "Homo_sapiens-Ignacius_graybullianus", # Microsyops_annectens
                                  "Anoplotherium_sp-Caenomeryx_filholi", # Mixtotherium_sp
                                  "Caenomeryx_filholi-Dichobune_leporina", # Mouillacitherium_elegans
                                  "Microchoerus_erinaceus", # Necrolemur_antiquus
                                  "Adinotherium_ovinum", # Nesodon
                                  "Adapis_parisiensis", # Notharctus_tenebrosus
                                  "Ornithorhynchus_anatinus", # Obdurodon_dicksoni
                                  "Lama_glama-Palaeolama_sp", # Oxydactylus_sp
                                  "Megacerops_coloradensis-Mesatirhinus_junius", # Palaeosyops_sp
                                  "Carpocyon_webbi", # Paracynarctus_sinclairi
                                  "Rapamys_atramontis-Castor_canadensis", # Paramys_delicatus #######
                                  "Ateles_geoffroyi-Homo_sapiens", # Parapithecus_grangeri
                                  "Manis_javanica-Smutsia_gigantea", # Patriomanis_americana
                                  "Equus_quagga-Adinotherium_ovinum", # Phenacodus_primaevus
                                  "Pekania_pennanti-Martes_martes", # Plionictis_sp
                                  "Mellivora_capensis-Enhydra_lutris", # Promartes_olcotti
                                  "Catonyx_tarijensis", # Proscelidodon_patrius
                                  "Myocastor_coypus-Echimys_chrysurus", # Prospaniomys_priscus
                                  "Homalodotherium_cunninghami-Protypotherium_australe", # Proterotherium
                                  "Caperea_marginata-Orcaella_heinsohni", # Protocetus_atavus
                                  "Sciurus_vulgaris-Tamiasciurus_hudsonicus", # Protosciurus_rachelae
                                  "Nesodon_imbricatus-Homalodotherium_cunninghami", # Protypotherium_australe
                                  "Smilodon_fatalis-Homotherium_serum", # Pseudaelurus_validus
                                  "Cormohipparion_occidentale", # Pseudhipparion_gratum
                                  "Paramys_delicatus", # Pseudotomus_horribilis
                                  "Protypotherium_australe", # Pseudotypotherium_pseudopachygnathum
                                  "Eusmilus_bidentatus-Urocyon_littoralis", # Quercygale_angustidens ###
                                  "Glis_glis-Castor_canadensis", # Rapamys_atramontis
                                  "Microchoerus_erinaceus-Necrolemur_antiquus", # Rooneyia_viejaensis
                                  "Notharctus_tenebrosus", # Smilodectes_gracilis
                                  "Megatherium_americanum", # Thalassocnus_natans
                                  "Barylambda_schmidti", # Titanoides_primaevus
                                  "Barylambda_schmidti-Titanoides_primaevus", # Trogosus_hillsii ######
                                  "Pseudotypotherium_pseudopachygnathum", # Typotheriopsis_sp
                                  "Trogosus_hillsii-Alcelaphus_buselaphus", # Uintatherium #####################
                                  "Merycochoerus_proprius", # Ustatochoerus ############
                                  "Promartes_olcotti" # Zodiolestes_daimonedlixensis
                      ))




# add column specifying if graft is poltomy (always FALSE for our purposes)
pairs_D$poly <- "FALSE" 

# convert to data frame
pairs_D2 <- data.frame(pairs_D)
View(pairs_D2)
# create list of tip ages of taxa to be grafted
tip.ages_D <- c(classifier$tip_age[(which(classifier$species %in% pairs_D2$bind))])
names(tip.ages_D) = pairs_D2$bind # add taxon names to tip age list
# View(classifier)

# Graft fossil taxa to full tree
merged_tree_D <- tree.merger(backbone=merged_tree_C,data=pairs_D2,
                           tip.ages=tip.ages_D,
                           plot=F)

# merged_tree_D <- tree.merger(backbone=mamTree_scaled_prefossil,data=pairs_D2,
                             # tip.ages=tip.ages_D, 
                             # plot=TRUE)



#write the merged tree with ALL extants plus our fossils
write.nexus(merged_tree_D , file = './input_data/merged_tree.nex')