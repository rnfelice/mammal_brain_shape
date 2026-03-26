library(RRphylo)
library(paleotree)
library(readr)
library(phytools)
# import species data
classifier <-read.csv("./input_data/mammal_data.csv", fileEncoding="latin1")

# get list of species in dataset
species <- classifier$species

#read the tip-dated full tree
#this is the maximum clade credibility tree with mean node ages
#in other words, we took the output of BEAST and summarized it to a single tree
#in tree annotator
tree_A1C1F1 <- read.nexus(file = './phylogeny_construction/BEAST2/supertree_mcc.trees')
#trim down to only the 466 tips for which we have shape data
tree_A1C1F1_final <-keep.tip(tree_A1C1F1, species)
root_age<-max(dateNodes(tree_A1C1F1))

clade_data <- read_csv("./phylogeny_construction/alternate_topologies.csv")

tips_clade_a <- strsplit(clade_data$clade_definition[1], split="-")[[1]]
tips_target_a <- strsplit(clade_data$clade_pos_2[1], split="-")[[1]]


tree_A2C1F1 <- move.lineage(tree_A1C1F1,
focal=getMRCA(tree_A1C1F1, tips_clade_a),
sister=getMRCA(tree_A1C1F1, tips_target_a),
poly=FALSE,rescale = FALSE,rootage=root_age)
tree_A2C1F1_final <-keep.tip(tree_A2C1F1, species)


cairo_pdf("./phylogeny_construction/diagnostic_figs/tree_A2C1F1(hyp2).pdf",width=,height=7)
plot(extract.clade(tree_A2C1F1_final,findMRCA(tree_A2C1F1_final,c("Uintatherium_anceps","Equus_quagga"))), cex = 0.6, lwd = 2)
title("Tree A2C1F1")
dev.off()

tips_clade_c <- strsplit(clade_data$clade_definition[2], split="-")[[1]]
tips_target_c <- strsplit(clade_data$clade_pos_2[2], split="-")[[1]]




tree_A2C2F1 <- move.lineage(tree_A2C1F1,
                            focal=getMRCA(tree_A2C1F1, tips_clade_c),
                            sister=getMRCA(tree_A2C1F1, tips_target_c),
                            poly=FALSE,rescale = FALSE,rootage=root_age)
tree_A2C2F1_final <-keep.tip(tree_A2C2F1, species)


cairo_pdf("./phylogeny_construction/diagnostic_figs/tree_A2C2F1(hyp3).pdf",width=,height=7)
plot(extract.clade(tree_A2C2F1_final,findMRCA(tree_A2C2F1_final,c("Uintatherium_anceps","Equus_quagga"))), cex = 0.6, lwd = 2)
title("Tree A2C2F1")
dev.off()




tips_clade_f <- strsplit(clade_data$clade_definition[3], split="-")[[1]]
tips_target_f <- strsplit(clade_data$clade_pos_2[3], split="-")[[1]]




tree_A2C2F2 <- move.lineage(tree_A2C2F1,
                            focal=getMRCA(tree_A2C2F1, tips_clade_f),
                            sister=getMRCA(tree_A2C2F1, tips_target_f),
                            poly=FALSE,rescale = FALSE,rootage=root_age)
tree_A2C2F2_final <-keep.tip(tree_A2C2F2, species)


cairo_pdf("./phylogeny_construction/diagnostic_figs/tree_A2C2F2(hyp4).pdf",width=,height=7)
plot(extract.clade(tree_A2C2F2_final,findMRCA(tree_A2C2F2_final,c("Uintatherium_anceps","Equus_quagga"))), cex = 0.6, lwd = 2)
title("Tree A2C2F2")
dev.off()


tree_A1C2F1 <- move.lineage(tree_A1C1F1,
                            focal=getMRCA(tree_A1C1F1, tips_clade_c),
                            sister=getMRCA(tree_A1C1F1, tips_target_c),
                            poly=FALSE,rescale = FALSE,rootage=root_age)
tree_A1C2F1_final <-keep.tip(tree_A1C2F1, species)


cairo_pdf("./phylogeny_construction/diagnostic_figs/tree_A1C2F1(hyp7).pdf",width=,height=7)
plot(extract.clade(tree_A1C2F1_final,findMRCA(tree_A1C2F1_final,c("Uintatherium_anceps","Equus_quagga"))), cex = 0.6, lwd = 2)
title("Tree A1C2F1")
dev.off()


tree_A1C1F2 <- move.lineage(tree_A1C1F1,
                            focal=getMRCA(tree_A1C1F1, tips_clade_f),
                            sister=getMRCA(tree_A1C1F1, tips_target_f),
                            poly=FALSE,rescale = FALSE,rootage=root_age)
tree_A1C1F2_final <-keep.tip(tree_A1C1F2, species)


cairo_pdf("./phylogeny_construction/diagnostic_figs/tree_A1C1F2(hyp6).pdf",width=,height=7)
plot(extract.clade(tree_A1C1F2_final,findMRCA(tree_A1C1F2_final,c("Uintatherium_anceps","Equus_quagga"))), cex = 0.6, lwd = 2)
title("Tree A1C1F2")
dev.off()


tree_A1C2F2 <- move.lineage(tree_A1C2F1,
                            focal=getMRCA(tree_A1C2F1, tips_clade_f),
                            sister=getMRCA(tree_A1C2F1, tips_target_f),
                            poly=FALSE,rescale = FALSE,rootage=root_age)
tree_A1C2F2_final <-keep.tip(tree_A1C2F2, species)


cairo_pdf("./phylogeny_construction/diagnostic_figs/tree_A1C2F2(hyp5).pdf",width=,height=7)
plot(extract.clade(tree_A1C2F2_final,findMRCA(tree_A1C2F2_final,c("Uintatherium_anceps","Equus_quagga"))), cex = 0.6, lwd = 2)
title("Tree A1C2F2")
dev.off()

tree_A2C1F2 <- move.lineage(tree_A2C1F1,
                            focal=getMRCA(tree_A2C1F1, tips_clade_f),
                            sister=getMRCA(tree_A2C1F1, tips_target_f),
                            poly=FALSE,rescale = FALSE,rootage=root_age)
tree_A2C1F2_final <-keep.tip(tree_A2C1F2, species)


cairo_pdf("./phylogeny_construction/diagnostic_figs/tree_A2C1F2(hyp8).pdf",width=,height=7)
plot(extract.clade(tree_A2C1F2_final,findMRCA(tree_A2C1F2_final,c("Uintatherium_anceps","Equus_quagga"))), cex = 0.6, lwd = 2)
title("Tree A2C1F2")
dev.off()

cairo_pdf("./phylogeny_construction/diagnostic_figs/tree_A1C1F1(hyp1).pdf",width=,height=7)
plot(extract.clade(tree_A1C1F1_final,findMRCA(tree_A1C1F1_final,c("Uintatherium_anceps","Equus_quagga"))), cex = 0.6, lwd = 2)
title("Tree A1C1F1")
dev.off()


#check
edge_list<-list(
tree1 = tree_A1C1F1_final$edge,
tree2 = tree_A2C1F1_final$edge,
tree3 = tree_A2C2F1_final$edge,
tree4 = tree_A2C2F2_final$edge,
tree5 = tree_A1C2F2_final$edge,
tree6 = tree_A1C1F2_final$edge,
tree7 = tree_A1C2F1_final$edge,
tree8 = tree_A2C1F2_final$edge)

combos <- t(as.data.frame(combn(names(edge_list), 2)))

for(i in 1:nrow(combos)){
 if(!identical(edge_list[combos[i,1]],edge_list[combos[i,2]])){
   print(paste0(combos[i,1], " and ", combos[i,2], " are different"))
 } else {
   print(paste0(combos[i,1], " and ", combos[i,2], " are identical"))
  
}
}
write.nexus(tree_A1C1F1_final, file = "./alternative_topologies/hyp1.nex")
write.nexus(tree_A2C1F1_final, file = "./alternative_topologies/hyp2.nex")
write.nexus(tree_A2C2F1_final, file = "./alternative_topologies/hyp3.nex")
write.nexus(tree_A2C2F2_final, file = "./alternative_topologies/hyp4.nex")
write.nexus(tree_A1C2F2_final, file = "./alternative_topologies/hyp5.nex")
write.nexus(tree_A1C1F2_final, file = "./alternative_topologies/hyp6.nex")
write.nexus(tree_A1C2F1_final, file = "./alternative_topologies/hyp7.nex")
write.nexus(tree_A2C1F2_final, file = "./alternative_topologies/hyp8.nex")

