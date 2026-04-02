library(geomorph)
library(ape)

# load the data
# aligned specimens
n_landmarks <- 122
n_dimensions <- 3
raw_data <- read.csv("./input_data/procrutes_aligned_coords.csv", row.names = 1)
Y.gpa.endo <- arrayspecs(raw_data, p = n_landmarks, k = n_dimensions)

# centroid sizes
Csize <- read.csv("./input_data/centroid_sizes.csv", row.names = 1)
Csize.endo <- setNames((Csize[, 1]), rownames(Csize))

# body mass and taxonomic data
classifier <- readr::read_csv(file = paste("./input_data/mammal_data.csv", sep = ""))


# Import trees
tree_hyp1 <- read.nexus("./phylogeny_construction/final_tree_hyp1.nex")
tree_hyp2 <- read.nexus("./phylogeny_construction/final_tree_hyp2.nex")
tree_hyp3 <- read.nexus("./phylogeny_construction/final_tree_hyp3.nex")
tree_hyp4 <- read.nexus("./phylogeny_construction/final_tree_hyp4.nex")
tree_hyp5 <- read.nexus("./phylogeny_construction/final_tree_hyp5.nex")
tree_hyp6 <- read.nexus("./phylogeny_construction/final_tree_hyp6.nex")
tree_hyp7 <- read.nexus("./phylogeny_construction/final_tree_hyp7.nex")
tree_hyp8 <- read.nexus("./phylogeny_construction/final_tree_hyp8.nex")


# Check phylogenetic signal, evolutionary rates and disparity of landmark data


# Check phylogenetic signal for each tree hypothesis
physignal_hyp1 <- physignal(phy = tree_hyp1, Y.gpa.endo)
physignal_hyp2 <- physignal(phy = tree_hyp2, Y.gpa.endo)
physignal_hyp3 <- physignal(phy = tree_hyp3, Y.gpa.endo)
physignal_hyp4 <- physignal(phy = tree_hyp4, Y.gpa.endo)
physignal_hyp5 <- physignal(phy = tree_hyp5, Y.gpa.endo)
physignal_hyp6 <- physignal(phy = tree_hyp6, Y.gpa.endo)
physignal_hyp7 <- physignal(phy = tree_hyp7, Y.gpa.endo)
physignal_hyp8 <- physignal(phy = tree_hyp8, Y.gpa.endo)


# Compare evolutionary rates of different brain regions

# define brain regions
neocortex <- c(2:3, 10:21, 52:59, 74:113)
olfactory <- c(1, 4, 22:29, 114:122)
cerebellum <- c(5, 9, 30:32, 45:51, 60:68)
brainstem <- c(6:8, 33:44, 69:73)

rate.groups <- c(1:dim(Y.gpa.endo)[[1]])
rate.groups[1:dim(Y.gpa.endo)[[1]]] <- 0
rate.groups[neocortex] <- 1
rate.groups[cerebellum] <- 2
rate.groups[brainstem] <- 3
rate.groups[olfactory] <- 4


# compare between-region rates for all trees
rate.comparison_hyp1 <- compare.multi.evol.rates(Y.gpa.endo, gp = rate.groups, phy = tree_hyp1)
rate.comparison_hyp2 <- compare.multi.evol.rates(Y.gpa.endo, gp = rate.groups, phy = tree_hyp2)
rate.comparison_hyp3 <- compare.multi.evol.rates(Y.gpa.endo, gp = rate.groups, phy = tree_hyp3)
rate.comparison_hyp4 <- compare.multi.evol.rates(Y.gpa.endo, gp = rate.groups, phy = tree_hyp4)
rate.comparison_hyp5 <- compare.multi.evol.rates(Y.gpa.endo, gp = rate.groups, phy = tree_hyp5)
rate.comparison_hyp6 <- compare.multi.evol.rates(Y.gpa.endo, gp = rate.groups, phy = tree_hyp6)
rate.comparison_hyp7 <- compare.multi.evol.rates(Y.gpa.endo, gp = rate.groups, phy = tree_hyp7)
rate.comparison_hyp8 <- compare.multi.evol.rates(Y.gpa.endo, gp = rate.groups, phy = tree_hyp8)


rate.comparison$pairwise.pvalue
summary(rate.comparison)
plot(rate.comparison)


# Calculate morphological disparity between groups


# morphological disparity for families
disp.clade <- morphol.disparity(coords ~ 1,
    groups = ~superorder, data = gdf,
    iter = 999, print.progress = FALSE
)
summary(disp.clade)
disp.clade$Procrustes.var
morphol.disparity(f1 = mass_phy, groups = ~family, data = gdf, iter = 999)


# morphological disparity for families, adjusting for centroid size
disp.clade.Csize <- morphol.disparity(coords ~ Csize,
    groups = ~superorder, data = gdf,
    iter = 999, print.progress = FALSE
)
summary(disp.clade.Csize)
disp.clade.Csize$Procrustes.var


# morphological disparity for families, adjusting for centroid size
disp.clade.mass <- morphol.disparity(coords ~ mass,
    groups = ~superorder, data = gdf,
    iter = 999, print.progress = FALSE
)
summary(disp.clade.mass)
disp.clade.mass$Procrustes.var

# disparity on results of phylogenetic allometry analysis
disp.clade.allom <- morphol.disparity(f1 = Allometry.Csize_hyp1, groups = ~superorder, data = gdf, iter = 999)

summary(disp.clade.allom)
disp.clade.allom$Procrustes.var


# Disparity of different regions

# Neocortex
gdf.neocortex <- geomorph.data.frame(coords = Y.gpa.endo[neocortex, , ], mass = classifier$Mass, Csize = Csize.endo, superorder = classifier$superorder)
View(gdf.neocortex)
names(gdf.neocortex) <- c("coords", "mass", "Csize", "superorder")

# morphological disparity for families
disp.clade.neocortex.allom <- morphol.disparity(
    f1 = Allometry.neocortex.hyp1, groups = ~superorder, data = gdf.neocortex,
    iter = 999, print.progress = FALSE
)
summary(disp.clade.neocortex.allom)
disp.clade.neocortex$Procrustes.var


# Cerebellum

# morphological disparity for families
disp.clade.cerebellum.allom <- morphol.disparity(
    f1 = Allometry.cerebellum.hyp1, groups = ~superorder, data = gdf.cerebellum,
    iter = 999, print.progress = FALSE
)
summary(disp.clade.cerebellum.allom)
disp.clade.neocortex$Procrustes.var


# Brainstem

# morphological disparity for families
disp.clade.brainstem.allom <- morphol.disparity(
    f1 = Allometry.brainstem.hyp1, groups = ~superorder, data = gdf.brainstem,
    iter = 999, print.progress = FALSE
)
summary(disp.clade.brainstem.allom)
disp.clade.neocortex$Procrustes.var


# Olfactory

# morphological disparity for families
disp.clade.olfactory.allom <- morphol.disparity(
    f1 = Allometry.olfactory.hyp1, groups = ~superorder, data = gdf.olfactory,
    iter = 999, print.progress = FALSE
)
summary(disp.clade.olfactory.allom)


plot(log(gdf$Csize), log(classifier$Mass), bg = col.gp, pch = 21, cex = 1.5)
sort(classifier$Mass)
