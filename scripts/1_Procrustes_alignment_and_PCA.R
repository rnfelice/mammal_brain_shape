## ------------------------------------------------------------------------------------- ##
##         CONTENTS AND DESCRIPTION OF THE TASKS PERFORMED IN THIS SCRIPT                ##
## ------------------------------------------------------------------------------------- ##
# Script for reflecting and aligning mammal endocasts                                     #
#                                                                                         #
# Import dataset                                                                          #
# Symmetrise landmarks                                                                    #
# Procrustes alignment and PCA                                                            #
# Plot landmark variance of Procrustes-aligned data                                       #
#---------------------------------------------------------------#                         #
# QUESTIONS AND DOUBTS CAN BE ADDRESSED TO:                     #                         #
# ===============================================================#                         #
# Andrew Knapp <ucfaakn@ucl.ac.uk>                              #                         #
#---------------------------------------------------------------#                         #
#                                                                                         #
#--------------#                                                                          #
# LAST UPDATE  #                                                                          #
# ==============#                                                                          #
# 2026-02-16   #                                                                          #
#--------------#                                                                          #
#                                                                                         #
## ------------------------------------------------------------------------------------- ##


#----------------#
# LOAD LIBRARIES #
#----------------#
library(rlang)
library(geomorph)
library(Morpho)
library(readr)
library(ggplot2)
library(png)
library(ape)
library(phytools)
library(Rvcg)
library(vegan)
library(permute)
library(lattice)
library(paleomorph)
library(EMMLi)
library(ff)
library(dplyr)
library(geiger)
library(tidyr)
library(DescTools)
library(ggConvexHull)
library(ggpubr)
library(viridis)
library(ggrepel)
library(SURGE)
library(phytools)
library(viridis)
library(rphylopic)

#-----------------------#
# SOURCE FUNCTIONS FILE #
#-----------------------#

# path to Functions.R- contrains additional functions used for handling landmark configurations.

source(paste("./scripts/utilities/Functions.R", sep = ""))


# load in pre-resampled and slid landmarks

# load landmarks
n_landmarks <- 122
n_dimensions <- 3
raw_data <- read.csv("./input_data/raw_landmark_coords.csv", row.names = 1)
arr <- arrayspecs(raw_data, p = n_landmarks, k = n_dimensions)


# Designate brain regions (new landmark scheme)
neocortex <- c(2:3, 10:21, 52:59, 74:113)
olfactory <- c(1, 4, 22:29, 114:122)
cerebellum <- c(5, 9, 30:32, 45:51, 60:68)
brainstem <- c(6:8, 33:44, 69:73)

# check landmarks
spheres3d(arr[, , 1], r = 0.5)
spheres3d(arr[olfactory, , 1], r = 0.51, col = "dodgerblue")
spheres3d(arr[neocortex, , 1], r = 0.51, col = "red")
spheres3d(arr[cerebellum, , 1], r = 0.51, col = "gold")
spheres3d(arr[brainstem, , 1], r = 0.51, col = "chartreuse3")


#----------------------------------------#
# REFLECTION OF LANDMARKS ACROSS MIDLINE #
#----------------------------------------#
# 1. Define matrix for vectors to get plane to reflect
#    for each specimen. nrow = number of specimens

mat.vectors.plane <- matrix(rep(0, length(dimnames(arr)[[3]]) * 3), nrow = length(dimnames(arr)[[3]]), ncol = 3)
rownames(mat.vectors.plane) <- as.character(dimnames(arr)[[3]])
colnames(mat.vectors.plane) <- c("v1", "v2", "v3")

#
# 1.2. Append info. Modify as suitable "v.all.specimens".
#      You can have as many as you want and add them to
#      the corresponding row, i.e. specimen
#####  These three numbers are the 3 landmarks that will designate the skull midline #####
v.all.specimens <- c(1, 2, 7)
#


# apply midline landmarks to all specimens
for (i in seq(1:length(dimnames(arr)[[3]]))) {
  mat.vectors.plane[i, ] <- v.all.specimens
}

mat.vectors.plane


# Make sure Functions.R has been loaded - it contains the reflect.lmks script.

# 2. Get arguments for function "reflect.lmks" and run function

v <- mat.vectors.plane # When you change position of lmks for the specimens
#### Here you need to list ALL landmarks and semilandmarks that lie on the midline. ####
midline <- c(1, 2, 4, 5, 6, 7, 8, 10:17, 25:33, 38:39) # With semilandmarks: positions are in c( 1:7, 34:108 )
X <- arr # This is the imported landmarks
comb.dataset <- reflect.lmks(X = X, v = v, midline = midline, plot.res = T)


# Symmetrize
#
#  Identify original lmks (right side) that match
#      the reflected ones, omiting midline
all.lmks <- 1:length(1:nrow(arr))
lmks.nomidl <- all.lmks[-midline]
right <- lmks.nomidl
left <- (length(all.lmks) + 1):dim(comb.dataset$original)[1]
#
# Create matrix with left and right paired LM
pairedLM <- cbind(left, right)

# Create array for symmetric coordinates
sym.comb.dataset <- array(dim = c(
  lmks = dim(comb.dataset$original)[1],
  coords = dim(comb.dataset$original)[2],
  sp = dim(comb.dataset$original)[3]
))
dimnames(sym.comb.dataset) <- list(
  c(unlist(dimnames(comb.dataset$original)[1])),
  c("x", "y", "z"),
  c(unlist(dimnames(comb.dataset$original)[3]))
)


# Run Morpho::symmetrize for all specimens
#      and save results in the array
for (i in seq(1:length(dimnames(arr)[[3]]))) {
  sym.comb.dataset[, , i] <- Morpho::symmetrize(comb.dataset$original[, , i], pairedLM)
}


reflected_patched_endocasts <- comb.dataset$original

# save landmarks post-reflection
# save(reflected_patched_endocasts, file = "./input_data/reflected_patched_endocasts_mammals.rda")


#-----------------------------------------------------------#
# TEST/RUN geomorph::gpagen || CREATE PROCRUSTES ALIGNMENT  #
#-----------------------------------------------------------#

# Perform Procrustes, using only the LMs in the vector above
Proc.endo <- Morpho::procSym(comb.dataset$original, reflect = F)

# Save Procrustes-aligned landmark data
save(Proc.endo, file = "./input_data/Proc.endo.mammals.rda")

# remove reflected half
Proc.endo$rotated <- (Proc.endo$rotated[1:length(all.lmks), , ])

# Subset out the Procrustes coordinates
Y.gpa.endo <- Proc.endo$rotated


# extract endocast centroid sizes
Csize.endo <- Proc.endo$size

# 2.2. Plot original half after Procrustes - ALL specimens superimposed!
plotAllSpecimens(Y.gpa.endo)

# save aligned specimens
write.csv(two.d.array(Y.gpa.endo), file = "./input_data/procrutes_aligned_coords.csv", quote = F)
# save endocast centroid sizes
write.csv(as.matrix(Csize.endo), file = "./input_data/centroid_sizes.csv", quote = F)


###################################################################

# Read classifier data- contains information on taxonomy, ecology, and other attributes of specimens. This will be used for plotting and analyses in subsequent scripts.
classifier <- readr::read_csv(file = paste("./input_data/mammal_data.csv", sep = ""))


# Get attributes of classifier
attributes(classifier)


# extract relevant information from classifier
order <- factor(classifier$order) # specimen sex
clade <- factor(classifier$clade)
superorder <- factor(classifier$superorder)
family <- factor(classifier$family)
status <- factor(classifier$status)


# Get to know how many group levels are and save them
# as "lev.group" object. This will allow you to colour points for the PCA
lev.order <- length(levels(order))
lev.family <- length(levels(family))
lev.clade <- length(levels(clade))
lev.superorder <- length(levels(superorder))

# Assign colours for all Superorders (assigned in alphabetical order)
col.gp.1 <- c(
  "#e1d69c", "#5fc0a4", "#1e798a", "#13805d", "#f2464e",
  "#7c0c17", "#570407", "#fc951d", "#efc35a"
)

names(col.gp.1) <- levels(superorder)

# 3. Assign one colour to each specimen
col.gp <- col.gp.1[match(superorder, names(col.gp.1))]


# -----------------------------#
# PCA AND PLOT                 #
# -----------------------------#


# 4. Run plain PCA
PCA.3D.endo <- geomorph::gm.prcomp(Y.gpa.endo)

### Plot histogram of proportion of PC variances
summary(PCA.3D.endo)

pvar <- (PCA.3D.endo$sdev^2) / (sum(PCA.3D.endo$sdev^2)) ### graph of PC proportional variance
names(pvar) <- seq(1:length(pvar))
barplot(pvar[1:20], xlab = "Principal Components", ylab = "Proportion of Variance", ylim = c(0, 1))


# save outputs
save(PCA.3D.endo, file = "./output_data/PCA.3D.endo.rda")


# 1. Set labels
xlab <- paste("Principal Component 1 (", signif((PCA.3D.endo$d / sum(PCA.3D.endo$d) * 100), 3)[1], "%)", sep = "")

ylab <- paste("Principal Component 2 (", signif((PCA.3D.endo$d / sum(PCA.3D.endo$d) * 100), 3)[2], "%)", sep = "")

zlab <- paste("Principal Component 3 (", signif((PCA.3D.endo$d / sum(PCA.3D.endo$d) * 100), 3)[3], "%)", sep = "")


# Basic PCAs
# create data frame of PCs for ggplot
PCs <- data.frame(PCA.3D.endo$x,
  csize = Csize.endo, clade = classifier$superorder, species = classifier$species,
  locomotion = classifier$Locomotion, social = classifier$Sociality, activity = classifier$Activity,
  diet = classifier$Diet, superorder = classifier$superorder
)


# plot PCs 1 and 2 using ggplot
P1 <- ggplot(PCs, mapping = aes(x = PCs[, 1], y = PCs[, 3], fill = log(csize))) +
  geom_point(col = "black", shape = 21, size = 5) +
  labs(x = xlab, y = ylab) +
  # geom_text_repel(aes(label = dimnames(comb.dataset$original)[[3]], fontface = "italic"), color = "black") +
  scale_fill_viridis() +
  # scale_colour_manual(values = col.gp)+
  theme_bw()
plot(P1)


# plot PCs 1 and 2 using ggplot
P1 <- ggplot(PCs, mapping = aes(x = PCs[, 1], y = PCs[, 2], colour = classifier$locomotion_strict, fill = classifier$locomotion_strict)) +
  geom_convexhull(alpha = 0.2, show.legend = FALSE) +
  geom_point(col = "black", shape = 21, size = 5) +
  geom_text_repel(aes(label = dimnames(comb.dataset$original)[[3]], fontface = "italic"), color = "black") +
  labs(x = xlab, y = ylab) +
  coord_equal() +
  # scale_fill_manual(values = col.gp)+
  # scale_colour_manual(values = col.gp)+
  theme_bw()

plot(P1)

P1 <- ggplot(PCs, mapping = aes(x = PCs[, 1], y = PCs[, 2], colour = classifier$superorder, fill = classifier$superorder)) +
  geom_convexhull(alpha = 0.2, show.legend = FALSE) +
  geom_point(col = "black", shape = 21, size = 3.5) +
  geom_text_repel(aes(label = dimnames(comb.dataset$original)[[3]], fontface = "italic"), color = "black") +
  labs(x = xlab, y = ylab) +
  scale_fill_manual(values = col.gp) +
  scale_colour_manual(values = col.gp) +
  coord_equal() +
  theme_bw()

P2 <- ggplot(PCs, mapping = aes(x = PCs[, 1], y = PCs[, 3], colour = classifier$superorder, fill = classifier$superorder)) +
  geom_convexhull(alpha = 0.2, show.legend = FALSE) +
  geom_point(col = "black", shape = 21, size = 3.5) +
  geom_text_repel(aes(label = dimnames(comb.dataset$original)[[3]], fontface = "italic"), color = "black") +
  labs(x = xlab, y = zlab) +
  scale_fill_manual(values = col.gp) +
  scale_colour_manual(values = col.gp) +
  coord_equal() +
  theme_bw()

# plot PCAs
ggarrange(P1, P2, ncol = 1, nrow = 2, common.legend = TRUE, legend = "right")

P1 <- ggplot(PCs, mapping = aes(x = PCs[, 1], y = PCs[, 2], colour = classifier$activity, fill = classifier$activity)) +
  geom_convexhull(alpha = 0.2, show.legend = FALSE) +
  geom_point(col = "black", shape = 21, size = 5) +
  geom_text_repel(aes(label = dimnames(comb.dataset$original)[[3]], fontface = "italic"), color = "black") +
  labs(x = xlab, y = ylab) +
  coord_equal() +
  # scale_fill_manual(values = col.gp)+
  # scale_colour_manual(values = col.gp)+
  theme_bw()

plot(P1)

library(ggConvexHull)
# plot PCs 1 and 3 using ggplot
P2 <- ggplot(PCs, mapping = aes(x = PCs[, 1], y = PCs[, 3], colour = clade, fill = clade)) +
  geom_convexhull(alpha = 0.2, show.legend = FALSE) +
  geom_point(col = "black", shape = 21, size = 5) +
  coord_equal() +
  labs(x = xlab, y = zlab) +
  scale_fill_manual(values = col.gp) +
  scale_colour_manual(values = col.gp) +
  theme_bw()
plot(P2)


# plot PCAs
ggarrange(P1, P2, ncol = 1, nrow = 2, common.legend = TRUE, legend = "right")


####################################

#----- Plot phylomorphospaces -----#

####################################

# Import trees
tree_hyp1 <- read.nexus("./phylogeny_construction/final_tree_hyp1.nex")
tree_hyp2 <- read.nexus("./phylogeny_construction/final_tree_hyp2.nex")
tree_hyp3 <- read.nexus("./phylogeny_construction/final_tree_hyp3.nex")
tree_hyp4 <- read.nexus("./phylogeny_construction/final_tree_hyp4.nex")
tree_hyp5 <- read.nexus("./phylogeny_construction/final_tree_hyp5.nex")
tree_hyp6 <- read.nexus("./phylogeny_construction/final_tree_hyp6.nex")
tree_hyp7 <- read.nexus("./phylogeny_construction/final_tree_hyp7.nex")
tree_hyp8 <- read.nexus("./phylogeny_construction/final_tree_hyp8.nex")

plot(tree_hyp8)

# Check that tree tip names match shape data names
shape.names <- dimnames(Y.gpa.endo)[[3]]
tree.names <- tree_hyp1$tip.label
setdiff(shape.names, tree.names)


# Plot phylomorphospace for endocast
## create matrix for use in phytools::fastAnc()
mat1 <- cbind(eval(substitute(PCs$Comp1), PCs), eval(substitute(PCs$Comp2), PCs))
mat2 <- cbind(eval(substitute(PCs$Comp1), PCs), eval(substitute(PCs$Comp3), PCs))
rownames(mat1) <- dimnames(Y.gpa.endo)[[3]]
rownames(mat2) <- dimnames(Y.gpa.endo)[[3]]
stopifnot(length(setdiff(tree_hyp1$tip.label, rownames(mat1))) == 0)
stopifnot(length(setdiff(tree_hyp1$tip.label, rownames(mat2))) == 0)

xAnc <- fastAnc(tree_hyp1, mat1[, 1])
yAnc <- fastAnc(tree_hyp1, mat1[, 2])
zAnc <- fastAnc(tree_hyp1, mat2[, 2])

all_node_coords.1 <-
  data.frame(
    # put PC values for all nodes and tips in a dataframe
    # tips go first in order of tip labels, then numerical order for nodes
    x = c(mat1[tree_hyp1$tip.label, 1], xAnc),
    y = c(mat1[tree_hyp1$tip.label, 2], yAnc),
    nodeid = 1:(tree_hyp1$Nnode + length(tree_hyp1$tip.label))
  )
all_node_coords.2 <-
  data.frame(
    # put PC values for all nodes and tips in a dataframe
    # tips go first in order of tip labels, then numerical order for nodes
    x = c(mat2[tree_hyp1$tip.label, 1], xAnc),
    y = c(mat2[tree_hyp1$tip.label, 2], zAnc),
    nodeid = 1:(tree_hyp1$Nnode + length(tree_hyp1$tip.label))
  )

# get edge list from tree object
edges.1 <- data.frame(tree_hyp1$edge)
names(edges.1) <- c("node1", "node2")
# translate tip/node numbers into PC coordinates in all_node_coords dataframe
edgecoords.1 <- merge(
  merge(edges.1, all_node_coords.1, by.x = "node1", by.y = "nodeid"),
  all_node_coords.1,
  by.x = "node2", by.y = "nodeid"
)

# get edge list from tree object
edges.2 <- data.frame(tree_hyp1$edge)
names(edges.2) <- c("node1", "node2")
# translate tip/node numbers into PC coordinates in all_node_coords dataframe
edgecoords.2 <- merge(
  merge(edges.2, all_node_coords.2, by.x = "node1", by.y = "nodeid"),
  all_node_coords.2,
  by.x = "node2", by.y = "nodeid"
)


# col.gp.1 <- c("orchid2", "light pink" , "blue3", "black","cadetblue1","chartreuse2" , "red1","dodgerblue",
#               "darkgreen", "darkorange2", "goldenrod1", "purple3", "yellow","darkred")
library(plotly)
PCs <- PCs %>% mutate(species = dimnames(Y.gpa.endo)[[3]])

phylopic_df <- data.frame(
  species = c(
    "Chlamyphorus_truncatus",
    "Homo_sapiens",
    "Mormoops_megalophylla",
    "Pongo_abelii",
    "Delphinus_delphis",
    "Odobenus_rosmarus",
    "Echinops_telfairi",
    "Kryptobaatar_dashzevegi"
  ),
  picnum = c(
    "06719808-32ca-4258-b2ec-9c7e1b1c24bb",
    "9fae30cd-fb59-4a81-a39c-e1826a35f612",
    "f137fec2-0f1b-4745-afac-a7f3343b88ce",
    "67144c22-93c2-4dc0-ba13-9f9dd2d223b9",
    "3caf4fbd-ca3a-48b4-925a-50fbe9acd887",
    "d2575005-1fcb-4a86-8c83-e3bda619adf2",
    "a89e224a-e374-4be0-a9bb-3481b16f87fb",
    "63b5a09a-e39c-43d2-9d92-0e3c800ea5b6"
  )
)


label_data <- dplyr::left_join(phylopic_df, PCs)
# label_data$species<-gsub("_"," ", label_data$species)
# limited dataset
PhM.whole.1 <- ggplot() + # PCs 1 and 2
  geom_segment(data = edgecoords.1, aes(x = x.x, xend = x.y, y = y.x, yend = y.y), size = 0.25) +
  geom_convexhull(data = PCs, aes(x = Comp1, y = Comp2, colour = superorder, fill = superorder), alpha = 0.2) + # ,show.legend = FALSE)+
  geom_point(data = PCs, aes(x = Comp1, y = Comp2, colour = superorder, fill = superorder, shape = status), colour = "black", shape = 21, stroke = .35, size = 2.7) + # ,show.legend = FALSE)+
  # geom_text_repel(data = label_data, aes(x=Comp1, y=Comp2, label=species, fontface = "italic"))+
  labs(x = xlab, y = ylab) +
  scale_fill_manual(values = col.gp) +
  scale_colour_manual(values = col.gp) +
  scale_x_continuous(expand = expansion(mult = c(0.09, 0.05))) +
  scale_y_continuous(expand = expansion(mult = c(0.10, 0.05))) +
  # geom_text(aes(label = shape.names), show.legend = FALSE)+
  theme_bw() +
  coord_fixed() +
  theme(
    legend.position = "left",
    legend.justification = c(0.5, -3),
    plot.margin = margin(.01, .9, 9, .05, unit = "lines")
  )
PhM.whole.1 <- PhM.whole.1 + geom_phylopic(data = label_data[c(2, 4), ], aes(
  x = Comp1 - .015,
  y = Comp2 - .01,
  color = superorder,
  uuid = picnum
), height = .07) +
  geom_phylopic(data = label_data[c(3), ], aes(
    x = Comp1 - .015,
    y = Comp2 - .03,
    color = superorder,
    uuid = picnum
  ), height = .04) +
  geom_phylopic(data = label_data[c(1), ], aes(
    x = Comp1 - .015,
    y = Comp2 - .02,
    color = superorder,
    uuid = picnum
  ), height = .03) +
  geom_phylopic(data = label_data[c(5), ], aes(
    x = Comp1 - .05,
    y = Comp2 - .01,
    color = superorder,
    uuid = picnum
  ), height = .035) +
  geom_phylopic(data = label_data[c(6), ], aes(
    x = Comp1 - .03,
    y = Comp2 - .05,
    color = superorder,
    uuid = picnum
  ), height = .035) +
  geom_phylopic(data = label_data[c(7), ], aes(
    x = Comp1 + .03,
    y = Comp2 - .01,
    color = superorder,
    uuid = picnum
  ), height = .03) +
  geom_phylopic(data = label_data[c(8), ], aes(
    x = Comp1 + .04,
    y = Comp2 - .01,
    color = superorder,
    uuid = picnum
  ), height = .03)


plot(PhM.whole.1)
ggsave("./Images/Figures/NewFig1.pdf", PhM.whole.1, height = 7, width = 7, units = "in")


PhM.whole.2 <- ggplot() + # PCs 1 and 2
  geom_segment(data = edgecoords.2, aes(x = x.x, xend = x.y, y = y.x, yend = y.y), size = 0.25) +
  geom_convexhull(data = PCs, aes(x = Comp1, y = Comp3, colour = superorder, fill = superorder), alpha = 0.2) + # ,show.legend = FALSE)+
  geom_point(data = PCs, aes(x = Comp1, y = Comp3, colour = superorder, fill = superorder, shape = status), colour = "black", shape = 21, stroke = .35, size = 2.7) + # ,show.legend = FALSE)+
  # geom_text_repel(data = label_data, aes(x=Comp1, y=Comp2, label=species, fontface = "italic"))+
  labs(x = xlab, y = zlab) +
  scale_fill_manual(values = col.gp) +
  scale_colour_manual(values = col.gp) +
  # geom_text(aes(label = shape.names), show.legend = FALSE)+
  theme_bw() +
  coord_fixed() +
  theme(
    legend.position = "left",
    legend.justification = c(0.5, -9),
    plot.margin = margin(.01, .9, 9, .05, unit = "lines")
  )
plot(PhM.whole.2)


ggsave("./Images/Figures/NewFig2.pdf", PhM.whole.2, height = 7, width = 7, units = "in")


# make interactive PCA plot
p <- ggplot(PCs, aes(x = Comp1, y = Comp2, label = species)) +
  geom_point(data = PCs, aes(x = Comp1, y = Comp2, colour = superorder, fill = superorder, shape = status), colour = "black", shape = 21, size = 3) +
  xlab("PC1") +
  ylab("PC2") +
  scale_fill_manual(values = col.gp) +
  scale_colour_manual(values = col.gp) +
  theme_bw()
p
ggplotly(p)

library(plotly)
p <- ggplot(PCs, aes(Comp1, Comp2))
p <- p + geom_point()

ggplotly(p)


#### Plot max and min shapes along PCs

template.mesh.pc <- vcgImport(paste("./ply_ASCII/Abrocoma_cinerea_endocast.ply",
  sep = ""
))

# Species_mesh that you already defined to vizualise the landmarks
Species_landmark <- arr[, , "Abrocoma_cinerea"] # same

# PC min and max
PC1max <- PCA.3D.endo$shapes$shapes.comp1$max
PC1min <- PCA.3D.endo$shapes$shapes.comp1$min
PC2max <- PCA.3D.endo$shapes$shapes.comp2$max
PC2min <- PCA.3D.endo$shapes$shapes.comp2$min
PC3max <- PCA.3D.endo$shapes$shapes.comp3$max
PC3min <- PCA.3D.endo$shapes$shapes.comp3$min


# plot template specimen with landmarks coloured by region
shade3d(atlas_olf$mesh, col = "cornflowerblue", specular = "black")
# spheres3d(arr[,,"Abrocoma_cinerea"], radius = 0.2, col = "grey")
spheres3d(arr[neocortex, , "Abrocoma_cinerea"], col = "#F28E2B", radius = 0.5)
spheres3d(arr[olfactory, , "Abrocoma_cinerea"], col = "#B07AA1", radius = 0.5)
spheres3d(arr[cerebellum, , "Abrocoma_cinerea"], col = "#59A14F", radius = 0.5)
spheres3d(arr[brainstem, , "Abrocoma_cinerea"], col = "#EDC948", radius = 0.5)

??spheres3d
str(PC1max)
# max PC1
warpskull <- tps3d(template.mesh.pc, Species_landmark, PC1max)
shade3d(warpskull, col = "cornflowerblue", specular = "black")
spheres3d(PC1max[neocortex, ], col = "#F28E2B", radius = 0.003)
spheres3d(PC1max[olfactory, ], col = "#B07AA1", radius = 0.003)
spheres3d(PC1max[cerebellum, ], col = "#59A14F", radius = 0.003)
spheres3d(PC1max[brainstem, ], col = "#EDC948", radius = 0.003)
rgl.snapshot(filename = "PC1max_lateral.png")
rgl.snapshot(filename = "PC1max_posterior.png")
# min PC1
warpskull <- tps3d(template.mesh.pc, Species_landmark, PC1min)
shade3d(warpskull, col = "cornflowerblue", specular = "black")
spheres3d(PC1min[neocortex, ], col = "#F28E2B", radius = 0.003)
spheres3d(PC1min[olfactory, ], col = "#B07AA1", radius = 0.003)
spheres3d(PC1min[cerebellum, ], col = "#59A14F", radius = 0.003)
spheres3d(PC1min[brainstem, ], col = "#EDC948", radius = 0.003)
rgl.snapshot(filename = "PC1min_lateral.png")
rgl.snapshot(filename = "PC1min_posterior.png")
# max PC2
warpskull <- tps3d(template.mesh.pc, Species_landmark, PC2max)
shade3d(warpskull, col = "cornflowerblue", specular = "black")
spheres3d(PC2max[neocortex, ], col = "#F28E2B", radius = 0.003)
spheres3d(PC2max[olfactory, ], col = "#B07AA1", radius = 0.003)
spheres3d(PC2max[cerebellum, ], col = "#59A14F", radius = 0.003)
spheres3d(PC2max[brainstem, ], col = "#EDC948", radius = 0.003)
rgl.snapshot(filename = "PC2max_lateral.png")
rgl.snapshot(filename = "PC2max_posterior.png")
# min PC2
warpskull <- tps3d(template.mesh.pc, Species_landmark, PC2min)
shade3d(warpskull, col = "cornflowerblue", specular = "black")
spheres3d(PC2min[neocortex, ], col = "#F28E2B", radius = 0.003)
spheres3d(PC2min[olfactory, ], col = "#B07AA1", radius = 0.003)
spheres3d(PC2min[cerebellum, ], col = "#59A14F", radius = 0.003)
spheres3d(PC2min[brainstem, ], col = "#EDC948", radius = 0.003)
rgl.snapshot(filename = "PC2min_lateral.png")
rgl.snapshot(filename = "PC2min_posterior.png")
# max PC3
warpskull <- tps3d(template.mesh.pc, Species_landmark, PC3max)
shade3d(warpskull, col = "cornflowerblue", specular = "black")
spheres3d(PC3max[neocortex, ], col = "#F28E2B", radius = 0.003)
spheres3d(PC3max[olfactory, ], col = "#B07AA1", radius = 0.003)
spheres3d(PC3max[cerebellum, ], col = "#59A14F", radius = 0.003)
spheres3d(PC3max[brainstem, ], col = "#EDC948", radius = 0.003)
rgl.snapshot(filename = "PC3max_lateral.png")
rgl.snapshot(filename = "PC2min_posterior.png")
# min PC3
warpskull <- tps3d(template.mesh.pc, Species_landmark, PC3min)
shade3d(warpskull, col = "cornflowerblue", specular = "black")
spheres3d(PC3min[neocortex, ], col = "#F28E2B", radius = 0.003)
spheres3d(PC3min[olfactory, ], col = "#B07AA1", radius = 0.003)
spheres3d(PC3min[cerebellum, ], col = "#59A14F", radius = 0.003)
spheres3d(PC3min[brainstem, ], col = "#EDC948", radius = 0.003)
rgl.snapshot(filename = "PC3min_lateral.png")
rgl.snapshot(filename = "PC2min_posterior.png")
