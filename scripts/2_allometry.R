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


# Allometry analyses

# Test whole brain shape allometry using each tree hypothesis
gdf <- geomorph.data.frame(coords = Y.gpa.endo, mass = classifier$Mass, Csize = Csize.endo, superorder = classifier$superorder)

names(gdf) <- c("coords", "mass", "Csize", "superorder")


# 1.2.1 phylogenetic allometry of whole skull shape against log body mass

# hyp 1
Allometry.Csize_hyp1 <- procD.pgls(coords ~ log(Csize),
    f2 = NULL, f3 = NULL, phy = tree_hyp1, logsz = FALSE, data = gdf,
    iter = 999, print.progress = FALSE
)
summary(Allometry.Csize_hyp1)

Allometry.mass_hyp1 <- procD.pgls(coords ~ log(mass),
    f2 = NULL, f3 = NULL, phy = tree_hyp1, logsz = FALSE, data = gdf,
    iter = 999, print.progress = FALSE
)
summary(Allometry.mass_hyp1)


# hyp 2
Allometry.Csize_hyp2 <- procD.pgls(coords ~ log(Csize),
    f2 = NULL, f3 = NULL, phy = tree_hyp2, logsz = FALSE, data = gdf,
    iter = 999, print.progress = FALSE
)
summary(Allometry.Csize_hyp2)

Allometry.mass_hyp2 <- procD.pgls(coords ~ log(mass),
    f2 = NULL, f3 = NULL, phy = tree_hyp2, logsz = FALSE, data = gdf,
    iter = 999, print.progress = FALSE
)
summary(Allometry.mass_hyp2)


# hyp 3
Allometry.Csize_hyp3 <- procD.pgls(coords ~ log(Csize),
    f2 = NULL, f3 = NULL, phy = tree_hyp3, logsz = FALSE, data = gdf,
    iter = 999, print.progress = FALSE
)
summary(Allometry.Csize_hyp3)

Allometry.mass_hyp3 <- procD.pgls(coords ~ log(mass),
    f2 = NULL, f3 = NULL, phy = tree_hyp3, logsz = FALSE, data = gdf,
    iter = 999, print.progress = FALSE
)
summary(Allometry.mass_hyp3)


# hyp 4
Allometry.Csize_hyp4 <- procD.pgls(coords ~ log(Csize),
    f2 = NULL, f3 = NULL, phy = tree_hyp4, logsz = FALSE, data = gdf,
    iter = 999, print.progress = FALSE
)
summary(Allometry.Csize_hyp4)

Allometry.mass_hyp4 <- procD.pgls(coords ~ log(mass),
    f2 = NULL, f3 = NULL, phy = tree_hyp4, logsz = FALSE, data = gdf,
    iter = 999, print.progress = FALSE
)
summary(Allometry.mass_hyp4)


# hyp 5
Allometry.Csize_hyp5 <- procD.pgls(coords ~ log(Csize),
    f2 = NULL, f3 = NULL, phy = tree_hyp5, logsz = FALSE, data = gdf,
    iter = 999, print.progress = FALSE
)
summary(Allometry.Csize_hyp5)

Allometry.mass_hyp5 <- procD.pgls(coords ~ log(mass),
    f2 = NULL, f3 = NULL, phy = tree_hyp5, logsz = FALSE, data = gdf,
    iter = 999, print.progress = FALSE
)
summary(Allometry.mass_hyp5)


# hyp 6
Allometry.Csize_hyp6 <- procD.pgls(coords ~ log(Csize),
    f2 = NULL, f3 = NULL, phy = tree_hyp6, logsz = FALSE, data = gdf,
    iter = 999, print.progress = FALSE
)
summary(Allometry.Csize_hyp6)

Allometry.mass_hyp6 <- procD.pgls(coords ~ log(mass),
    f2 = NULL, f3 = NULL, phy = tree_hyp6, logsz = FALSE, data = gdf,
    iter = 999, print.progress = FALSE
)
summary(Allometry.mass_hyp6)


# hyp 7
Allometry.Csize_hyp7 <- procD.pgls(coords ~ log(Csize),
    f2 = NULL, f3 = NULL, phy = tree_hyp7, logsz = FALSE, data = gdf,
    iter = 999, print.progress = FALSE
)
summary(Allometry.Csize_hyp7)

Allometry.mass_hyp7 <- procD.pgls(coords ~ log(mass),
    f2 = NULL, f3 = NULL, phy = tree_hyp7, logsz = FALSE, data = gdf,
    iter = 999, print.progress = FALSE
)
summary(Allometry.mass_hyp7)


# hyp 8
Allometry.Csize_hyp8 <- procD.pgls(coords ~ log(Csize),
    f2 = NULL, f3 = NULL, phy = tree_hyp8, logsz = FALSE, data = gdf,
    iter = 999, print.progress = FALSE
)
summary(Allometry.Csize_hyp8)

Allometry.mass_hyp8 <- procD.pgls(coords ~ log(mass),
    f2 = NULL, f3 = NULL, phy = tree_hyp8, logsz = FALSE, data = gdf,
    iter = 999, print.progress = FALSE
)
summary(Allometry.mass_hyp8)


# Neocortex allometry
gdf.neocortex <- geomorph.data.frame(coords = Y.gpa.endo[neocortex, , ], mass = classifier$Mass, Csize = Csize.endo, superorder = classifier$superorder)

names(gdf.neocortex) <- c("coords", "mass", "Csize", "superorder")

# 1.2.1 phylogenetic allometry of whole skull shape against log body mass

# hyp 1 Csize
Allometry.neocortex.hyp1 <- procD.pgls(coords ~ log(Csize),
    f2 = NULL, f3 = NULL, phy = tree_hyp1, logsz = FALSE, data = gdf.neocortex,
    iter = 999, print.progress = FALSE
)
summary(Allometry.neocortex.hyp1)

# hyp 1 mass
Allometry.neocortex.mass <- procD.pgls(coords ~ log(mass),
    f2 = NULL, f3 = NULL, phy = tree_hyp1, logsz = FALSE, data = gdf.neocortex,
    iter = 999, print.progress = FALSE
)
summary(Allometry.neocortex.mass)


# Cerebellum allometry
gdf.cerebellum <- geomorph.data.frame(coords = Y.gpa.endo[cerebellum, , ], mass = classifier$Mass, Csize = Csize.endo, superorder = classifier$superorder)

names(gdf.cerebellum) <- c("coords", "mass", "Csize", "superorder")

# 1.2.1 phylogenetic allometry of whole skull shape against log body mass

# hyp 1
Allometry.cerebellum.hyp1 <- procD.pgls(coords ~ log(Csize),
    f2 = NULL, f3 = NULL, phy = tree_hyp1, logsz = FALSE, data = gdf.cerebellum,
    iter = 999, print.progress = FALSE
)
summary(Allometry.cerebellum.hyp1)

# hyp 1 mass
Allometry.cerebellum.mass <- procD.pgls(coords ~ log(mass),
    f2 = NULL, f3 = NULL, phy = tree_hyp1, logsz = FALSE, data = gdf.cerebellum,
    iter = 999, print.progress = FALSE
)
summary(Allometry.cerebellum.mass)


# Brainstem allometry
gdf.brainstem <- geomorph.data.frame(coords = Y.gpa.endo[brainstem, , ], mass = classifier$Mass, Csize = Csize.endo, superorder = classifier$superorder)

names(gdf.brainstem) <- c("coords", "mass", "Csize", "superorder")

# 1.2.1 phylogenetic allometry of whole skull shape against log body mass

# hyp 1
Allometry.brainstem.hyp1 <- procD.pgls(coords ~ log(Csize),
    f2 = NULL, f3 = NULL, phy = tree_hyp1, logsz = FALSE, data = gdf.brainstem,
    iter = 999, print.progress = FALSE
)
summary(Allometry.brainstem.hyp1)

# hyp 1 mass
Allometry.brainstem.mass <- procD.pgls(coords ~ log(mass),
    f2 = NULL, f3 = NULL, phy = tree_hyp1, logsz = FALSE, data = gdf.brainstem,
    iter = 999, print.progress = FALSE
)
summary(Allometry.brainstem.mass)


# Olfactory allometry
gdf.olfactory <- geomorph.data.frame(coords = Y.gpa.endo[olfactory, , ], mass = classifier$Mass, Csize = Csize.endo, superorder = classifier$superorder)

names(gdf.olfactory) <- c("coords", "mass", "Csize", "superorder")

# 1.2.1 phylogenetic allometry of whole skull shape against log body mass

# hyp 1
Allometry.olfactory.hyp1 <- procD.pgls(coords ~ log(Csize),
    f2 = NULL, f3 = NULL, phy = tree_hyp1, logsz = FALSE, data = gdf.olfactory,
    iter = 999, print.progress = FALSE
)
summary(Allometry.olfactory.hyp1)

# hyp 1 mass
Allometry.olfactory.mass <- procD.pgls(coords ~ log(mass),
    f2 = NULL, f3 = NULL, phy = tree_hyp1, logsz = FALSE, data = gdf.olfactory,
    iter = 999, print.progress = FALSE
)
summary(Allometry.olfactory.mass)
