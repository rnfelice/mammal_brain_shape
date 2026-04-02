# 3. Convergence testing
# This script tests convergence among groups of species
library(geomorph)
library(convevol)
library(ape)
library(tidyverse)


# load the data and trees
# aligned specimens
n_landmarks <- 122
n_dimensions <- 3
raw_data <- read.csv("./input_data/procrutes_aligned_coords.csv", row.names = 1)
Y.gpa.endo <- arrayspecs(raw_data, p = n_landmarks, k = n_dimensions)

PCA.3D.endo <- geomorph::gm.prcomp(Y.gpa.endo)
tree_hyp1 <- read.nexus("./phylogeny_construction/final_tree_hyp1.nex")
# 3.1 Designate groups for convergence testing

burrowers <- c(
  "Notoryctes_typhlops", "Amblysomus_hottentotus", "Talpa_europaea",
  "Chlamyphorus_truncatus"
)

aquatic <- c("Enhydra_lutris", "Neomonachus_tropicalis", "Inia_geoffrensis")


# run convSigCt convergence check, using first 3 PC scores

# Forelimb burrowers

conv.burrow <- convSigCt(final_tree, PCA.3D.endo$x[, c(1:3)],
  focaltaxa = burrowers, nsim = 100
)

plotCt(output = conv.burrow, phy = tree_hyp1, focaltaxa = burrowers)


# Aquatic taxa

conv.aquatic <- convSigCt(tree_hyp1, PCA.3D.endo$x[, c(1:3)],
  focaltaxa = aquatic, nsim = 100
)

plotCt(output = conv.aquatic, phy = tree_hyp1, focaltaxa = aquatic)


# Plot endocast meshes

# Burrowers

# 1.1.3 Import mesh file for Notoryctes_typhlops
Notoryctes.mesh <- vcgImport(paste("./ply_ASCII/Notoryctes_typhlops_endocast.ply",
  sep = ""
))
# plot mesh
shade3d(Notoryctes.mesh, col = "cornflowerblue", specular = "black")


# 1.1.3 Import mesh file for Amblysomus_hottentotus
Amblysomus.mesh <- vcgImport(paste("./ply_ASCII/Amblysomus_hottentotus_endocast.ply",
  sep = ""
))
# plot mesh
shade3d(Amblysomus.mesh, col = "cornflowerblue", specular = "black")


# 1.1.3 Import mesh file for Chlamyphorus_truncatus
Chlamyphorus.mesh <- vcgImport(paste("./ply_ASCII/Chlamyphorus_truncatus_endocast.ply",
  sep = ""
))
# plot mesh
shade3d(Chlamyphorus.mesh, col = "cornflowerblue", specular = "black")


# 1.1.3 Import mesh file for Talpa_europaea
Talpa.mesh <- vcgImport(paste("./ply_ASCII/Talpa_europaea_endocast.ply",
  sep = ""
))
# plot mesh
shade3d(Talpa.mesh, col = "cornflowerblue", specular = "black")


# Swimmers

# 1.1.3 Import mesh file for Enhydra_lutris
Enhydra.mesh <- vcgImport(paste("./ply_ASCII/Enhydra_lutris_endocast.ply",
  sep = ""
))
# plot mesh
shade3d(Enhydra.mesh, col = "cornflowerblue", specular = "black")


# 1.1.3 Import mesh file for Neomonachus_tropicalis
Neomonachus.mesh <- vcgImport(paste("./ply_ASCII/Neomonachus_tropicalis_endocast.ply",
  sep = ""
))
# plot mesh
shade3d(Neomonachus.mesh, col = "cornflowerblue", specular = "black")


# 1.1.3 Import mesh file for Inia_geoffrensis
Inia.mesh <- vcgImport(paste("./ply_ASCII/Inia_geoffrensis_endocast.ply",
  sep = ""
))
# plot mesh
shade3d(Inia.mesh, col = "cornflowerblue", specular = "black")


# Plot PCA highlighting putative convergent taxa
conv <- c(burrowers, aquatic)

col.conv.1 <- c(
  "#f2464e", "#e1d69c", "#1e798a", "#efc35a",
  "#1e798a", "#1e798a", "#1e798a"
)

names(col.conv.1) <- conv

col.conv <- col.conv.1[match(conv, names(col.conv.1))]


# plot PCA with convergent taxa
PhM.conv.1 <- ggplot() + # PCs 1 and 2
  geom_segment(data = edgecoords.1, aes(x = x.x, xend = x.y, y = y.x, yend = y.y), col = "grey", size = 0.5) +
  # geom_convexhull(data = PCs, aes(x=Comp1, y = Comp2, colour = superorder, fill = superorder), alpha = 0.2)+#,show.legend = FALSE)+
  geom_point(data = PCs, aes(x = Comp1, y = Comp2, fill = col.conv), fill = "grey", colour = "grey", shape = 21, size = 3) + # ,show.legend = FALSE)+
  labs(x = xlab, y = ylab) +
  # scale_fill_manual(values = col.gp)+
  # scale_colour_manual(values = col.gp)+
  # geom_text(aes(label = shape.names), show.legend = FALSE)+
  theme_bw()
plot(PhM.conv.1)
