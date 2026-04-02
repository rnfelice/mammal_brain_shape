library(geomorph)
library(tidyverse)
library(ggConvexHull)
library(ggpubr)

# load the data
n_landmarks <- 122
n_dimensions <- 3
raw_data <- read.csv("./input_data/procrutes_aligned_coords.csv", row.names = 1)
Y.gpa.endo <- arrayspecs(raw_data, p = n_landmarks, k = n_dimensions)

Csize <- read.csv("./input_data/centroid_sizes.csv", row.names = 1)
Csize.endo <- setNames((Csize[, 1]), rownames(Csize))

classifier <- readr::read_csv(file = paste("./input_data/mammal_data.csv", sep = ""))

# Designate brain regions
neocortex <- c(2:3, 10:21, 52:59, 74:113)
olfactory <- c(1, 4, 22:29, 114:122)
cerebellum <- c(5, 9, 30:32, 45:51, 60:68)
brainstem <- c(6:8, 33:44, 69:73)

# Scripts for plotting PCAs of separate brain regions
PCA.3D.neocortex <- geomorph::gm.prcomp(Y.gpa.endo[neocortex, , ])
PCA.3D.olfactory <- geomorph::gm.prcomp(Y.gpa.endo[olfactory, , ])
PCA.3D.cerebellum <- geomorph::gm.prcomp(Y.gpa.endo[cerebellum, , ])
PCA.3D.brainstem <- geomorph::gm.prcomp(Y.gpa.endo[brainstem, , ])


# 1. Set labels
xlab.neocortex <- paste("Principal Component 1 (", signif((PCA.3D.neocortex$d / sum(PCA.3D.neocortex$d) * 100), 3)[1], "%)", sep = "")
ylab.neocortex <- paste("Principal Component 2 (", signif((PCA.3D.neocortex$d / sum(PCA.3D.neocortex$d) * 100), 3)[2], "%)", sep = "")


# 1. Set labels
xlab.olfactory <- paste("Principal Component 1 (", signif((PCA.3D.olfactory$d / sum(PCA.3D.olfactory$d) * 100), 3)[1], "%)", sep = "")
ylab.olfactory <- paste("Principal Component 2 (", signif((PCA.3D.olfactory$d / sum(PCA.3D.olfactory$d) * 100), 3)[2], "%)", sep = "")


# 1. Set labels
xlab.cerebellum <- paste("Principal Component 1 (", signif((PCA.3D.cerebellum$d / sum(PCA.3D.cerebellum$d) * 100), 3)[1], "%)", sep = "")
ylab.cerebellum <- paste("Principal Component 2 (", signif((PCA.3D.cerebellum$d / sum(PCA.3D.cerebellum$d) * 100), 3)[2], "%)", sep = "")


# 1. Set labels
xlab.brainstem <- paste("Principal Component 1 (", signif((PCA.3D.brainstem$d / sum(PCA.3D.brainstem$d) * 100), 3)[1], "%)", sep = "")
ylab.brainstem <- paste("Principal Component 2 (", signif((PCA.3D.brainstem$d / sum(PCA.3D.brainstem$d) * 100), 3)[2], "%)", sep = "")


PC.neocortex <- data.frame(PCA.3D.neocortex$x, csize = Csize.endo, clade = classifier$superorder)
PC.olfactory <- data.frame(PCA.3D.olfactory$x, csize = Csize.endo, clade = classifier$superorder)
PC.cerebellum <- data.frame(PCA.3D.cerebellum$x, csize = Csize.endo, clade = classifier$superorder)
PC.brainstem <- data.frame(PCA.3D.brainstem$x, csize = Csize.endo, clade = classifier$superorder)


# prep colors for plotting
# Assign colours for all Superorders (assigned in alphabetical order)
superorder <- factor(classifier$superorder)
col.gp.1 <- c(
  "#e1d69c", "#5fc0a4", "#1e798a", "#13805d", "#f2464e",
  "#7c0c17", "#570407", "#fc951d", "#efc35a"
)

names(col.gp.1) <- levels(superorder)

# 3. Assign one colour to each specimen
col.gp <- col.gp.1[match(superorder, names(col.gp.1))]


P.neocortex <- ggplot(PC.neocortex, mapping = aes(x = PC.neocortex[, 1], y = PC.neocortex[, 2], colour = superorder, fill = superorder)) +
  geom_convexhull(alpha = 0.2, show.legend = FALSE) +
  geom_point(col = "black", shape = 21, size = 3.5) +
  # geom_text_repel(aes(label = dimnames(comb.dataset$original)[[3]], fontface = "italic"), color = "black") +
  labs(x = xlab.neocortex, y = ylab.neocortex) +
  scale_fill_manual(values = col.gp) +
  scale_colour_manual(values = col.gp) +
  ggtitle("Neocortex") +
  theme_bw()
plot(P.neocortex)


P.olfactory <- ggplot(PC.olfactory, mapping = aes(x = PC.olfactory[, 1], y = PC.olfactory[, 2], colour = superorder, fill = superorder)) +
  geom_convexhull(alpha = 0.2, show.legend = FALSE) +
  geom_point(col = "black", shape = 21, size = 3.5) +
  # geom_text_repel(aes(label = dimnames(comb.dataset$original)[[3]], fontface = "italic"), color = "black") +
  labs(x = xlab.olfactory, y = ylab.olfactory) +
  scale_fill_manual(values = col.gp) +
  scale_colour_manual(values = col.gp) +
  ggtitle("Olfactory") +
  theme_bw()
plot(P.olfactory)


P.cerebellum <- ggplot(PC.cerebellum, mapping = aes(x = PC.cerebellum[, 1], y = PC.cerebellum[, 2], colour = superorder, fill = superorder)) +
  geom_convexhull(alpha = 0.2, show.legend = FALSE) +
  geom_point(col = "black", shape = 21, size = 3.5) +
  # geom_text_repel(aes(label = dimnames(comb.dataset$original)[[3]], fontface = "italic"), color = "black") +
  labs(x = xlab.cerebellum, y = ylab.cerebellum) +
  scale_fill_manual(values = col.gp) +
  scale_colour_manual(values = col.gp) +
  ggtitle("Cerebellum") +
  theme_bw()
plot(P.cerebellum)


P.brainstem <- ggplot(PC.brainstem, mapping = aes(x = PC.brainstem[, 1], y = PC.brainstem[, 2], colour = superorder, fill = superorder)) +
  geom_convexhull(alpha = 0.2, show.legend = FALSE) +
  geom_point(col = "black", shape = 21, size = 3.5) +
  # geom_text_repel(aes(label = dimnames(comb.dataset$original)[[3]], fontface = "italic"), color = "black") +
  labs(x = xlab.brainstem, y = ylab.brainstem) +
  scale_fill_manual(values = col.gp) +
  scale_colour_manual(values = col.gp) +
  ggtitle("Brainstem") +
  theme_bw()
plot(P.brainstem)

ggarrange(P.neocortex, P.olfactory, P.cerebellum, P.brainstem, ncol = 2, nrow = 2, common.legend = TRUE, legend = "right")


#### Plot max and min shapes along PCs

draw_segments <- function(links, landmarks, col = "black", lwd = 1, antialiasing = T) {
  for (i in 1:nrow(links)) {
    point1 <- landmarks[links[i, 1], ]
    point2 <- landmarks[links[i, 2], ]
    points <- rbind(point1, point2)
    segments3d(points, col = col, lwd = lwd, line_antialias = antialiasing)
  }
}

olfactory_links <- read.csv("./input_data/plotting_links/olfactory.links.csv", header = F)
neocortex_links <- read.csv("./input_data/plotting_links/neocortex.links.csv", header = F)
cerebellum_links <- read.csv("./input_data/plotting_links/cerebellum.links.csv", header = F)
brainstem_links <- read.csv("./input_data/plotting_links/brainstem.links.csv", header = F)


# Neocortex

# PC min and max
PC1max.neocortex <- PCA.3D.neocortex$shapes$shapes.comp1$max
PC1min.neocortex <- PCA.3D.neocortex$shapes$shapes.comp1$min
PC2max.neocortex <- PCA.3D.neocortex$shapes$shapes.comp2$max
PC2min.neocortex <- PCA.3D.neocortex$shapes$shapes.comp2$min


# maybe:
# open3d(antialias=4)
spheres3d(PC1max.neocortex, col = "#009E73", radius = 0.005)
draw_segments(links = neocortex_links, landmarks = PC1max.neocortex, lwd = 5)

## Save RGL parameters to a list object
pp1 <- par3d(no.readonly = TRUE)
## Save the list to a text file
dput(pp1, file = "neocortex.R", control = "all")

## Then, in a later session, to recreate the plot just as you had it:

library(rgl)
pp <- dget("neocortex.R")
plot3d(iris)
open3d()
par3d(pp1)


spheres3d(PC1min.neocortex, col = "#009E73", radius = 0.005)
draw_segments(links = neocortex_links, landmarks = PC1min.neocortex, lwd = 5)
spheres3d(PC2max.neocortex, col = "#009E73", radius = 0.005)
draw_segments(links = neocortex_links, landmarks = PC2max.neocortex, lwd = 5)
spheres3d(PC2min.neocortex, col = "#009E73", radius = 0.005)
draw_segments(links = neocortex_links, landmarks = PC2min.neocortex, lwd = 5)

# Olfactory bulb

# PC min and max
PC1max.olfactory <- PCA.3D.olfactory$shapes$shapes.comp1$max
PC1min.olfactory <- PCA.3D.olfactory$shapes$shapes.comp1$min
PC2max.olfactory <- PCA.3D.olfactory$shapes$shapes.comp2$max
PC2min.olfactory <- PCA.3D.olfactory$shapes$shapes.comp2$min


spheres3d(PC1max.olfactory, col = "#F0E442", radius = 0.003)
draw_segments(links = olfactory_links, landmarks = PC1max.olfactory, lwd = 5)
spheres3d(PC1min.olfactory, col = "#F0E442", radius = 0.003)
draw_segments(links = olfactory_links, landmarks = PC1min.olfactory, lwd = 5)
spheres3d(PC2max.olfactory, col = "#F0E442", radius = 0.003)
draw_segments(links = olfactory_links, landmarks = PC2max.olfactory, lwd = 5)
spheres3d(PC2min.olfactory, col = "#F0E442", radius = 0.003)
draw_segments(links = olfactory_links, landmarks = PC2min.olfactory, lwd = 5)

# Cerebellum

# PC min and max
PC1max.cerebellum <- PCA.3D.cerebellum$shapes$shapes.comp1$max
PC1min.cerebellum <- PCA.3D.cerebellum$shapes$shapes.comp1$min
PC2max.cerebellum <- PCA.3D.cerebellum$shapes$shapes.comp2$max
PC2min.cerebellum <- PCA.3D.cerebellum$shapes$shapes.comp2$min

spheres3d(PC1max.cerebellum, col = "#56B4E9", radius = 0.005)
draw_segments(links = cerebellum_links, landmarks = PC1max.cerebellum, lwd = 5)
spheres3d(PC1min.cerebellum, col = "#56B4E9", radius = 0.005)
draw_segments(links = cerebellum_links, landmarks = PC1min.cerebellum, lwd = 5)
spheres3d(PC2max.cerebellum, col = "#56B4E9", radius = 0.005)
draw_segments(links = cerebellum_links, landmarks = PC2max.cerebellum, lwd = 5)
spheres3d(PC2min.cerebellum, col = "#56B4E9", radius = 0.005)
draw_segments(links = cerebellum_links, landmarks = PC2min.cerebellum, lwd = 5)


# Brainstem

# PC min and max
PC1max.brainstem <- PCA.3D.brainstem$shapes$shapes.comp1$max
PC1min.brainstem <- PCA.3D.brainstem$shapes$shapes.comp1$min
PC2max.brainstem <- PCA.3D.brainstem$shapes$shapes.comp2$max
PC2min.brainstem <- PCA.3D.brainstem$shapes$shapes.comp2$min


spheres3d(PC1max.brainstem, col = "#D55E00", radius = 0.003)
draw_segments(links = brainstem_links, landmarks = PC1max.brainstem, lwd = 5)
spheres3d(PC1min.brainstem, col = "#D55E00", radius = 0.003)
draw_segments(links = brainstem_links, landmarks = PC1min.brainstem, lwd = 5)
spheres3d(PC2max.brainstem, col = "#D55E00", radius = 0.003)
draw_segments(links = brainstem_links, landmarks = PC2max.brainstem, lwd = 5)
spheres3d(PC2min.brainstem, col = "#D55E00", radius = 0.003)
draw_segments(links = brainstem_links, landmarks = PC2min.brainstem, lwd = 5)


PCA <- geomorph::gm.prcomp(Y.gpa.endo)
spheres3d(PCA$shapes$shapes.comp1$min[neocortex, ], col = "#009E73", radius = 0.004)
draw_segments(links = neocortex_links, landmarks = PCA$shapes$shapes.comp1$min[neocortex, ], lwd = 5)
spheres3d(PCA$shapes$shapes.comp1$min[olfactory, ], col = "#F0E442", radius = 0.004)
draw_segments(links = olfactory_links, landmarks = PCA$shapes$shapes.comp1$min[olfactory, ], lwd = 5)
spheres3d(PCA$shapes$shapes.comp1$min[cerebellum, ], col = "#56B4E9", radius = 0.004)
draw_segments(links = cerebellum_links, landmarks = PCA$shapes$shapes.comp1$min[cerebellum, ], lwd = 5)
spheres3d(PCA$shapes$shapes.comp1$min[brainstem, ], col = "#D55E00", radius = 0.004)
draw_segments(links = brainstem_links, landmarks = PCA$shapes$shapes.comp1$min[brainstem, ], lwd = 5)
spheres3d(PCA$shapes$shapes.comp1$min[brainstem, ], col = "#D55E00", radius = 0.004)

clear3d()

spheres3d(PCA$shapes$shapes.comp2$min[neocortex, ], col = "#009E73", radius = 0.004)
draw_segments(links = neocortex_links, landmarks = PCA$shapes$shapes.comp2$min[neocortex, ], lwd = 5)
spheres3d(PCA$shapes$shapes.comp2$min[olfactory, ], col = "#F0E442", radius = 0.004)
draw_segments(links = olfactory_links, landmarks = PCA$shapes$shapes.comp2$min[olfactory, ], lwd = 5)
spheres3d(PCA$shapes$shapes.comp2$min[cerebellum, ], col = "#56B4E9", radius = 0.004)
draw_segments(links = cerebellum_links, landmarks = PCA$shapes$shapes.comp2$min[cerebellum, ], lwd = 5)
spheres3d(PCA$shapes$shapes.comp2$min[brainstem, ], col = "#D55E00", radius = 0.004)
draw_segments(links = brainstem_links, landmarks = PCA$shapes$shapes.comp2$min[brainstem, ], lwd = 5)
spheres3d(PCA$shapes$shapes.comp2$min[brainstem, ], col = "#D55E00", radius = 0.004)
