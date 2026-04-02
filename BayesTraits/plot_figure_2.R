library(BayesTraitR)
library(ape)
library(viridis)
library(ggtree)
library(deeptime)
library(phytools)
library(viridis)
library(tidyverse)
library(patchwork)
library(rphylopic)
library(cowplot)
library(extrafont)



source("./scripts/utilities/BayesTraitsPlottingFunctions.R")

time_tree <- read.nexus("./phylogeny_construction/final_tree_hyp1.nex")


data <- read.csv("./input_data/mammal_data.csv", row.names = 1)
cladenames <- levels(as.factor(data$superorder))


clade_ancestral_nodes <- data.frame(superorder = cladenames, MRCA = NA)
for (i in 1:nrow(clade_ancestral_nodes)) {
  these_tips <- rownames(data)[which(data$superorder == clade_ancestral_nodes$superorder[i])]
  if (length(these_tips) > 1) {
    this_anc <- getMRCA(time_tree, these_tips)
  } else {
    this_anc <- NA
  }
  clade_ancestral_nodes$MRCA[i] <- this_anc
}
# get phylopics and associate with clades

phylopic_df <- data.frame(
  cladename = c(
    "Primates",
    "Whales",
    "Artiodactyla",
    "Perissodactyla",
    "Carnivora",
    "Chiroptera",
    "Eulipotyphla",
    "Muroidea",
    "Echimyidae",
    "Lagomorpha",
    "Pilosa",
    "Proboscidea",
    "Diprotodontia",
    "Monotremata",
    "Phocidae"
  ),
  target_node = c(
    getMRCA(time_tree, c("Homo_sapiens", "Gorilla_gorilla")),
    getMRCA(time_tree, c("Monodon_monoceros", "Protocetus_atavus")),
    getMRCA(time_tree, rownames(data)[which(data$order == "Artiodactyla")]),
    getMRCA(time_tree, rownames(data)[which(data$order == "Perissodactyla")]),
    getMRCA(time_tree, rownames(data)[which(data$order == "Carnivora")]),
    getMRCA(time_tree, rownames(data)[which(data$order == "Chiroptera")]),
    getMRCA(time_tree, rownames(data)[which(data$order == "Eulipotyphla")]),
    getMRCA(time_tree, c("Cricetus_cricetus", "Mus_musculus")),
    getMRCA(time_tree, c("Capromys_pilorides", "Hystrix_indica")),
    getMRCA(time_tree, rownames(data)[which(data$order == "Lagomorpha")]),
    getMRCA(time_tree, rownames(data)[which(data$order == "Pilosa")]),
    getMRCA(time_tree, rownames(data)[which(data$order == "Proboscidea")]),
    getMRCA(time_tree, rownames(data)[which(data$order == "Diprotodontia")]),
    getMRCA(time_tree, rownames(data)[which(data$order == "Monotremata")]),
    getMRCA(time_tree, c("Neomonachus_tropicalis", "Odobenus_rosmarus"))
  ),
  image = c(
    "9fae30cd-fb59-4a81-a39c-e1826a35f612",
    "f36c9daa-a102-42dd-88ac-a126753943d2",
    "f25267a2-90b4-4ca7-af95-a425f3564bbe",
    "81caf94e-5cbe-4e5e-8101-545abea2bfc0",
    "a8be3867-83d4-4614-9fed-9bcd2546dd5b",
    "21180755-3394-40bf-93eb-810954c0f7ba",
    "f83c6893-f0ed-4ec4-b558-aa774c5c9b5b",
    "f7d6d04c-73fa-4bf3-8c94-48134e6857b9",
    "5ee3b11f-9cd8-4aaf-a59e-ac1ddc18d5e8",
    "630a20e2-2d73-49ec-8cff-80f212b638c6",
    "96f7ce4a-633b-42a8-961d-c7a12c8684ce",
    "80db1004-bc9f-4318-84e7-cdd9639a1f3e",
    "0904270e-b105-46e0-b81f-c4911d47d467",
    "b406c409-2735-4a3d-a7aa-8afe0b6e72dc",
    "41d9f21e-7e0b-42d2-a088-e154bfa63984"
  ),
  superorder = c(
    "Euarchontoglires",
    "Laurasiatheria",
    "Laurasiatheria",
    "Laurasiatheria",
    "Laurasiatheria",
    "Laurasiatheria",
    "Laurasiatheria",
    "Euarchontoglires",
    "Euarchontoglires",
    "Euarchontoglires",
    "Xenarthra",
    "Afrotheria",
    "Marsupialia",
    "Monotremata",
    "Laurasiatheria"
  ),
  shape_size = c(
    23,
    8.5,
    22,
    14,
    14,
    11,
    12,
    16,
    18,
    22,
    18,
    20,
    20,
    18,
    10
  )
)
clade_ancestral_nodes <- na.omit(clade_ancestral_nodes)

primate_tips <- rownames(data)[which(data$order == "Primates")]
primateMRCA <- getMRCA(time_tree, primate_tips)
whippomorpha_tips <- c("Hippopotamus_amphibius", "Monodon_monoceros")
whippomorphaMRCA <- getMRCA(time_tree, whippomorpha_tips)
whippomorpha_tips2 <- extract.clade(time_tree, whippomorphaMRCA)$tip.label


col.gp.1 <- c(
  "Afrotheria" = "#e1d69c",
  "Euarchontoglires" = "#5fc0a4",
  "Laurasiatheria" = "#1e798a",
  "Leptictida" = "#13805d",
  "Marsupialia" = "#f2464e",
  "Monotremata" = "#7c0c17",
  "Multituberculata" = "#570407",
  "Sparassodonta" = "#fc951d",
  "Xenarthra" = "#efc35a"
)


rjpp_results.shape <- process_PPP("./BayesTraits/postprocessing/whole_brain/hyp_1/whole_brain_lambda_results.txt",
  phy = ladderize(time_tree),
  save_summary_trees = F,
  col.palette = "viridis"
)


data_with_shape_rates <- rjpp_results.shape$all.res[[1]]
data_with_shape_rates <- data_with_shape_rates %>%
  dplyr::rename(., node = child.node) %>%
  dplyr::rename(., Mean_Rate_Shape = Mean.SigV)
tree_with_shape_data <- treeio::full_join(time_tree, y = data_with_shape_rates, by = "node")

p_shape <- ggtree(tree_with_shape_data, layout = "rectangular", ladderize = T, lineend = "square", aes(color = log(Mean_Rate_Shape))) +
  # geom_tiplab(size=1,color='black')+
  scale_color_viridis_c(name = "log(rate of evolution\nof brain shape evolution)", option = "turbo", direction = 1) +
  geom_nodepoint(aes(subset = node %in% filter(tree_with_shape_data@extraInfo, Pct.time.scaled >= 50 & Mean.Scalar >= 2)$node),
    col = "black",
    shape = 18,
    size = 2.7
  )
p_shape <- revts(p_shape)
p_shape
primate_tip_numbers <- which(p_shape$data$label %in% primate_tips)
whippomorpha_tip_numbers <- which(p_shape$data$label %in% whippomorpha_tips2)
p_shape$data$label[which(p_shape$data$isTip == TRUE)] <- sub("_", " ", p_shape$data$label[which(p_shape$data$isTip == TRUE)])
# extract subclades before going further
primate_tree_plot <- viewClade(p_shape + xlim(-202, 15) + aes(size = .001), primateMRCA) +
  geom_tiplab(aes(subset = node %in% primate_tip_numbers), size = 1, color = "black", fontface = 3, align = T) +
  theme(legend.position = "none", plot.margin = unit(c(0.5, 1, 0.5, 0.7), "lines"))
whale_tree_plot <- viewClade(p_shape + xlim(-202, 15), whippomorphaMRCA) +
  geom_tiplab(aes(subset = node %in% whippomorpha_tip_numbers), size = 1, color = "black", fontface = 3, align = T) +
  aes(linewidth = .1) +
  theme(legend.position = "none", plot.margin = unit(c(0.5, 1, 0.5, 0.7), "lines"))
plot_grid(whale_tree_plot, primate_tree_plot)

p_shape <- p_shape + ggnewscale::new_scale_color()
p_shape <- p_shape +
  geom_cladelab(
    data = clade_ancestral_nodes,
    mapping = aes(
      node = MRCA,
      label = superorder,
      color = superorder
    ), textcolour = "black", offset = 3, offset.text = 22, barsize = 2, fontsize = 0, extend = .5, hjust = .5
  ) +
  # add node labels here:,p_shape@data Pct.time.scaled >= 50 & Mean.Scalar >= 2)
  scale_color_manual(values = col.gp.1, guide = "none") +
  theme(
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 8),
    plot.margin = unit(c(.5, 1, .7, 0.2), "lines"),
    legend.position = "none",
    text = element_text(family = "Arial")
  ) +
  # uncomment these to highlight key cades
  # geom_hilight(node=whippomorphaMRCA, fill=NA,color="red")+
  #  geom_hilight(node=primateMRCA, fill=NA,color="red")+
  xlim(-202, 7)
p_shape <- p_shape + coord_geo(neg = T, abbrv = T, size = 3, bord = c(), height = unit(1, "line"), ylim = c(-2, Ntip(time_tree)))
p_shape
# get node coordinates from the ggtree object
node_coords <- p_shape$data %>%
  filter(node %in% phylopic_df$target_node) %>%
  left_join(phylopic_df, by = c("node" = "target_node"))
# need to adjust the position of a few of the node pictures

node_coords <- node_coords %>%
  mutate(y_adjusted = y) %>%
  mutate(y_adjusted = ifelse(cladename %in% c("Lagomorpha", "Chiroptera"), y + 6, y_adjusted)) %>%
  mutate(y_adjusted = ifelse(cladename %in% c("Echimyidae", "Pilosa"), y + 14, y_adjusted)) %>%
  mutate(y_adjusted = ifelse(cladename %in% c("Proboscidea"), y - 10, y_adjusted)) %>%
  mutate(y_adjusted = ifelse(cladename %in% c("Monotremata"), y - 15, y_adjusted)) %>%
  mutate(y_adjusted = ifelse(cladename %in% c("Diprotodontia"), y - 20, y_adjusted))

# p_shape <- p_shape +
#  geom_phylopic(
#    data = node_coords,
#    aes(x = 24 , y = y, uuid = image,size=shape_size, color = superorder),
#    alpha = 1,
#  )+scale_size_identity()


p_shape_histo <- ggplot(
  tree_with_shape_data@extraInfo,
  aes(x = log(Mean_Rate_Shape), fill = after_stat(x))
) +
  geom_histogram(bins = 20, color = "black", linewidth = .2) +
  scale_fill_viridis_c(
    option = "turbo",
    direction = 1
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 9),
    panel.background = element_blank(),
    text = element_text(family = "Arial")
  ) +
  xlab("log(rate of\nbrain shape\nevolution)")
p_shape_composite <- ggdraw() +
  draw_plot(p_shape + annotate("text", x = -198, y = 460, label = "bold(A)", parse = TRUE, size = 8)) +
  draw_plot(p_shape_histo,
    x = 0.1,
    y = .65,
    width = .35,
    height = .2
  )
# plot the size rates just to get its dimensions

rjpp_results.size <- process_PPP("./BayesTraits/brainsize/linear_BM/data/brain_body_BM_processed.txt",
  phy = ladderize(time_tree),
  save_summary_trees = F,
  col.palette = "viridis"
)
data_with_size_rates <- rjpp_results.size$all.res[[1]]
data_with_size_rates <- data_with_size_rates %>%
  dplyr::rename(., node = child.node) %>%
  dplyr::rename(., Mean_Rate_Size = Mean.SigV)
tree_with_size_data <- treeio::full_join(time_tree, y = data_with_size_rates, by = "node")

plottemp <- ggtree(tree_with_size_data)
tree_with_data_just_data <- plottemp$data

## reverse x-axis and
## set offset to make the tree on the right-hand side of the first tree
tree_with_data_just_data$x <- max(tree_with_data_just_data$x) - tree_with_data_just_data$x + max(p_shape$data$x) + 1


data(periods)
periods2 <- periods
periods2$max_age <- periods2$max_age - max(tree_with_data_just_data$x)
periods2$min_age <- periods2$min_age - max(tree_with_data_just_data$x)
p_size_histo <- ggplot(tree_with_data_just_data, aes(x = log(Mean_Rate_Size), fill = after_stat(x))) +
  geom_histogram(bins = 20, color = "black", linewidth = .2) +
  scale_fill_viridis_c(
    option = "turbo",
    direction = 1
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 9)
  ) +
  xlab("log(rate of\nrelative brain size\nevolution")

p_size <- revts(ggtree(
  tree_with_data_just_data,
  aes(color = log(Mean_Rate_Size)),
  lineend = "square"
)) +
  ggnewscale::new_scale_fill() +
  scale_color_viridis_c(name = "log(rate of evolution\nof brain size scaling)", option = "turbo", direction = 1) +
  theme(
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 8),
    plot.margin = unit(c(.5, 1, .7, 0.2), "lines"),
    legend.position = "none",
    text = element_text(family = "Arial")
  ) + ggnewscale::new_scale_color() +
  geom_cladelab(
    data = clade_ancestral_nodes,
    mapping = aes(
      node = MRCA,
      label = superorder,
      color = superorder
    ), textcolour = "black", offset = -max(tree_with_data_just_data$x) - 3, offset.text = -5, barsize = 2, extend = .5, fontsize = 0, align = T
  ) +
  scale_color_manual(values = col.gp.1, guide = "none") +
  xlim(-202, 7) +
  coord_geo(dat = periods2, neg = F, size = 3, bord = c(), height = unit(1, "line"), abbrv = T, ylim = c(-2, Ntip(time_tree)))


p_size_composite <- ggdraw() +
  draw_plot(p_size + annotate("text", x = -4, y = 460, label = "bold(B)", parse = TRUE, size = 8)) +
  draw_plot(p_size_histo,
    x = .47,
    y = .65,
    width = .35,
    height = .2
  )

p_strip <- revts(ggtree(time_tree, layout = "rectangular", color = "transparent", ladderize = T)) +
  geom_cladelab(
    data = clade_ancestral_nodes,
    mapping = aes(
      node = MRCA,
      label = superorder,
      color = superorder
    ),
    textcolour = "black", offset = 0, offset.text = -2,
    barsize = 0, barcolor = "transparent", fontsize = 2.5, extend = 0.5, hjust = 0.5
  ) +
  geom_phylopic(
    data = node_coords,
    aes(x = -2, y = y_adjusted, uuid = image, color = superorder, height = shape_size),
    alpha = 1, na.rm = TRUE
  ) +
  scale_height_continuous(range = c(10, 22)) +
  scale_color_manual(values = col.gp.1, guide = "none") +
  scale_fill_manual(values = col.gp.1, guide = "none") +
  theme_void() +
  theme(legend.position = "none", plot.margin = unit(c(.5, 0, .7, 0), "lines")) +
  xlim(-10, 10) +
  coord_geo(neg = T, abbrv = T, size = 0, bord = c(), height = unit(0, "line"), ylim = c(-2, Ntip(time_tree)))

big_trees_composite <- plot_grid(
  p_shape + annotate("text", x = -192, y = 460, label = "bold(A)", parse = TRUE, size = 8),
  p_strip,
  p_size + annotate("text", x = 0, y = 460, label = "bold(B)", parse = TRUE, size = 8),
  nrow = 1,
  rel_widths = c(10, 2, 10),
  align = "h",
  axis = "tb"
)

legend_plot <- ggplot() +
  geom_point(
    data = data.frame(x = 1, y = 1),
    aes(x = x, y = y),
    shape = 18, size = 2.7, color = "black"
  ) +
  annotate("text",
    x = 1.1, y = 1, label = "significant rate shift",
    hjust = 0, size = 3, family = "Arial"
  ) +
  theme_void() +
  xlim(0.9, 3) +
  ylim(0.5, 1.5)

big_trees_composite <- ggdraw() +
  draw_plot(big_trees_composite) +
  draw_plot(p_size_histo,
    x = .75,
    y = .65,
    width = .2,
    height = .2
  ) +
  draw_plot(p_shape_histo,
    x = .02,
    y = .65,
    width = .2,
    height = .2
  ) +
  draw_plot(legend_plot, x = 0.05, y = 0.60, width = 0.2, height = 0.05)


big_trees_composite
p_shape2 <- ggtree(tree_with_shape_data, layout = "rectangular", ladderize = T, lineend = "square", linewidth = .7, aes(color = log(Mean_Rate_Shape))) +
  # geom_tiplab(size=1,color='black')+
  scale_color_viridis_c(name = "log(rate of evolution\nof brain shape evolution)", option = "turbo", direction = 1) +
  geom_nodepoint(aes(subset = node %in% filter(tree_with_shape_data@extraInfo, Pct.time.scaled >= 50 & Mean.Scalar >= 2)$node),
    fill = "black",
    color = "white",
    shape = 23,
    size = 2.7
  )
p_shape2 <- revts(p_shape2)
p_shape2
primate_tip_numbers <- which(p_shape2$data$label %in% primate_tips)
whippomorpha_tip_numbers <- which(p_shape2$data$label %in% whippomorpha_tips2)
p_shape2$data$label[which(p_shape2$data$isTip == TRUE)] <- sub("_", " ", p_shape2$data$label[which(p_shape2$data$isTip == TRUE)])
# extract subclades before going further
primate_tree_plot <- viewClade(p_shape2 + xlim(-202, 15), primateMRCA) +
  geom_tiplab(aes(subset = node %in% primate_tip_numbers), size = 1.6, color = "black", fontface = 3, align = T) +
  theme(legend.position = "none", plot.margin = unit(c(0.5, 1, 0.5, 1), "lines"))
whale_tree_plot <- viewClade(p_shape2 + xlim(-202, 15), whippomorphaMRCA) +
  geom_tiplab(aes(subset = node %in% whippomorpha_tip_numbers), size = 2, color = "black", fontface = 3, align = T) +
  aes(linewidth = .1) +
  theme(legend.position = "none", plot.margin = unit(c(0.5, 1, 0.5, 1), "lines"))

sub_trees_composite <- plot_grid(whale_tree_plot, primate_tree_plot)

super_composite <- plot_grid(big_trees_composite, sub_trees_composite, ncol = 1, nrow = 2, rel_heights = c(1, .3))

cairo_pdf("~/Downloads/primatewhales.pdf", width = 7, height = 10)
super_composite
dev.off()

cairo_pdf("./Images/Figures/Figure_2.pdf", width = 7, height = 10)
big_trees_composite
dev.off()
