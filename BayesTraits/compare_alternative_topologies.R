library(BayesTraitR)
library(ape)
library(viridis)
library(ggtree)
library(deeptime)
library(phytools)
library(viridis)
library(tidyverse)
library(patchwork)
library(showtext)

source("./scripts/utilities/BayesTraitsPlottingFunctions.R")

time_tree_1 <- read.nexus("./phylogeny_construction/final_tree_hyp1.nex")


rjpp_results.wholebrain_1 <- process_PPP("./BayesTraits/postprocessing/whole_brain/hyp_1/whole_brain_lambda_results.txt",
  phy = ladderize(time_tree_1),
  save_summary_trees = F,
  col.palette = "viridis"
)
time_tree_2 <- read.nexus("./phylogeny_construction/final_tree_hyp2.nex")


rjpp_results.wholebrain_2 <- process_PPP("./BayesTraits/postprocessing/whole_brain/hyp_2/whole_brain_lambda_results.txt",
  phy = ladderize(time_tree_2),
  save_summary_trees = F,
  col.palette = "viridis"
)
time_tree_3 <- read.nexus("./phylogeny_construction/final_tree_hyp3.nex")


rjpp_results.wholebrain_3 <- process_PPP("./BayesTraits/postprocessing/whole_brain/hyp_3/whole_brain_lambda_results.txt",
  phy = ladderize(time_tree_3),
  save_summary_trees = F,
  col.palette = "viridis"
)


time_tree_4 <- read.nexus("./phylogeny_construction/final_tree_hyp4.nex")


rjpp_results.wholebrain_4 <- process_PPP("./BayesTraits/postprocessing/whole_brain/hyp_4/whole_brain_lambda_results.txt",
  phy = ladderize(time_tree_4),
  save_summary_trees = F,
  col.palette = "viridis"
)

time_tree_5 <- read.nexus("./phylogeny_construction/final_tree_hyp5.nex")


rjpp_results.wholebrain_5 <- process_PPP("./BayesTraits/postprocessing/whole_brain/hyp_5/whole_brain_lambda_results.txt",
  phy = ladderize(time_tree_5),
  save_summary_trees = F,
  col.palette = "viridis"
)
time_tree_6 <- read.nexus("./phylogeny_construction/final_tree_hyp6.nex")


rjpp_results.wholebrain_6 <- process_PPP("./BayesTraits/postprocessing/whole_brain/hyp_6/whole_brain_lambda_results.txt",
  phy = ladderize(time_tree_6),
  save_summary_trees = F,
  col.palette = "viridis"
)

time_tree_7 <- read.nexus("./phylogeny_construction/final_tree_hyp7.nex")


rjpp_results.wholebrain_7 <- process_PPP("./BayesTraits/postprocessing/whole_brain/hyp_7/whole_brain_lambda_results.txt",
  phy = ladderize(time_tree_7),
  save_summary_trees = F,
  col.palette = "viridis"
)

time_tree_8 <- read.nexus("./phylogeny_construction/final_tree_hyp8.nex")


rjpp_results.wholebrain_8 <- process_PPP("./BayesTraits/postprocessing/whole_brain/hyp_8/whole_brain_lambda_results.txt",
  phy = ladderize(time_tree_8),
  save_summary_trees = F,
  col.palette = "viridis"
)


# get a relative rate at time for each dataset
RaT.wholebrain1 <- rate.at.time.df(timeslices = 0.1, obj = rjpp_results.wholebrain_1$all.res[[1]], plot = F, relative.rates = "mean")
RaT.wholebrain2 <- rate.at.time.df(timeslices = 0.1, obj = rjpp_results.wholebrain_2$all.res[[1]], plot = F, relative.rates = "mean")
RaT.wholebrain3 <- rate.at.time.df(timeslices = 0.1, obj = rjpp_results.wholebrain_3$all.res[[1]], plot = F, relative.rates = "mean")
RaT.wholebrain4 <- rate.at.time.df(timeslices = 0.1, obj = rjpp_results.wholebrain_4$all.res[[1]], plot = F, relative.rates = "mean")
RaT.wholebrain5 <- rate.at.time.df(timeslices = 0.1, obj = rjpp_results.wholebrain_5$all.res[[1]], plot = F, relative.rates = "mean")
RaT.wholebrain6 <- rate.at.time.df(timeslices = 0.1, obj = rjpp_results.wholebrain_6$all.res[[1]], plot = F, relative.rates = "mean")
RaT.wholebrain7 <- rate.at.time.df(timeslices = 0.1, obj = rjpp_results.wholebrain_7$all.res[[1]], plot = F, relative.rates = "mean")
RaT.wholebrain8 <- rate.at.time.df(timeslices = 0.1, obj = rjpp_results.wholebrain_8$all.res[[1]], plot = F, relative.rates = "mean")

# extract the rates through time


rate.mean.wholebrain1 <- extract.stat(RaT.wholebrain1, stat = "mean", plot = F, range = "confidence")
rate.mean.wholebrain1 <- rate.mean.wholebrain1 %>% mutate(tree = "Hypothesis 1")
rate.mean.wholebrain2 <- extract.stat(RaT.wholebrain2, stat = "mean", plot = F, range = "confidence")
rate.mean.wholebrain2 <- rate.mean.wholebrain2 %>% mutate(tree = "Hypothesis 2")
rate.mean.wholebrain3 <- extract.stat(RaT.wholebrain3, stat = "mean", plot = F, range = "confidence")
rate.mean.wholebrain3 <- rate.mean.wholebrain3 %>% mutate(tree = "Hypothesis 3")
rate.mean.wholebrain4 <- extract.stat(RaT.wholebrain4, stat = "mean", plot = F, range = "confidence")
rate.mean.wholebrain4 <- rate.mean.wholebrain4 %>% mutate(tree = "Hypothesis 4")
rate.mean.wholebrain5 <- extract.stat(RaT.wholebrain5, stat = "mean", plot = F, range = "confidence")
rate.mean.wholebrain5 <- rate.mean.wholebrain5 %>% mutate(tree = "Hypothesis 5")
rate.mean.wholebrain6 <- extract.stat(RaT.wholebrain6, stat = "mean", plot = F, range = "confidence")
rate.mean.wholebrain6 <- rate.mean.wholebrain6 %>% mutate(tree = "Hypothesis 6")
rate.mean.wholebrain7 <- extract.stat(RaT.wholebrain7, stat = "mean", plot = F, range = "confidence")
rate.mean.wholebrain7 <- rate.mean.wholebrain7 %>% mutate(tree = "Hypothesis 7")
rate.mean.wholebrain8 <- extract.stat(RaT.wholebrain8, stat = "mean", plot = F, range = "confidence")
rate.mean.wholebrain8 <- rate.mean.wholebrain8 %>% mutate(tree = "Hypothesis 8")


rates_through_time <- bind_rows(
  rate.mean.wholebrain1,
  rate.mean.wholebrain2,
  rate.mean.wholebrain3,
  rate.mean.wholebrain4,
  rate.mean.wholebrain5,
  rate.mean.wholebrain6,
  rate.mean.wholebrain7,
  rate.mean.wholebrain8
)

annotation <- data.frame(
  x = c(52, 72),
  y = c(3, 3),
  label = c("K-PG/nBoundary", "PETM")
)


region_colors <- palette.colors(n = 8, palette = "Tableau 10")


rate_through_time_plot <- ggplot(rates_through_time, aes(x = -1 * (time), color = tree, fill = tree)) +
  scale_x_continuous(
    trans = scales::reverse_trans(),
    breaks = seq(100, 0, by = -10),
    limits = c(100, 0)
  ) +
  # uncomment for confidence intervals
  # geom_ribbon(
  #   aes(ymin = `5%`, ymax = `95%`),
  #   alpha = 0.3
  # ) +
  geom_line(
    aes(y = rate)
  ) +
  theme_minimal() +
  ylim(0, 4) +
  scale_color_manual(
    name = "",
    values = region_colors
  ) +
  geom_vline(xintercept = 66) +
  geom_vline(xintercept = 55.8, linetype = "dashed") +
  coord_geo(center_end_labels = T, height = unit(1, "lines"), size = 4) +
  labs(x = "Ma before present", y = "Relative Rate of Phenotypic Evolution") +
  annotate(
    geom = "text", label = "PETM",
    x = 51, y = 3.3, size = 4.5
  ) +
  annotate(
    geom = "text", label = "K-Pg\nMass\nExtinction",
    x = 74, y = 3.3, size = 4.5
  ) +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.8, 0.75),
    legend.key.width = unit(1, "cm")
  )


cairo_pdf("./Images/Figures/Figure_rates_through_times_supp.pdf", width = 7, height = 4)
rate_through_time_plot
dev.off()
############
# plot trees
############
# first build clade labels for the tree plot

data <- read.csv("./input_data/mammal_data.csv", row.names = 1)
cladenames <- levels(as.factor(data$order))

clade_ancestral_nodes <- data.frame(order = cladenames, MRCA = NA)
for (i in 1:nrow(clade_ancestral_nodes)) {
  these_tips <- rownames(data)[which(data$order == clade_ancestral_nodes$order[i])]
  if (length(these_tips) > 1) {
    this_anc <- getMRCA(time_tree_1, these_tips)
  } else {
    this_anc <- NA
  }
  clade_ancestral_nodes$MRCA[i] <- this_anc
}

grey_colors.1 <- c(
  "Artiodactyla" = "#3C3C3C",
  "Perissodactyla" = "#717171",
  "Notoungulata" = "#3C3C3C",
  "Pantodonta" = "#717171",
  "Carnivora" = "#3C3C3C",
  "Hyaenodonta" = "#717171",
  "Pholidota" = "#3C3C3C",
  "Chiroptera" = "#717171",
  "Eulipotyphla" = "#3C3C3C",
  "Rodentia" = "#717171",
  "Lagomorpha" = "#3C3C3C",
  "Scandentia" = "#717171",
  "Primates" = "#3C3C3C",
  "Pilosa" = "#717171",
  "Cingulata" = "#3C3C3C",
  "Afrosoricida" = "#717171",
  "Sirenia" = "#3C3C3C",
  "Proboscidea" = "#717171",
  "Diprotodontia" = "#3C3C3C",
  "Dasyuromorphia" = "#717171",
  "Peramelemorphia" = "#3C3C3C",
  "Didelphimorphia" = "#717171",
  "Monotremata" = "#3C3C3C"
)

clade_ancestral_nodes <- na.omit(clade_ancestral_nodes)
# remove clades we dont want to plot
clade_ancestral_nodes <- clade_ancestral_nodes[-which(clade_ancestral_nodes$order == "Plesiadapiformes"), ]


###########
###########

rate_data_wholebrain_1 <- rjpp_results.wholebrain_1$all.res[[1]]
rate_data_wholebrain_1 <- rate_data_wholebrain_1 %>%
  dplyr::rename(., node = child.node) %>%
  dplyr::rename(., Mean_Rate = Mean.SigV)
tree_with_wholebrain_1 <- treeio::full_join(time_tree_1, y = rate_data_wholebrain_1, by = "node")

rate_data_wholebrain_2 <- rjpp_results.wholebrain_2$all.res[[1]]
rate_data_wholebrain_2 <- rate_data_wholebrain_2 %>%
  dplyr::rename(., node = child.node) %>%
  dplyr::rename(., Mean_Rate = Mean.SigV)
tree_with_wholebrain_2 <- treeio::full_join(time_tree_2, y = rate_data_wholebrain_2, by = "node")

rate_data_wholebrain_3 <- rjpp_results.wholebrain_3$all.res[[1]]
rate_data_wholebrain_3 <- rate_data_wholebrain_3 %>%
  dplyr::rename(., node = child.node) %>%
  dplyr::rename(., Mean_Rate = Mean.SigV)
tree_with_wholebrain_3 <- treeio::full_join(time_tree_3, y = rate_data_wholebrain_3, by = "node")

rate_data_wholebrain_4 <- rjpp_results.wholebrain_4$all.res[[1]]
rate_data_wholebrain_4 <- rate_data_wholebrain_4 %>%
  dplyr::rename(., node = child.node) %>%
  dplyr::rename(., Mean_Rate = Mean.SigV)
tree_with_wholebrain_4 <- treeio::full_join(time_tree_4, y = rate_data_wholebrain_4, by = "node")

rate_data_wholebrain_5 <- rjpp_results.wholebrain_5$all.res[[1]]
rate_data_wholebrain_5 <- rate_data_wholebrain_5 %>%
  dplyr::rename(., node = child.node) %>%
  dplyr::rename(., Mean_Rate = Mean.SigV)
tree_with_wholebrain_5 <- treeio::full_join(time_tree_5, y = rate_data_wholebrain_5, by = "node")

rate_data_wholebrain_6 <- rjpp_results.wholebrain_6$all.res[[1]]
rate_data_wholebrain_6 <- rate_data_wholebrain_6 %>%
  dplyr::rename(., node = child.node) %>%
  dplyr::rename(., Mean_Rate = Mean.SigV)
tree_with_wholebrain_6 <- treeio::full_join(time_tree_6, y = rate_data_wholebrain_6, by = "node")

rate_data_wholebrain_7 <- rjpp_results.wholebrain_7$all.res[[1]]
rate_data_wholebrain_7 <- rate_data_wholebrain_7 %>%
  dplyr::rename(., node = child.node) %>%
  dplyr::rename(., Mean_Rate = Mean.SigV)
tree_with_wholebrain_7 <- treeio::full_join(time_tree_7, y = rate_data_wholebrain_7, by = "node")

rate_data_wholebrain_8 <- rjpp_results.wholebrain_8$all.res[[1]]
rate_data_wholebrain_8 <- rate_data_wholebrain_8 %>%
  dplyr::rename(., node = child.node) %>%
  dplyr::rename(., Mean_Rate = Mean.SigV)
tree_with_wholebrain_8 <- treeio::full_join(time_tree_8, y = rate_data_wholebrain_8, by = "node")


p_wholebrain_1 <- ggtree(tree_with_wholebrain_1,
  layout = "rectangular",
  ladderize = T,
  lineend = "square", linewidth = 0.3, aes(color = scale(log(Mean_Rate), scale = F))
) +
  scale_color_viridis_c(
    name = "log(relative evolutionary rate)", option = "turbo", direction = 1,
    breaks = function(x) c(min(x), max(x)),
    labels = c("Low", "High")
  )
p_wholebrain_1 <- revts(p_wholebrain_1)
p_wholebrain_1 <- p_wholebrain_1 + ggnewscale::new_scale_color()
p_wholebrain_1 <- p_wholebrain_1 +
  # geom_cladelab(
  #   data = clade_ancestral_nodes,
  #   mapping = aes(
  #     node = MRCA,
  #     label = order,
  #     color = order
  #   ), textcolour = "black", offset = 3, offset.text = 3, barsize=2,extend=.5,angle=65,fontsize = 2.5, align=T)+
  geom_nodepoint(aes(subset = node %in% filter(tree_with_wholebrain_1@extraInfo, Pct.time.scaled >= 50 & Mean.Scalar >= 2)$node),
    col = "black",
    shape = 18,
    size = 3
  ) +
  scale_color_manual(values = grey_colors.1, guide = "none") +
  theme(
    legend.text = element_text(family = "Arial"),
    legend.title = element_text(size = 8, hjust = 0.5, family = "Arial"),
    plot.margin = unit(c(0, .3, 0.3, 0.3), "lines"),
    legend.position = "inside",
    legend.position.inside = c(0.3, 0.7),
    legend.direction = "horizontal",
    legend.key.height = unit(.7, "lines"),
    legend.title.position = "top",
    legend.background = element_blank()
  ) +
  geom_point(shape = 18, size = 2.7, color = "black", x = -152, y = 300) +
  annotate("text",
    x = -150, y = 300, label = "significant rate shift",
    hjust = 0, size = 3, family = "Arial"
  ) +
  xlim(-202, 20)
p_wholebrain_1
p_wholebrain_1$data$label[which(p_wholebrain_1$data$isTip == TRUE)] <- sub("_", " ", p_wholebrain_1$data$label[which(p_wholebrain_1$data$isTip == TRUE)])
p_wholebrain_1 <- p_wholebrain_1 + geom_tiplab(align = T, fontface = 3, size = .7, linesize = .2)

p_wholebrain_2 <- ggtree(tree_with_wholebrain_2,
  layout = "rectangular",
  ladderize = T,
  lineend = "square", linewidth = 0.3, aes(color = scale(log(Mean_Rate), scale = F))
) +
  scale_color_viridis_c(
    name = "log(relative evolutionary rate)", option = "turbo", direction = 1,
    breaks = function(x) c(min(x), max(x)),
    labels = c("Low", "High")
  )
p_wholebrain_2 <- revts(p_wholebrain_2)
p_wholebrain_2 <- p_wholebrain_2 + ggnewscale::new_scale_color()
p_wholebrain_2 <- p_wholebrain_2 +
  # geom_cladelab(
  #   data = clade_ancestral_nodes,
  #   mapping = aes(
  #     node = MRCA,
  #     label = order,
  #     color = order
  #   ), textcolour = "black", offset = 3, offset.text = 3, barsize=2,extend=.5,angle=65,fontsize = 2.5, align=T)+
  geom_nodepoint(aes(subset = node %in% filter(tree_with_wholebrain_2@extraInfo, Pct.time.scaled >= 50 & Mean.Scalar >= 2)$node),
    col = "black",
    shape = 18,
    size = 1
  ) +
  scale_color_manual(values = grey_colors.1, guide = "none") +
  theme(
    legend.text = element_text(family = "Arial"),
    legend.title = element_text(size = 8, hjust = 0.5, family = "Arial"),
    plot.margin = unit(c(0, .3, 0.3, 0.3), "lines"),
    legend.position = "inside",
    legend.position.inside = c(0.3, 0.7),
    legend.direction = "horizontal",
    legend.key.height = unit(.7, "lines"),
    legend.title.position = "top",
    legend.background = element_blank()
  ) +
  geom_point(shape = 18, size = 2.7, color = "black", x = -152, y = 300) +
  annotate("text",
    x = -150, y = 300, label = "significant rate shift",
    hjust = 0, size = 3, family = "Arial"
  ) +
  xlim(-202, 20)
p_wholebrain_2
p_wholebrain_2$data$label[which(p_wholebrain_2$data$isTip == TRUE)] <- sub("_", " ", p_wholebrain_2$data$label[which(p_wholebrain_2$data$isTip == TRUE)])
p_wholebrain_2 <- p_wholebrain_2 + geom_tiplab(align = T, fontface = 3, size = .7, linesize = .2)

p_wholebrain_3 <- ggtree(tree_with_wholebrain_3,
  layout = "rectangular",
  ladderize = T,
  lineend = "square", linewidth = 0.3, aes(color = scale(log(Mean_Rate), scale = F))
) +
  scale_color_viridis_c(
    name = "log(relative evolutionary rate)", option = "turbo", direction = 1,
    breaks = function(x) c(min(x), max(x)),
    labels = c("Low", "High")
  )
p_wholebrain_3 <- revts(p_wholebrain_3)
p_wholebrain_3 <- p_wholebrain_3 + ggnewscale::new_scale_color()
p_wholebrain_3 <- p_wholebrain_3 +
  # geom_cladelab(
  #   data = clade_ancestral_nodes,
  #   mapping = aes(
  #     node = MRCA,
  #     label = order,
  #     color = order
  #   ), textcolour = "black", offset = 3, offset.text = 3, barsize=2,extend=.5,angle=65,fontsize = 2.5, align=T)+
  geom_nodepoint(aes(subset = node %in% filter(tree_with_wholebrain_3@extraInfo, Pct.time.scaled >= 50 & Mean.Scalar >= 2)$node),
    col = "black",
    shape = 18,
    size = 1
  ) +
  scale_color_manual(values = grey_colors.1, guide = "none") +
  theme(
    legend.text = element_text(family = "Arial"),
    legend.title = element_text(size = 8, hjust = 0.5, family = "Arial"),
    plot.margin = unit(c(0, .3, 0.3, 0.3), "lines"),
    legend.position = "inside",
    legend.position.inside = c(0.3, 0.7),
    legend.direction = "horizontal",
    legend.key.height = unit(.7, "lines"),
    legend.title.position = "top",
    legend.background = element_blank()
  ) +
  geom_point(shape = 18, size = 2.7, color = "black", x = -152, y = 300) +
  annotate("text",
    x = -150, y = 300, label = "significant rate shift",
    hjust = 0, size = 3, family = "Arial"
  ) +
  xlim(-202, 20)
p_wholebrain_3
p_wholebrain_3$data$label[which(p_wholebrain_3$data$isTip == TRUE)] <- sub("_", " ", p_wholebrain_3$data$label[which(p_wholebrain_3$data$isTip == TRUE)])
p_wholebrain_3 <- p_wholebrain_3 + geom_tiplab(align = T, fontface = 3, size = .7, linesize = .2)

p_wholebrain_4 <- ggtree(tree_with_wholebrain_4,
  layout = "rectangular",
  ladderize = T,
  lineend = "square", linewidth = 0.3, aes(color = scale(log(Mean_Rate), scale = F))
) +
  scale_color_viridis_c(
    name = "log(relative evolutionary rate)", option = "turbo", direction = 1,
    breaks = function(x) c(min(x), max(x)),
    labels = c("Low", "High")
  )
p_wholebrain_4 <- revts(p_wholebrain_4)
p_wholebrain_4 <- p_wholebrain_4 + ggnewscale::new_scale_color()
p_wholebrain_4 <- p_wholebrain_4 +
  # geom_cladelab(
  #   data = clade_ancestral_nodes,
  #   mapping = aes(
  #     node = MRCA,
  #     label = order,
  #     color = order
  #   ), textcolour = "black", offset = 3, offset.text = 3, barsize=2,extend=.5,angle=65,fontsize = 2.5, align=T)+
  geom_nodepoint(aes(subset = node %in% filter(tree_with_wholebrain_4@extraInfo, Pct.time.scaled >= 50 & Mean.Scalar >= 2)$node),
    col = "black",
    shape = 18,
    size = 1
  ) +
  scale_color_manual(values = grey_colors.1, guide = "none") +
  theme(
    legend.text = element_text(family = "Arial"),
    legend.title = element_text(size = 8, hjust = 0.5, family = "Arial"),
    plot.margin = unit(c(0, .3, 0.3, 0.3), "lines"),
    legend.position = "inside",
    legend.position.inside = c(0.3, 0.7),
    legend.direction = "horizontal",
    legend.key.height = unit(.7, "lines"),
    legend.title.position = "top",
    legend.background = element_blank()
  ) +
  geom_point(shape = 18, size = 2.7, color = "black", x = -152, y = 300) +
  annotate("text",
    x = -150, y = 300, label = "significant rate shift",
    hjust = 0, size = 3, family = "Arial"
  ) +
  xlim(-202, 20)
p_wholebrain_4
p_wholebrain_4$data$label[which(p_wholebrain_4$data$isTip == TRUE)] <- sub("_", " ", p_wholebrain_4$data$label[which(p_wholebrain_4$data$isTip == TRUE)])
p_wholebrain_4 <- p_wholebrain_4 + geom_tiplab(align = T, fontface = 3, size = .7, linesize = .2)

p_wholebrain_5 <- ggtree(tree_with_wholebrain_5,
  layout = "rectangular",
  ladderize = T,
  lineend = "square", linewidth = 0.3, aes(color = scale(log(Mean_Rate), scale = F))
) +
  scale_color_viridis_c(
    name = "log(relative evolutionary rate)", option = "turbo", direction = 1,
    breaks = function(x) c(min(x), max(x)),
    labels = c("Low", "High")
  )
p_wholebrain_5 <- revts(p_wholebrain_5)
p_wholebrain_5 <- p_wholebrain_5 + ggnewscale::new_scale_color()
p_wholebrain_5 <- p_wholebrain_5 +
  # geom_cladelab(
  #   data = clade_ancestral_nodes,
  #   mapping = aes(
  #     node = MRCA,
  #     label = order,
  #     color = order
  #   ), textcolour = "black", offset = 3, offset.text = 3, barsize=2,extend=.5,angle=65,fontsize = 2.5, align=T)+
  geom_nodepoint(aes(subset = node %in% filter(tree_with_wholebrain_5@extraInfo, Pct.time.scaled >= 50 & Mean.Scalar >= 2)$node),
    col = "black",
    shape = 18,
    size = 1
  ) +
  scale_color_manual(values = grey_colors.1, guide = "none") +
  theme(
    legend.text = element_text(family = "Arial"),
    legend.title = element_text(size = 8, hjust = 0.5, family = "Arial"),
    plot.margin = unit(c(0, .3, 0.3, 0.3), "lines"),
    legend.position = "inside",
    legend.position.inside = c(0.3, 0.7),
    legend.direction = "horizontal",
    legend.key.height = unit(.7, "lines"),
    legend.title.position = "top",
    legend.background = element_blank()
  ) +
  geom_point(shape = 18, size = 2.7, color = "black", x = -152, y = 300) +
  annotate("text",
    x = -150, y = 300, label = "significant rate shift",
    hjust = 0, size = 3, family = "Arial"
  ) +
  xlim(-202, 20)
p_wholebrain_5
p_wholebrain_5$data$label[which(p_wholebrain_5$data$isTip == TRUE)] <- sub("_", " ", p_wholebrain_5$data$label[which(p_wholebrain_5$data$isTip == TRUE)])
p_wholebrain_5 <- p_wholebrain_5 + geom_tiplab(align = T, fontface = 3, size = .7, linesize = .2)

p_wholebrain_6 <- ggtree(tree_with_wholebrain_6,
  layout = "rectangular",
  ladderize = T,
  lineend = "square", linewidth = 0.3, aes(color = scale(log(Mean_Rate), scale = F))
) +
  scale_color_viridis_c(
    name = "log(relative evolutionary rate)", option = "turbo", direction = 1,
    breaks = function(x) c(min(x), max(x)),
    labels = c("Low", "High")
  )
p_wholebrain_6 <- revts(p_wholebrain_6)
p_wholebrain_6 <- p_wholebrain_6 + ggnewscale::new_scale_color()
p_wholebrain_6 <- p_wholebrain_6 +
  # geom_cladelab(
  #   data = clade_ancestral_nodes,
  #   mapping = aes(
  #     node = MRCA,
  #     label = order,
  #     color = order
  #   ), textcolour = "black", offset = 3, offset.text = 3, barsize=2,extend=.5,angle=65,fontsize = 2.5, align=T)+
  geom_nodepoint(aes(subset = node %in% filter(tree_with_wholebrain_6@extraInfo, Pct.time.scaled >= 50 & Mean.Scalar >= 2)$node),
    col = "black",
    shape = 18,
    size = 1
  ) +
  scale_color_manual(values = grey_colors.1, guide = "none") +
  theme(
    legend.text = element_text(family = "Arial"),
    legend.title = element_text(size = 8, hjust = 0.5, family = "Arial"),
    plot.margin = unit(c(0, .3, 0.3, 0.3), "lines"),
    legend.position = "inside",
    legend.position.inside = c(0.3, 0.7),
    legend.direction = "horizontal",
    legend.key.height = unit(.7, "lines"),
    legend.title.position = "top",
    legend.background = element_blank()
  ) +
  geom_point(shape = 18, size = 2.7, color = "black", x = -152, y = 300) +
  annotate("text",
    x = -150, y = 300, label = "significant rate shift",
    hjust = 0, size = 3, family = "Arial"
  ) +
  xlim(-202, 20)
p_wholebrain_6
p_wholebrain_6$data$label[which(p_wholebrain_6$data$isTip == TRUE)] <- sub("_", " ", p_wholebrain_6$data$label[which(p_wholebrain_6$data$isTip == TRUE)])
p_wholebrain_6 <- p_wholebrain_6 + geom_tiplab(align = T, fontface = 3, size = .7, linesize = .2)

p_wholebrain_7 <- ggtree(tree_with_wholebrain_7,
  layout = "rectangular",
  ladderize = T,
  lineend = "square", linewidth = 0.3, aes(color = scale(log(Mean_Rate), scale = F))
) +
  scale_color_viridis_c(
    name = "log(relative evolutionary rate)", option = "turbo", direction = 1,
    breaks = function(x) c(min(x), max(x)),
    labels = c("Low", "High")
  )
p_wholebrain_7 <- revts(p_wholebrain_7)
p_wholebrain_7 <- p_wholebrain_7 + ggnewscale::new_scale_color()
p_wholebrain_7 <- p_wholebrain_7 +
  # geom_cladelab(
  #   data = clade_ancestral_nodes,
  #   mapping = aes(
  #     node = MRCA,
  #     label = order,
  #     color = order
  #   ), textcolour = "black", offset = 3, offset.text = 3, barsize=2,extend=.5,angle=65,fontsize = 2.5, align=T)+
  geom_nodepoint(aes(subset = node %in% filter(tree_with_wholebrain_7@extraInfo, Pct.time.scaled >= 50 & Mean.Scalar >= 2)$node),
    col = "black",
    shape = 18,
    size = 1
  ) +
  scale_color_manual(values = grey_colors.1, guide = "none") +
  theme(
    legend.text = element_text(family = "Arial"),
    legend.title = element_text(size = 8, hjust = 0.5, family = "Arial"),
    plot.margin = unit(c(0, .3, 0.3, 0.3), "lines"),
    legend.position = "inside",
    legend.position.inside = c(0.3, 0.7),
    legend.direction = "horizontal",
    legend.key.height = unit(.7, "lines"),
    legend.title.position = "top",
    legend.background = element_blank()
  ) +
  geom_point(shape = 18, size = 2.7, color = "black", x = -152, y = 300) +
  annotate("text",
    x = -150, y = 300, label = "significant rate shift",
    hjust = 0, size = 3, family = "Arial"
  ) +
  xlim(-202, 20)
p_wholebrain_7
p_wholebrain_7$data$label[which(p_wholebrain_7$data$isTip == TRUE)] <- sub("_", " ", p_wholebrain_7$data$label[which(p_wholebrain_7$data$isTip == TRUE)])
p_wholebrain_7 <- p_wholebrain_7 + geom_tiplab(align = T, fontface = 3, size = .7, linesize = .2)

p_wholebrain_8 <- ggtree(tree_with_wholebrain_8,
  layout = "rectangular",
  ladderize = T,
  lineend = "square", linewidth = 0.3, aes(color = scale(log(Mean_Rate), scale = F))
) +
  scale_color_viridis_c(
    name = "log(relative evolutionary rate)", option = "turbo", direction = 1,
    breaks = function(x) c(min(x), max(x)),
    labels = c("Low", "High")
  )
p_wholebrain_8 <- revts(p_wholebrain_8)
p_wholebrain_8 <- p_wholebrain_8 + ggnewscale::new_scale_color()
p_wholebrain_8 <- p_wholebrain_8 +
  # geom_cladelab(
  #   data = clade_ancestral_nodes,
  #   mapping = aes(
  #     node = MRCA,
  #     label = order,
  #     color = order
  #   ), textcolour = "black", offset = 3, offset.text = 3, barsize=2,extend=.5,angle=65,fontsize = 2.5, align=T)+
  geom_nodepoint(aes(subset = node %in% filter(tree_with_wholebrain_8@extraInfo, Pct.time.scaled >= 50 & Mean.Scalar >= 2)$node),
    col = "black",
    shape = 18,
    size = 1
  ) +
  scale_color_manual(values = grey_colors.1, guide = "none") +
  theme(
    legend.text = element_text(family = "Arial"),
    legend.title = element_text(size = 8, hjust = 0.5, family = "Arial"),
    plot.margin = unit(c(0, .3, 0.3, 0.3), "lines"),
    legend.position = "inside",
    legend.position.inside = c(0.3, 0.7),
    legend.direction = "horizontal",
    legend.key.height = unit(.7, "lines"),
    legend.title.position = "top",
    legend.background = element_blank()
  ) +
  geom_point(shape = 18, size = 2.7, color = "black", x = -152, y = 300) +
  annotate("text",
    x = -150, y = 300, label = "significant rate shift",
    hjust = 0, size = 3, family = "Arial"
  ) +
  xlim(-202, 20)
p_wholebrain_8
p_wholebrain_8$data$label[which(p_wholebrain_8$data$isTip == TRUE)] <- sub("_", " ", p_wholebrain_8$data$label[which(p_wholebrain_8$data$isTip == TRUE)])
p_wholebrain_8 <- p_wholebrain_8 + geom_tiplab(align = T, fontface = 3, size = .7, linesize = .2)

# Export:


showtext_auto()
showtext_opts(dpi = 300)

for (i in 1:8) {
  plot_obj <- get(paste0("p_wholebrain_", i))

  cairo_pdf(
    filename = paste0("./Images/Figures/Figure_wholebrain_hyp_", i, ".pdf"),
    width = 7, height = 9
  )
  print(plot_obj)
  dev.off()
  ggsave(
    filename = paste0("./Images/Figures/Figure_wholebrain_hyp_", i, ".tif"), plot_obj,
    width = 7, height = 9, units = "in"
  )
}

showtext_auto(FALSE)
