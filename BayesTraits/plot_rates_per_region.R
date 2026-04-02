library(BayesTraitR)
library(ape)
library(viridis)
library(ggtree)
library(deeptime)
library(phytools)
library(viridis)
library(tidyverse)
library(patchwork)
library(cowplot)

source("./scripts/utilities/BayesTraitsPlottingFunctions.R")

time_tree <- read.nexus("./phylogeny_construction/final_tree_hyp1.nex")

rjpp_results.size <- process_PPP("./BayesTraits/postprocessing/brainsize/linear_BM/brain_body_linear_BM_processed.txt",
  phy = ladderize(time_tree),
  save_summary_trees = F,
  col.palette = "viridis"
)
rjpp_results.wholebrain <- process_PPP("./BayesTraits/postprocessing/whole_brain/hyp_1/whole_brain_lambda_results.txt",
  phy = ladderize(time_tree),
  save_summary_trees = F,
  col.palette = "viridis"
)
rjpp_results.cerebellum <- process_PPP("./BayesTraits/postprocessing/cerebellum/hyp_1/cerebellum_lambda_results.txt",
  phy = ladderize(time_tree),
  save_summary_trees = F,
  col.palette = "viridis"
)
rjpp_results.neocortex <- process_PPP("./BayesTraits/postprocessing/neocortex/hyp_1/neocortex_lambda_results.txt",
  phy = ladderize(time_tree),
  save_summary_trees = F,
  col.palette = "viridis"
)
rjpp_results.olfactory <- process_PPP("./BayesTraits/postprocessing/olfactory/hyp_1/olfactory_lambda_results.txt",
  phy = ladderize(time_tree),
  save_summary_trees = F,
  col.palette = "viridis"
)
rjpp_results.brainstem <- process_PPP("./BayesTraits/postprocessing/brainstem/hyp_1/brainstem_lambda_results.txt",
  phy = ladderize(time_tree),
  save_summary_trees = F,
  col.palette = "viridis"
)


# get a relative rate at time for each dataset
RaT.brainsize <- rate.at.time.df(timeslices = 0.1, obj = rjpp_results.size$all.res[[1]], plot = F, relative.rates = "mean")
RaT.wholebrain <- rate.at.time.df(timeslices = 0.1, obj = rjpp_results.wholebrain$all.res[[1]], plot = F, relative.rates = "mean")
RaT.cerebellum <- rate.at.time.df(timeslices = 0.1, obj = rjpp_results.cerebellum$all.res[[1]], plot = F, relative.rates = "mean")
RaT.neocortex <- rate.at.time.df(timeslices = 0.1, obj = rjpp_results.neocortex$all.res[[1]], plot = F, relative.rates = "mean")
RaT.olfactory <- rate.at.time.df(timeslices = 0.1, obj = rjpp_results.olfactory$all.res[[1]], plot = F, relative.rates = "mean")
RaT.brainstem <- rate.at.time.df(timeslices = 0.1, obj = rjpp_results.brainstem$all.res[[1]], plot = F, relative.rates = "mean")

# extract the rates through time

rate.mean.brainsize <- extract.stat(RaT.brainsize, stat = "mean", plot = F, range = "confidence")
rate.mean.brainsize <- rate.mean.brainsize %>% mutate(trait = "Relative Brain Size")
rate.mean.wholebrain <- extract.stat(RaT.wholebrain, stat = "mean", plot = F, range = "confidence")
rate.mean.wholebrain <- rate.mean.wholebrain %>% mutate(trait = "Whole Endocast Shape")
rate.mean.cerebellum <- extract.stat(RaT.cerebellum, stat = "mean", plot = F, range = "confidence")
rate.mean.cerebellum <- rate.mean.cerebellum %>% mutate(trait = "Cerebellum Shape")
rate.mean.neocortex <- extract.stat(RaT.neocortex, stat = "mean", plot = F, range = "confidence")
rate.mean.neocortex <- rate.mean.neocortex %>% mutate(trait = "Neocortex Shape")
rate.mean.olfactory <- extract.stat(RaT.olfactory, stat = "mean", plot = F, range = "confidence")
rate.mean.olfactory <- rate.mean.olfactory %>% mutate(trait = "Olfactory Bulb Shape")
rate.mean.brainstem <- extract.stat(RaT.brainstem, stat = "mean", plot = F, range = "confidence")
rate.mean.brainstem <- rate.mean.brainstem %>% mutate(trait = "Brainstem Shape")
rates_through_time <- bind_rows(
  rate.mean.brainsize,
  rate.mean.wholebrain,
  rate.mean.cerebellum,
  rate.mean.neocortex,
  rate.mean.olfactory,
  rate.mean.brainstem
)

annotation <- data.frame(
  x = c(52, 72),
  y = c(3, 3),
  label = c("K-PG/nBoundary", "PETM")
)

palette.colors(n = 10, palette = "Tableau 10")
# "#4E79A7"  "#E15759" "#76B7B2" "#59A14F" "#EDC948"
region_colors <- c("black", "#E15759", "#B07AA1", "#F28E2B", "#59A14F", "#EDC948")
trait_names <- c(
  "Relative Brain Size",
  "Whole Endocast Shape",
  "Olfactory Bulb Shape",
  "Neocortex Shape",
  "Cerebellum Shape",
  "Brainstem Shape"
)


names(region_colors) <- trait_names


single_RTT_plot <- ggplot(rates_through_time, aes(x = -1 * (time), color = trait, linetype = trait, linewidth = trait, fill = trait)) +
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
    breaks = c(
      "Relative Brain Size",
      "Whole Endocast Shape",
      "Olfactory Bulb Shape",
      "Neocortex Shape",
      "Cerebellum Shape",
      "Brainstem Shape"
    ),
    values = region_colors
  ) +
  scale_linetype_manual(
    name = "",
    breaks = c(
      "Relative Brain Size",
      "Whole Endocast Shape",
      "Olfactory Bulb Shape",
      "Neocortex Shape",
      "Cerebellum Shape",
      "Brainstem Shape"
    ),
    values = c(1, 1, 1, 1, 1, 1)
  ) +
  scale_linewidth_manual(
    name = "",
    breaks = c(
      "Relative Brain Size",
      "Whole Endocast Shape",
      "Olfactory Bulb Shape",
      "Neocortex Shape",
      "Cerebellum Shape",
      "Brainstem Shape"
    ),
    values = c(1.2, 1.2, 0.8, 0.8, 0.8, 0.8)
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
  annotate(
    geom = "text", label = "A",
    x = 97, y = 3.8,
    fontface = "bold", size = 4
  ) +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.8, 0.75),
    legend.key.width = unit(1, "cm")
  )


highlight_data <- rates_through_time %>%
  mutate(facet_trait = trait)

tag_data <- data.frame(
  facet_trait = factor(trait_names, levels = trait_names),
  label = c("B", "C", "D", "E", "F", "G")
)

geom_text(
  data = tag_data,
  aes(x = 95, y = 3.8, label = label),
  inherit.aes = FALSE,
  fontface = "bold", size = 6
)

full_RTT_data <- rates_through_time %>%
  crossing(facet_trait = unique(rates_through_time$trait)) %>%
  mutate(
    is_focal = trait == facet_trait,
    line_color = ifelse(is_focal, region_colors[trait], "grey70")
  ) %>%
  mutate(facet_trait = factor(facet_trait, levels = c(
    "Relative Brain Size",
    "Whole Endocast Shape",
    "Olfactory Bulb Shape",
    "Neocortex Shape",
    "Cerebellum Shape",
    "Brainstem Shape"
  )))

facet_plot <- ggplot(full_RTT_data, aes(x = -1 * (time), y = rate, group = trait)) +
  scale_x_continuous(
    trans = scales::reverse_trans(),
    breaks = seq(100, 0, by = -10),
    limits = c(100, 0),
    labels = c("", "90", "80", "70", "60", "50", "40", "30", "20", "10", "0")
  ) +
  geom_line(
    data = filter(full_RTT_data, !is_focal),
    aes(linetype = trait, linewidth = trait),
    color = "grey70"
  ) +
  geom_line(
    data = filter(full_RTT_data, is_focal),
    aes(color = line_color, linetype = trait, linewidth = trait)
  ) +
  scale_color_identity() +
  facet_wrap(~facet_trait) +
  theme_minimal() +
  ylim(0, 4) +
  scale_linetype_manual(
    name = "Trait",
    breaks = c(
      "Relative Brain Size",
      "Whole Endocast Shape",
      "Olfactory Bulb Shape",
      "Neocortex Shape",
      "Cerebellum Shape",
      "Brainstem Shape"
    ),
    values = c(1, 1, 1, 1, 1, 1)
  ) +
  scale_linewidth_manual(
    name = "Trait",
    breaks = c(
      "Relative Brain Size",
      "Whole Endocast Shape",
      "Olfactory Bulb Shape",
      "Neocortex Shape",
      "Cerebellum Shape",
      "Brainstem Shape"
    ),
    values = c(0.6, 0.6, 0.6, 0.6, 0.6, 0.6)
  ) +
  geom_vline(xintercept = 66) +
  geom_vline(xintercept = 55.8, linetype = "dashed") +
  coord_geo(center_end_labels = T, height = unit(1, "lines"), size = 4) +
  geom_text(
    data = tag_data,
    aes(x = 95, y = 3.6, label = label),
    inherit.aes = FALSE,
    fontface = "bold", size = 4
  ) +
  labs(x = "Ma before present", y = "Relative Rate of Phenotypic Evolution") +
  theme(
    legend.position = "none",
    legend.key.width = unit(2, "cm")
  )

# join together the full plot and faceted plots
composite_plot <- single_RTT_plot / facet_plot
# save
ggsave(
  filename = "./Images/Figures/Figure_3.pdf",
  plot = composite_plot,
  width = 7,
  height = 7,
  units = "in"
)
ggsave(
  filename = "./Images/Figures/Figure_3.tif",
  plot = composite_plot,
  width = 7,
  height = 7,
  units = "in",
  dpi = 600
)


####### plot the trees


data <- read.csv("./input_data/mammal_data.csv", row.names = 1)
cladenames <- levels(as.factor(data$order))


clade_ancestral_nodes <- data.frame(order = cladenames, MRCA = NA)
for (i in 1:nrow(clade_ancestral_nodes)) {
  these_tips <- rownames(data)[which(data$order == clade_ancestral_nodes$order[i])]
  if (length(these_tips) > 1) {
    this_anc <- getMRCA(time_tree, these_tips)
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

rate_data_olfactory <- rjpp_results.olfactory$all.res[[1]]
rate_data_olfactory <- rate_data_olfactory %>%
  dplyr::rename(., node = child.node) %>%
  dplyr::rename(., Mean_Rate = Mean.SigV)
tree_with_olfactory <- treeio::full_join(time_tree, y = rate_data_olfactory, by = "node")

rate_data_neocortex <- rjpp_results.neocortex$all.res[[1]]
rate_data_neocortex <- rate_data_neocortex %>%
  dplyr::rename(., node = child.node) %>%
  dplyr::rename(., Mean_Rate = Mean.SigV)
tree_with_neocortex <- treeio::full_join(time_tree, y = rate_data_neocortex, by = "node")

rate_data_cerebellum <- rjpp_results.cerebellum$all.res[[1]]
rate_data_cerebellum <- rate_data_cerebellum %>%
  dplyr::rename(., node = child.node) %>%
  dplyr::rename(., Mean_Rate = Mean.SigV)
tree_with_cerebellum <- treeio::full_join(time_tree, y = rate_data_cerebellum, by = "node")

rate_data_brainstem <- rjpp_results.brainstem$all.res[[1]]
rate_data_brainstem <- rate_data_brainstem %>%
  dplyr::rename(., node = child.node) %>%
  dplyr::rename(., Mean_Rate = Mean.SigV)
tree_with_brainstem <- treeio::full_join(time_tree, y = rate_data_brainstem, by = "node")

p_olfactory <- ggtree(tree_with_olfactory, layout = "rectangular", ladderize = T, lineend = "square", linewidth = 0.3, aes(color = scale(log(Mean_Rate), scale = F))) +
  coord_flip() +
  scale_color_viridis_c(
    name = "log(relative evolutionary rate)", option = "turbo", direction = 1,
    breaks = function(x) c(min(x), max(x)),
    labels = c("Low", "High")
  )
p_olfactory <- revts(p_olfactory)
p_olfactory <- p_olfactory + ggnewscale::new_scale_color()
p_olfactory <- p_olfactory +
  # geom_cladelab(
  #   data = clade_ancestral_nodes,
  #   mapping = aes(
  #     node = MRCA,
  #     label = order,
  #     color = order
  #   ), textcolour = "black", offset = 3, offset.text = 3, barsize=2,extend=.5,angle=65,fontsize = 2.5, align=T)+
  geom_nodepoint(aes(subset = node %in% filter(tree_with_olfactory@extraInfo, Pct.time.scaled >= 50 & Mean.Scalar >= 2)$node),
    col = "black",
    shape = 18,
    size = 2.7
  ) +
  scale_color_manual(values = grey_colors.1, guide = "none") +
  theme(
    legend.text = element_text(family = "Arial"),
    legend.title = element_text(size = 8, hjust = 0.5, family = "Arial"),
    plot.margin = unit(c(0, .3, 0.3, 0.3), "lines"),
    legend.position = "inside",
    legend.position.inside = c(0.7, 0.3),
    legend.direction = "horizontal",
    legend.key.height = unit(.7, "lines"),
    legend.title.position = "top",
    legend.background = element_blank()
  ) +
  geom_point(shape = 18, size = 2.7, color = "black", x = -140, y = 165) +
  annotate("text",
    x = -140, y = 170, label = "significant rate shift",
    hjust = 0, size = 3, family = "Arial"
  ) +
  xlim(-202, 20)
p_olfactory


p_neocortex <- ggtree(tree_with_neocortex, layout = "rectangular", ladderize = T, lineend = "square", linewidth = 0.3, aes(color = scale(log(Mean_Rate), scale = F))) +
  coord_flip() +
  scale_color_viridis_c(name = "log(rate of evolution\nof neocortex shape)", option = "turbo", direction = 1)
p_neocortex <- revts(p_neocortex)
p_neocortex <- p_neocortex + ggnewscale::new_scale_color()
p_neocortex <- p_neocortex +
  # geom_cladelab(
  #   data = clade_ancestral_nodes,
  #   mapping = aes(
  #     node = MRCA,
  #     label = order,
  #     color = order
  #   ), textcolour = "black", offset = 3, offset.text = 3, barsize=2,extend=.5,angle=65,fontsize = 2.5, align=T)+
  geom_nodepoint(aes(subset = node %in% filter(tree_with_neocortex@extraInfo, Pct.time.scaled >= 50 & Mean.Scalar >= 2)$node),
    col = "black",
    shape = 18,
    size = 2.7
  ) +
  scale_color_manual(values = grey_colors.1, guide = "none") +
  theme(
    legend.position = "none",
    plot.margin = unit(c(0, .3, 0.3, 0.3), "lines")
  ) +
  xlim(-202, 20)
p_neocortex

p_cerebellum <- ggtree(tree_with_cerebellum, layout = "rectangular", ladderize = T, lineend = "square", linewidth = 0.3, aes(color = scale(log(Mean_Rate), scale = F))) +
  coord_flip() +
  scale_color_viridis_c(name = "log(rate of evolution\nof cerebellum shape)", option = "turbo", direction = 1)
p_cerebellum <- revts(p_cerebellum)
p_cerebellum <- p_cerebellum + ggnewscale::new_scale_color()
p_cerebellum <- p_cerebellum +
  # geom_cladelab(
  #   data = clade_ancestral_nodes,
  #   mapping = aes(
  #     node = MRCA,
  #     label = order,
  #     color = order
  #   ), textcolour = "black", offset = 3, offset.text = 3, barsize=2,extend=.5,angle=65,fontsize = 2.5, align=T)+
  geom_nodepoint(aes(subset = node %in% filter(tree_with_cerebellum@extraInfo, Pct.time.scaled >= 50 & Mean.Scalar >= 2)$node),
    col = "black",
    shape = 18,
    size = 2.7
  ) +
  scale_color_manual(values = grey_colors.1, guide = "none") +
  theme(
    legend.position = "none",
    plot.margin = unit(c(0, .3, 0.3, 0.3), "lines"),
  ) +
  xlim(-202, 20)
p_cerebellum

p_brainstem <- ggtree(tree_with_brainstem, layout = "rectangular", ladderize = T, lineend = "square", linewidth = 0.3, aes(color = scale(log(Mean_Rate), scale = F))) +
  coord_flip() +
  scale_color_viridis_c(name = "log(rate of evolution\nof brainstem shape)", option = "turbo", direction = 1)
p_brainstem <- revts(p_brainstem)
p_brainstem <- p_brainstem + ggnewscale::new_scale_color()
p_brainstem <- p_brainstem +
  # geom_cladelab(
  #   data = clade_ancestral_nodes,
  #   mapping = aes(
  #     node = MRCA,
  #     label = order,
  #     color = order
  #   ), textcolour = "black", offset = 3, offset.text = 3, barsize=2,extend=.5,angle=65,fontsize = 2.5, align=T)+
  geom_nodepoint(aes(subset = node %in% filter(tree_with_brainstem@extraInfo, Pct.time.scaled >= 50 & Mean.Scalar >= 2)$node),
    col = "black",
    shape = 18,
    size = 2.7
  ) +
  scale_color_manual(values = grey_colors.1, guide = "none") +
  theme(
    legend.position = "none",
    plot.margin = unit(c(0, .3, 0.3, 0.3), "lines")
  ) +
  xlim(-202, 20)
p_brainstem


p_strip <- revts(ggtree(tree_with_olfactory, layout = "rectangular", color = "transparent", ladderize = T) + coord_flip()) +
  geom_cladelab(
    data = clade_ancestral_nodes,
    mapping = aes(
      node = MRCA,
      label = order,
      color = order
    ),
    textcolour = "black", offset = 0,
    offset.text = 2, barsize = 2,
    extend = .5, angle = 65, fontsize = 1.5, align = T
  ) +
  scale_color_manual(values = grey_colors.1, guide = "none") +
  theme_void() +
  theme(legend.position = "none", plot.margin = unit(c(0, .3, -.01, 0.3), "lines")) +
  xlim(0, 20)
p_strip

###########
###########


plotgrid1 <- plot_grid(
  p_strip,
  p_olfactory,
  p_neocortex,
  p_cerebellum,
  p_brainstem,
  nrow = 5,
  rel_heights = c(3.5, 10, 10, 10, 10),
  align = "v",
  axis = "lr"
)

plotgrid1
cairo_pdf(filename = "./Images/Figures/Figure_S14.pdf", width = 7, height = 10)
plotgrid1
dev.off()
