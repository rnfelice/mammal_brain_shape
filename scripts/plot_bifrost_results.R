#plot bifrost results
library(bifrost)
library(ape)
library(phytools)
library(mvMORPH)
library(geomorph)
library(viridis)
library(deeptime)
source("./scripts/utilities/BayesTraitsPlottingFunctions.R")
load("~/Downloads/bifrost_results_temp/bifrost_res_july-allometry-mass.rda")
load("~/Downloads/bifrost_results_temp/bifrost_res_july-allometry.rda")
load("~/Downloads/bifrost_results_temp/bifrost_res_july-highmin.rda")
load("~/Downloads/bifrost_results_temp/bifrost_res_july-lowmin.rda")
load("~/Downloads/bifrost_results_temp/bifrost_res_july.rda")




generateViridisColorScale2 <- function (params, option = "D") 
{
  sorted_indices <- order(params)
  sorted_params <- params[sorted_indices]
  normalized_sorted_params <- (sorted_params - min(sorted_params))/(max(sorted_params) - 
                                                                      min(sorted_params))
  colors <- viridis(length(normalized_sorted_params), option = option)
  named_sorted_colors <- setNames(colors, names(sorted_params))
  param_color_mapping <- setNames(sorted_params, names(sorted_params))
  return(list(NamedColors = named_sorted_colors, ParamColorMapping = param_color_mapping))
}
plotSimmap(
  bifrost_tree,
  direction='rightwards',
  fsize = .2,
  colors = generateViridisColorScale2(log(res$model_no_uncertainty$param), option="turbo")$NamedColors
)


branch_rates_bifrost1 <- rateMap(
  list(res),
  summary = "branch",
  log=F,
  weights="ic"
)

branch_rates_bifrost_lowmin <- rateMap(
  list(res.lowmin),
  summary = "branch",
  log=F,
  weights="ic"
)

branch_rates_bifrost_allom <- rateMap(
  list(res.allom),
  summary = "branch",
  log=F,
  weights="ic"
)

branch_rates_bifrost_highmin <- rateMap(
  list(res.highmin),
  summary = "branch",
  log=F,
  weights="ic"
)

branch_rates_bifrost_allom_mass <- rateMap(
  list(res.allom.logmass),
  summary = "branch",
  log=F,
  weights="ic"
)

results_bifrost_df<-rate.at.time.df.bifrost(timeslices = .1, obj=as_tibble(branch_rates_bifrost1$intervals),plot=T,relative.rates="mean")
results_bifrost_lowmin_df<-rate.at.time.df.bifrost(timeslices = .1, obj=as_tibble(branch_rates_bifrost_lowmin$intervals),plot=T,relative.rates="mean")
results_bifrost_highmin_df<-rate.at.time.df.bifrost(timeslices = .1, obj=as_tibble(branch_rates_bifrost_highmin$intervals),plot=T,relative.rates="mean")
results_bifrost_allom_df<-rate.at.time.df.bifrost(timeslices = .1, obj=as_tibble(branch_rates_bifrost_allom$intervals),plot=T,relative.rates="mean")
results_bifrost_allom_mass_df<-rate.at.time.df.bifrost(timeslices = .1, obj=as_tibble(branch_rates_bifrost_allom_mass$intervals),plot=T,relative.rates="mean")

rate.mean.bifrost <- extract.stat(results_bifrost_df, stat = "mean", plot = F, range = "confidence")
rate.mean.bifrost <- rate.mean.bifrost %>% mutate(trait = "Bifrost- moderate min cladesize")

rate.mean.bifrost_highmin <- extract.stat(results_bifrost_highmin_df, stat = "mean", plot = F, range = "confidence")
rate.mean.bifrost_highmin <- rate.mean.bifrost_highmin %>% mutate(trait = "Bifrost- high min cladesize")

rate.mean.bifrost_lowmin <- extract.stat(results_bifrost_lowmin_df, stat = "mean", plot = F, range = "confidence")
rate.mean.bifrost_lowmin <- rate.mean.bifrost_lowmin %>% mutate(trait = "Bifrost- low min cladesize")

rate.mean.bifrost_allom <- extract.stat(results_bifrost_allom_df, stat = "mean", plot = F, range = "confidence")
rate.mean.bifrost_allom <- rate.mean.bifrost_allom %>% mutate(trait = "Bifrost- with centroid size")

rate.mean.bifrost_allom_mass <- extract.stat(results_bifrost_allom_mass_df, stat = "mean", plot = F, range = "confidence")
rate.mean.bifrost_allom_mass <- rate.mean.bifrost_allom_mass %>% mutate(trait = "Bifrost- with mass")


######################
##plotthe tree with branch rates
#####################

# Build the branch-rate map used beneath the shift-node annotations.
brain_rate_map <- rateMap(res, progress = FALSE)
brain_transitions <- shift_transitions(res, include_root = FALSE)
# Recreate the manuscript's Fisher-binned log-rate scale and color palette.
brain_rate_breaks <- classInt::classIntervals(
  log(res$model_no_uncertainty$param),
  n = 10,
  style = "fisher"
)$brks

brain_rate_breaks[c(1, length(brain_rate_breaks))] <-
  range(brain_rate_map$intervals$value)

brain_rate_palette <- grDevices::colorRampPalette(
  rev(RColorBrewer::brewer.pal(11, "RdYlBu"))
)(length(brain_rate_breaks) - 1L)

# Apply those fixed categories to a reusable tree view.
brain_shift_tree <- rateMapView(
  brain_rate_map,
  color_mode = "category",
  category_breaks = brain_rate_breaks,
  palette = brain_rate_palette,
  legend_title = "Log fitted BMM rate"
)
brain_shift_marks <- shift_node_marks(
  brain_transitions,
  ic_weights = res$ic_weights,
  marker_base_cex = 1.5,
  marker_scale_factor = 0.08,
  marker_transform = "square_root",
  letter_order = "child_state"
)

brain_shift_marks$summary
plot(
  brain_shift_tree,
  type = "arc",
  fsize=.2,
  show_tip_labels = T,
  lwd = c(1.20, 5.0),
  legend = 42,
  legend_fsize = 0.72,
  legend_digits = 3,
  arc_height = 0.75,
  mar = c(0.2, 0.2, 0.2, 0.2)
)
brain_shift_marks$increase_key
# Add only supported rate-increase markers to the active tree plot.
plot(
  brain_shift_marks,
  rate_changes = "increase",
  show_low_support = FALSE,
  show_legend = FALSE
)



#load some bayestraits results for comparison

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

# get a relative rate at time for each dataset
RaT.brainsize <- rate.at.time.df(timeslices = 0.1, obj = rjpp_results.size$all.res[[1]], plot = F, relative.rates = "mean")
RaT.wholebrain <- rate.at.time.df(timeslices = 0.1, obj = rjpp_results.wholebrain$all.res[[1]], plot = F, relative.rates = "mean")

rate.mean.brainsize <- extract.stat(RaT.brainsize, stat = "mean", plot = F, range = "confidence")
rate.mean.brainsize <- rate.mean.brainsize %>% mutate(trait = "Relative Brain Size")
rate.mean.wholebrain <- extract.stat(RaT.wholebrain, stat = "mean", plot = F, range = "confidence")
rate.mean.wholebrain <- rate.mean.wholebrain %>% mutate(trait = "Whole Endocast Shape")


rates_through_time <- bind_rows(
  rate.mean.brainsize,
  rate.mean.wholebrain,
  rate.mean.bifrost,
  rate.mean.bifrost_lowmin,
  rate.mean.bifrost_highmin,
  rate.mean.bifrost_allom,
  rate.mean.bifrost_allom_mass,
)

annotation <- data.frame(
  x = c(52, 72),
  y = c(3, 3),
  label = c("K-PG/nBoundary", "PETM")
)

palette.colors(n = 10, palette = "Tableau 10")
# "#4E79A7"  "#E15759" "#76B7B2" "#59A14F" "#EDC948"
region_colors <- c("black", "#E15759", "#B07AA1","#4E79A7",  "#EDC948", "#76B7B2", "#59A14F")
trait_names <- c(
  "Relative Brain Size",
  "Whole Endocast Shape",
  "Bifrost- moderate min cladesize",
  "Bifrost- low min cladesize",
 "Bifrost- high min cladesize",
 "Bifrost- with centroid size",
 "Bifrost- with mass"
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
  "Bifrost- moderate min cladesize",
  "Bifrost- low min cladesize",
 "Bifrost- high min cladesize",
 "Bifrost- with centroid size",
 "Bifrost- with mass"
)
,
    values = region_colors
  ) +
  scale_linetype_manual(
    name = "",
    breaks = c(
  "Relative Brain Size",
  "Whole Endocast Shape",
  "Bifrost- moderate min cladesize",
  "Bifrost- low min cladesize",
 "Bifrost- high min cladesize",
 "Bifrost- with centroid size",
 "Bifrost- with mass"
)
,
    values = c(1, 1, 1, 1, 1, 1, 1)
  ) +
  scale_linewidth_manual(
    name = "",
    breaks = c(
  "Relative Brain Size",
  "Whole Endocast Shape",
  "Bifrost- moderate min cladesize",
  "Bifrost- low min cladesize",
 "Bifrost- high min cladesize",
 "Bifrost- with centroid size",
 "Bifrost- with mass"
)
,
    values = c(1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2)
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

single_RTT_plot

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
