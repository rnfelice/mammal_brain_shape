
pak::pak("jakeberv/bifrost")

library(bifrost)
library(ape)
library(phytools)
library(mvMORPH)


set.seed(1)

# load phylogeny, trait data and landmark data
tree_hyp1 <- read.nexus('./phylogeny_construction/final_tree_hyp1.nex') # tree hypothesis 1
classifier <- read.csv('./input_data/mammal_data.csv')

# procrustes-aligned landmark data
n_landmarks <- 122
n_dimensions <- 3
raw_data <- read.csv("./input_data/procrutes_aligned_coords.csv", row.names = 1)

data<-raw_data[tree_hyp1$tip.label,]
geiger::name.check(tree_hyp1, data)


base <- phytools::paintBranches(
  tree_hyp1,
  edge      = unique(tree_hyp1$edge[, 2]),
  state     = "0",
  anc.state = "0"
)

identical(rownames(data), base$tip.label)

# Run bifrost's greed# Run bifrost's greed# Run bifrost's greedy search for shifts
res <- searchOptimalConfiguration(
  baseline_tree              = base,
  trait_data                 = data,
  formula                    = "trait_data ~ 1",
  min_descendant_tips        = 10,
  num_cores                  = 4,
  shift_acceptance_threshold = 20,  # conservative GIC threshold
  IC                         = "GIC",
  plot                       = TRUE,
  store_model_fit_history    = FALSE,
  method = "EmpBayes"
)


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
  ladderize.simmap(res$tree_no_uncertainty_untransformed),
  direction='rightwards',
  fsize = .2,
  colors = generateViridisColorScale2(res$model_no_uncertainty$param, option="turbo")$NamedColors
)
