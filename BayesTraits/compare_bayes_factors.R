library(tidyverse)

source("./scripts/utilities/BayesTraitsModelComparison.R")


# list the stones files- these include the marginal likelihood of each model
bt_output_folder <- "/Volumes/FeliceLabHD/BT_mammal_brains_26_3/"

stones_files <- list.files(
  bt_output_folder,
  recursive = TRUE,
  full.names = FALSE,
  pattern = "*Stones.txt"
)

# sort the files by dataset, evolutionary model, and run
stones_file_table <- data.frame(filepath = stones_files) %>%
  mutate(
    brain_region = str_extract(filepath, "^[^/]+"),
    tree         = str_extract(filepath, "(?<=hyp_)[1-8]"),
    rate_type    = str_extract(filepath, "(?<=SCORES_)(single_rate|var_rate)"),
    evo_model    = str_extract(filepath, "(?<=tree)(BM|OU|delta|kappa|lambda)"),
    run          = str_extract(filepath, "(?<=_)([ab])(?=_run)"),
    filename     = basename(filepath)
  )

# separate the relative brain size analyses from the brain shape analyses
stones_file_table_size <- stones_file_table %>% filter(brain_region == "brain_size")

stones_file_table <- stones_file_table %>% filter(brain_region != "brain_size")

##########
# get model likelihoods for brain shape analyses
##########
# clean up the names of the models

stones_file_table <- stones_file_table %>% filter(run == "a")
stones_file_table <- stones_file_table %>%
  mutate(rate_type = case_match(
    rate_type,
    "var_rate" ~ "Variable Rates",
    "single_rate" ~ "Single Rate"
  )) %>%
  mutate(evo_model = case_match(
    evo_model,
    "lambda" ~ "Lambda",
    "delta" ~ "Delta",
    "kappa" ~ "Kappa",
    .default = evo_model
  ))

# extract the marginal likelihoods from the stones files

stones_results <- stones_file_table %>%
  group_by(brain_region, tree) %>%
  group_map(~ {
    stone_paths <- paste0(bt_output_folder, .x$filepath)
    labels <- paste(.x$rate_type, .x$evo_model, sep = "\n")
    list(
      brain_region = .y$brain_region,
      tree = .y$tree,
      stones = getStones(stone_paths, labels = labels)
    )
  })


# set the output folder:
figure_folder <- "./Images/Figures/BF_plots/"

# make plots and output them
for (i in 1:length(stones_results)) {
  g <- plot(stones_results[[i]]$stones)
  g <- g + ggtitle("Bayes Factor Comparison",
    subtitle = paste0(
      "Dataset: ",
      stringr::str_to_title(gsub("_", " ", stones_results[[i]]$brain_region)),
      "; Phylogenetic Hypothesis: Tree ",
      stones_results[[i]]$tree
    )
  )
  g <- g +
    theme_minimal(base_size = 8) +
    theme(
      legend.position = "bottom",
      axis.text.x.top = element_text(
        angle = -45,
        hjust = 1,
        vjust = 0,
        size = 9
      ),
      axis.text.y = element_text(size = 9),
      plot.title = element_text(size = 14),
      plot.subtitle = element_text(size = 12),
      legend.text = element_text(size = 9),
      legend.title = element_text(size = 9)
    )
  # save:
  ggsave(paste0(figure_folder, "Figure_S", 20 + i, ".pdf"),
    plot = g,
    height = 7,
    width = 7,
    units = "in"
  )
  ggsave(paste0(figure_folder, "Figure_S", 20 + i, ".tif"),
    plot = g,
    height = 7,
    width = 7,
    units = "in"
  )
}


##########
# get model likelihoods for relative brain size analyses
##########
# clean up the names of the models
stones_file_table_size <- stones_file_table_size %>%
  mutate(
    tree = 1,
    rate_type = "Variable Rates",
    evo_model = str_extract(filepath, "(?<=_)(BM|OU|delta|kappa|lambda)(?=/Brain)"),
    run = str_extract(filepath, "(?<=BT-)(\\d+)(?=\\.txt\\.Stones)"),
    regression_type = str_extract(filepath, "(?<=/)[^/_]+(?=_)")
  )
stones_file_table_size <- stones_file_table_size %>%
  mutate(evo_model = case_match(
    evo_model,
    "lambda" ~ "Lambda",
    "delta" ~ "Delta",
    "kappa" ~ "Kappa",
    .default = evo_model
  )) %>%
  mutate(regression_type = case_match(
    regression_type,
    "clade" ~ "Clade-specific slopes",
    "linear" ~ "Linear Regression",
    "quadratic" ~ "Quadratic Regression"
  ))

# extract the marginal likelihoods from the stones files

stones_file_table_size <- stones_file_table_size %>% filter(run == "001")
stones_results_size <- stones_file_table_size %>%
  group_by(brain_region, tree) %>%
  group_map(~ {
    stone_paths <- paste0(bt_output_folder, .x$filepath)
    labels <- paste(.x$regression_type, .x$evo_model, sep = "\n")
    list(
      brain_region = .y$brain_region,
      tree = .y$tree,
      stones = getStones(stone_paths, labels = labels)
    )
  })

# plot and save:

g <- plot(stones_results_size[[1]]$stones)
g <- g + ggtitle("Bayes Factor Comparison",
  subtitle = paste0("Dataset: Relative Brain Size Evolution; Phylogenetic Hypothesis: Tree 1")
)
g <- g +
  theme_minimal(base_size = 8) +
  theme(
    legend.position = "bottom",
    axis.text.x.top = element_text(
      angle = -45,
      hjust = 1,
      vjust = 0,
      size = 9
    ),
    axis.text.y = element_text(size = 9),
    plot.title = element_text(size = 14),
    plot.subtitle = element_text(size = 12),
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 9),
    plot.margin = margin(0.3, 10, 0.3, 0.3, unit = "lines")
  )

g
ggsave(paste0(figure_folder, "Figure_S", 20 + length(stones_results) + 1, ".pdf"),
  plot = g,
  height = 7,
  width = 7,
  units = "in"
)
ggsave(paste0(figure_folder, "Figure_S", 20 + length(stones_results) + 1, ".tif"),
  plot = g,
  height = 7,
  width = 7,
  units = "in"
)
