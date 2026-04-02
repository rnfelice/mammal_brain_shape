library(BayesTraitR)
library(mcmc)
library(coda)
library(tidyverse)
# get log files
bt_output_folder <- "/Volumes/FeliceLabHD/bt_mammal_brains_26_3/brainsize"
all_files <- list.files(bt_output_folder, recursive = T)

log_files <- all_files[str_detect(all_files, "\\.txt\\.Log\\.txt$")]

log_file_table <- data.frame(filepath = log_files) %>%
  mutate(
    evo_model = str_extract(filepath, "^[^/]+"),
    run = str_extract(filepath, "(?<=-)(\\d+)(?=\\.txt)"),
    filename = basename(filepath)
  )

cols_wanted <- c("Lh", "Alpha", "Alpha 1", "Sigma^2 1", "Beta 1", "Beta 2", "Var", "No RJ Local Branch", "No RJ Local Node", "Lambda", "Kappa", "Delta", "OU")

convergence_results_size <- log_file_table %>%
  group_by(evo_model) %>%
  summarise(
    mpsrf = {
      mcmc_list <- map(filepath, ~ {
        readBTlog(paste0(bt_output_folder, "/", .x)) %>%
          select(any_of(cols_wanted)) %>%
          as.mcmc()
      })
      gelman.diag(mcmc.list(mcmc_list))$mpsrf
    }
  )

convergence_results

# shape analyses
bt_output_folder <- "/Volumes/FeliceLabHD/BT_mammal_brains_26/"

log_files <- list.files(bt_output_folder, recursive = TRUE, full.names = FALSE, pattern = "*Log.txt$")

# remove relative brain size analyses
log_files <- log_files[-which(str_detect(log_files, "brainsize"))]

log_file_table <- data.frame(filepath = log_files) %>%
  mutate(
    brain_region = str_extract(filepath, "^[^/]+"),
    rate_type    = str_extract(filepath, "(?<=SCORES_)(single_rate|var_rate)"),
    evo_model    = str_extract(filepath, "(?<=tree)(BM|OU|delta|kappa|lambda)"),
    run          = str_extract(filepath, "(?<=_)([ab])(?=_run)"),
    filename     = basename(filepath)
  )

convergence_results_shape <- log_file_table %>%
  group_by(brain_region, rate_type, evo_model) %>%
  summarise(
    mpsrf = {
      mcmc_list <- map(filepath, ~ {
        readBTlog(paste0(bt_output_folder, "/", .x)) %>%
          select(any_of(cols_wanted)) %>%
          as.mcmc()
      })
      gelman.diag(mcmc.list(mcmc_list))$mpsrf
    }
  )


log1 <- readBTlog(paste0(bt_output_folder, "/", "cerebellum/pPCscores/Phylo_PC_SCORES_var_rate_treeOU_a_run.txt.Log.txt")) |>
  select(any_of(cols_wanted)) %>%
  as.mcmc()
log2 <- readBTlog(paste0(bt_output_folder, "/", "cerebellum/pPCscores/Phylo_PC_SCORES_var_rate_treeOU_b_run.txt.Log.txt")) |>
  select(any_of(cols_wanted)) %>%
  as.mcmc()
list1 <- mcmc.list(as.mcmc(log1[1:98000, ]), as.mcmc(log2[1:98000, ]))
gelman.diag(list1)
