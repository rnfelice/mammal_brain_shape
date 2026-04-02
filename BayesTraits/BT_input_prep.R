# Prepping files for BayesTraits Rates Analysis ----------------------------------------------
library(ape)
library(phytools)
library(paleotree)
library(tidyverse)
library(geomorph)

# load data
# load landmarks
n_landmarks <- 122
n_dimensions <- 3
raw_data <- read.csv("./input_data/procrutes_aligned_coords.csv", row.names = 1)
Y.gpa.endo <- arrayspecs(raw_data, p = n_landmarks, k = n_dimensions)
# number if different tree hypotheses
ntrees <- 8

# define modules
module_list <- list(
  whole_brain = c(1:122),
  neocortex = c(2:3, 10:21, 52:59, 74:113),
  olfactory = c(1, 4, 22:29, 114:122),
  cerebellum = c(5, 9, 30:32, 45:51, 60:68),
  brainstem = c(6:8, 33:44, 69:73)
)

# definte function to write bayestraits comman files:
# the function then sets the prior for alpha to be a uniform distribution from the min to the max of that PC score (trait)
write_cmd_file <- function(x, varrate = TRUE, run, model = "BM", folder) {
  filename_cmd <- paste0(
    folder,
    "/scripts/command_files/BT_Control_",
    model,
    "_",
    run,
    "_control.cmd"
  )
  if (varrate == TRUE) {
    filename_cmd <- gsub("Control_", "Control_varrates_", filename_cmd)
  }
  sink(filename_cmd)
  cat(
    "7",
    "2",
    "stones 500 5000",
    "iterations 1000000000",
    "burnin 20000000",
    "sample 10000",
    "PriorAll gamma 5 12",
    sep = "\n"
  )
  cat("\n")
  for (k in 1:ncol(x)) {
    cat(paste0(
      "Prior Alpha-",
      k,
      " uniform ",
      signif(min(x[, k]), 5),
      " ",
      signif(max(x[, k]), 5)
    ))
    cat("\n")
  }
  if (model != "BM") {
    cat(model)
  }
  cat("\n")
  if (varrate == TRUE) {
    cat("varrates")
    cat("\n")
    cat("RJLockModel")
    cat("\n")
  }
  cat("run")
  sink()
}


list_of_models <- c("BM", "OU", "lambda", "kappa", "delta")

runs <- c("a", "b") # do two runs of each so that you can check convergence at the end


for (treenum in 1:ntrees) {
  this_tree <- read.nexus(
    file = paste0("./phylogeny_construction/final_tree_hyp", treenum, ".nex")
  )
  for (modulenum in 1:length(module_list)) {
    # define output folder for this tree and brain region
    outputfolder <- paste0(
      "./BayesTraits/",
      names(module_list)[modulenum],
      "/hyp_",
      treenum
    )
    # do phylogenetic PCA
    phypc <- gm.prcomp(
      A = Y.gpa.endo[module_list[[modulenum]], , ],
      phy = this_tree,
      GLS = TRUE
    )

    # keeping pc axes for 90% of variance, multiply by 1000 to make sure tiny values dont underflow in BayesTraits
    pcscores <- phypc$x[, c(
      1:which(cumsum(phypc$d / sum(phypc$d)) > 0.90)[1]
    )] *
      1000

    for (j in 1:length(list_of_models)) {
      for (i in 1:length(runs)) {
        # once for single rate models
        write.table(
          pcscores,
          file = paste0(
            outputfolder,
            "/pPCscores/Phylo_PC_SCORES_",
            "single_rate_tree",
            list_of_models[j],
            "_",
            runs[i],
            "_run",
            ".txt"
          ),
          quote = FALSE,
          col.names = FALSE
        )
        write_cmd_file(
          x = pcscores,
          model = list_of_models[j],
          run = runs[i],
          varrate = FALSE,
          folder = outputfolder
        )
        # and again for variable rates
        write.table(
          pcscores,
          file = paste0(
            outputfolder,
            "/pPCscores/Phylo_PC_SCORES_",
            "var_rate_tree",
            list_of_models[j],
            "_",
            runs[i],
            "_run",
            ".txt"
          ),
          quote = FALSE,
          col.names = FALSE
        )
        write_cmd_file(
          x = pcscores,
          model = list_of_models[j],
          run = runs[i],
          varrate = TRUE,
          folder = outputfolder
        )
      }
    }
  }
}
