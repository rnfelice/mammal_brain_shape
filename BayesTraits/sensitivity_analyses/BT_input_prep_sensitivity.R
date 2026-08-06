# Prepping files for BayesTraits Rates Analysis ----------------------------------------------
library(ape)
library(phytools)
library(paleotree)
library(tidyverse)
library(geomorph)
library(geiger)
# devtools::install_github("joannabaker/BayesTraitR")
library(BayesTraitR)

#load shape data
n_landmarks <- 122
n_dimensions <- 3
raw_data <- read.csv("./input_data/procrutes_aligned_coords.csv", row.names = 1)
Y.gpa.endo <- arrayspecs(raw_data, p = n_landmarks, k = n_dimensions)

#load size data
Csize <- read.csv("./input_data/centroid_sizes.csv", row.names = 1)

data <- read.csv("./input_data/mammal_data.csv", row.names = 1)
# confirm same order
data <- data[rownames(Csize), ]
bodysize <- data[, "Mass", drop = F]
clade <- data[, "order", drop = F]
regression_data <- data.frame(spp.names = rownames(Csize), lnCsize = log(Csize[, 1]), lnMass = log(bodysize[, 1]), lnMassSq = log(bodysize[, 1])^2, order = clade)

# define modules
module_list <- list(
#for sensitivity analysis we will focus on whole-brain only
  whole_brain = c(1:122)) 

#definte function to write bayestraits comman files:
#the function then sets the prior for alpha to be a uniform distribution from the min to the max of that PC score (trait)
write_cmd_file <- function(x, varrate = TRUE, run, model = "BM", folder) {
  filename_cmd <- paste0(
    folder,
    "/command_files/BT_Control_",
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


tree_files<-list.files("./BayesTraits/sensitivity_analyses/Alternative_Trees/",
pattern=".nex")
tree_names<-gsub(".*tree_|\\.nex.*", "", tree_files)

for (treenum in 1:length(tree_files)) {
  this_tree <- read.nexus(
    file = paste0("./BayesTraits/sensitivity_analyses/Alternative_Trees/", tree_files[treenum])
  )
  #drop tips not in dataset
  this_tree<-drop.tip(this_tree, name.check(phy=this_tree,data=two.d.array(Y.gpa.endo))$tree_not_data)

  
  for (modulenum in 1:length(module_list)) {
    #define output folder for this tree and brain region
    outputfolder <- paste0(
      "./BayesTraits/sensitivity_analyses/BT_input/",
      names(module_list)[modulenum],
      "/",
      tree_names[treenum]
    )
    
    write.nexus(this_tree, file=paste0(outputfolder, "/trees/", "tree_",tree_names[treenum],".nex" ))
    #do phylogenetic PCA
    phypc <- gm.prcomp(
      A = Y.gpa.endo[module_list[[modulenum]], , ],
      phy = this_tree,
      GLS = TRUE
    )

    #keeping pc axes for 90% of variance, multiply by 1000 to make sure tiny values dont underflow in BayesTraits
    pcscores <- phypc$x[, c(
      1:which(cumsum(phypc$d / sum(phypc$d)) > 0.90)[1]
    )] *
      1000

    for (j in 1:length(list_of_models)) {
      for (i in 1:length(runs)) {
        #once for single rate models
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
        #and again for variable rates
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


#create body size analysis scripts
 bodysize_folder <- paste0(
  "BayesTraits/sensitivity_analyses/BT_input/brainsize/",
  tree_names[treenum]
)
  
  
createBTjob(
  cols = c("lnCsize", "lnMass"),
  dataset = regression_data,
  tree = this_tree,
  jobname = "Brain_Body_BT",
  bi = 100000000,
  it = 1000000000,
  sa = 100000,
  model = 9,
  MCMC = T,
  reps = 2,
  optarg = c(
      "varrates",
      "stones 500 10000",
      "RJLockModel",
      # following Baker et al 2025, we use a "wide and uninformative normal prior centred on zero with a standard deviation of 2.5 for all regression parameters"
      BTBetapriors(data = regression_data, prior = "normal", pars = c(0, 2.5))
  ),
  outdir = paste0(bodysize_folder, "/linear"),
  fm = as.formula(lnCsize ~ lnMass)
)
#linear regression+lambda
createBTjob(
  cols = c("lnCsize", "lnMass"),
  dataset = regression_data,
  tree = this_tree,
  jobname = "Brain_Body_BT",
  bi = 100000000,
  it = 1000000000,
  sa = 100000,
  model = 9,
  MCMC = T,
  reps = 2,
  optarg = c(
      "varrates",
      "lambda",
      "stones 500 10000",
      "RJLockModel",
      # following Baker et al 2025, we use a "wide and uninformative normal prior centred on zero with a standard deviation of 2.5 for all regression parameters"
      BTBetapriors(data = regression_data, prior = "normal", pars = c(0, 2.5))
  ),
  outdir = outdir = paste0(bodysize_folder, "/lambda"),
  fm = as.formula(lnCsize ~ lnMass)
)
  
  #quadratic regression
  
  createBTjob(
    cols = c("lnCsize", "lnMass", "lnMassSq"),
    dataset = regression_data,
    tree = this_tree,
    jobname = "Brain_Body_BT",
    bi = 100000000,
    it = 1000000000,
    sa = 100000,
    model = 9,
    MCMC = T,
    reps = 2,
    optarg = c(
        "varrates",
        "stones 500 10000",
        "RJLockModel",
        # following Baker et al 2025, we use a "wide and uninformative normal prior centred on zero with a standard deviation of 2.5 for all regression parameters"
        BTBetapriors(data = regression_data, prior = "normal", pars = c(0, 2.5))
    ),
    outdir = outdir = paste0(bodysize_folder, "/quadratic"),
    fm = as.formula(lnCsize ~ lnMass + lnMassSq)
)

# now repeat to do BayesTraits analysis with order as a covariate
createBTjob(
    cols = c("lnCsize", "lnMass", "order"),
    dataset = regression_data,
    tree = this_tree,
    jobname = "Brain_Body_BT",
    bi = 100000000,
    it = 1000000000,
    sa = 100000,
    model = 9,
    MCMC = T,
    reps = 2,
    optarg = c(
        "varrates",
        "stones 500 10000",
        "RJLockModel",
        # following Baker et al 2025, we use a "wide and uninformative normal prior centred on zero with a standard deviation of 2.5 for all regression parameters"
        BTBetapriors(data = regression_data, prior = "normal", pars = c(0, 2.5))
    ),
    outdir = outdir = paste0(bodysize_folder, "/clade_specific_slopes"),
    fm = as.formula(lnCsize ~ lnMass + order)
)

}

####################
#Brain size analyses
####################






# read tree

