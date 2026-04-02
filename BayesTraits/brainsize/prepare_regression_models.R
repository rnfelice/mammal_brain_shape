# devtools::install_github("joannabaker/BayesTraitR")
library(BayesTraitR)

setwd("~/Library/CloudStorage/OneDrive-UniversityCollegeLondon/Project/data/Mammals/analysis/BayesTraits/")

Csize <- read.csv("./input_data/centroid_sizes.csv", row.names = 1)

data <- read.csv("../input_data/mammal_data.csv", row.names = 1)
# confirm same order
data <- data[rownames(Csize), ]
bodysize <- data[, "Mass", drop = F]
clade <- data[, "order", drop = F]

# read tree
tree1 <- read.nexus("./whole_brain/hyp_1/trees/final_tree_hyp1.nex")
regression_data <- data.frame(spp.names = rownames(Csize), lnCsize = log(Csize[, 1]), lnMass = log(bodysize[, 1]), lnMassSq = log(bodysize[, 1])^2, order = clade)

createBTjob(
    cols = c("lnCsize", "lnMass", "lnMassSq"),
    dataset = regression_data,
    tree = tree1,
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
    outdir = "./brainsize/quadratic",
    fm = as.formula(lnCsize ~ lnMass + lnMassSq)
)

# now repeat to do BayesTraits analysis with order as a covariate
createBTjob(
    cols = c("lnCsize", "lnMass", "order"),
    dataset = regression_data,
    tree = tree1,
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
    outdir = "./brainsize/clade_specific_slopes",
    fm = as.formula(lnCsize ~ lnMass + order)
)
