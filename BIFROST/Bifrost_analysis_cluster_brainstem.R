#bifrost test


library(bifrost)
library(ape)
library(phytools)
library(mvMORPH)
library(geomorph)

cat("=== Parallelism diagnostics ===\n")
cat("SLURM_CPUS_PER_TASK env var:", Sys.getenv("SLURM_CPUS_PER_TASK", "not set"), "\n")
cat("parallelly::availableCores():", parallelly::availableCores(), "\n")
cat("future::nbrOfWorkers() (pre-plan, expect 1):", future::nbrOfWorkers(), "\n")
cat("===============================\n")

set.seed(1999)


# load phylogeny, trait data and landmark data
tree_hyp1 <- read.nexus('/home/rfelice/me/bifrost/data/final_tree_hyp1.nex') # tree hypothesis 1
classifier <- read.csv('/home/rfelice/me/bifrost/data/mammal_data.csv')

# procrustes-aligned landmark data
n_landmarks <- 122
n_dimensions <- 3
raw_data <- read.csv("/home/rfelice/me/bifrost/data/procrutes_aligned_coords.csv", row.names = 1)

data<-raw_data[tree_hyp1$tip.label,]
geiger::name.check(tree_hyp1, data)


base <- phytools::paintBranches(
  tree_hyp1,
  edge      = unique(tree_hyp1$edge[, 2]),
  state     = "0",
  anc.state = "0"
)

identical(rownames(data), base$tip.label)
library(future)

options(future.globals.maxSize = 10 * 1024^3)

{
    n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
    Sys.setenv(
      OMP_NUM_THREADS        = "1",
      OPENBLAS_NUM_THREADS   = "1",
      MKL_NUM_THREADS        = "1",
      VECLIB_MAXIMUM_THREADS = "1"
    )
}

# no future.batchtools / plan() needed — bifrost sets its own plan via num_cores

res <- searchOptimalConfiguration(
  baseline_tree              = base,
  trait_data                 = as.matrix(data),
  formula                    = "trait_data ~ 1",
  min_descendant_tips        = 4,
  num_cores                  = n_cores,
  shift_acceptance_threshold = 20,  # conservative GIC threshold
  IC                         = "GIC",
  plot                       = TRUE,
  store_model_fit_history    = FALSE,
  method = "EmpBayes"
)

save(res,
 file="/home/rfelice/me/bifrost/data/bifrost_res_july.rda")
