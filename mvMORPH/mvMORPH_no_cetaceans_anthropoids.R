# Perform mvMORPH analysis after removing cetaceans and anthropoids from dataset
# tree_hyp1

library(phytools)
library(readr)
library(snow)
library(geomorph)
library(mvMORPH)
library(tidyverse)

# load phylogeny, trait data and landmark data
tree_hyp1 <- read.nexus("./phylogeny_construction/final_tree_hyp1.nex") # tree hypothesis 1
classifier <- read.csv("./input_data/mammal_data.csv")
Y.gpa.endo <- load("./input_data/Y.gpa.endo.R") # procrustes-aligned landmark data
dataset <- two.d.array(Y.gpa.endo)

# Remove cetaceans and anthropoids from tree ('trim' indicates dataset with these clades removed)
# Find MRCA of each clade

# Cetaceans
findMRCA(tree_hyp1, c("Orcaella_heinsohni", "Protocetus_atavus"))
plot(extract.clade(tree_hyp1, 665), cex = 0.6)
trim_cetaceans <- (extract.clade(tree_hyp1, 665))

# Anthropoids
findMRCA(tree_hyp1, c("Parapithecus_grangeri", "Cercopithecus_ascanius"))
plot(extract.clade(tree_hyp1, 888), cex = 0.6)
trim_simians <- (extract.clade(tree_hyp1, 888))

# create vector of all simians and cetaceans to drop
trim_list <- c(trim_cetaceans$tip.label, trim_simians$tip.label)

# Drop tips of cetaceans and simians
tree_hyp1_trimmed <- drop.tip(tree_hyp1, tip = trim_list)

plot(tree_hyp1_trimmed, cex = 0.3, type = "fan")


# Edit classifier to removed trimmed taxa
classifier_trim <- classifier %>% filter(!(species %in% trim_list))

# check new data frame
View(classifier_trim)


# Perform mvMORPH analyses

# Sociality

# remove taxa with unknown sociality
unknown_social_trim <- classifier_trim %>%
  filter(is.na(Sociality)) %>%
  pull(species)
known_social_data_trim <- classifier_trim %>% filter(!is.na(Sociality))
tree.social.trim_hyp1 <- drop.tip(tree_hyp1_trimmed, unknown_social)
Y.gpa.endo.social.trim <- Y.gpa.endo[, , which(dimnames(Y.gpa.endo)[[3]] %in% known_social_data_trim$species)]
dim(Y.gpa.endo.social.trim)

writeNexus(tree.social.trim_hyp1, file = "social_data_trim_tree_hyp1.nex")
write.csv(known_social_data_trim, file = "known_social_data_trim.csv")
save(Y.gpa.endo.social.trim, file = "Y.gpa.endo.social.trim.R")


##########################
#### mvMORPH analysis ####
##########################

# Known data only
# Remove unknowns from dataset

##### Sociality - known data only #####

social.known.trim <- data.frame(known_social_data_trim$Sociality)
rownames(social.known.trim) <- known_social_data_trim$species
cat_devo.known.social.trim <- as.factor(social.known.trim[tree.social.trim_hyp1$tip.label, ])
names(cat_devo.known.social.trim) <- tree.social.trim_hyp1$tip.label

social.known.trim <- as.vector(social.known.trim)

simmap.known.social.trim <- make.simmap(tree.social.trim_hyp1, cat_devo.known.social.trim,
  model = "ARD", nsim = 100
)

fit.known.social.trim <- fitMk(tree.social.trim_hyp1, cat_devo.known.social.trim, model = "ARD")

# plotSimmap(tree_simmap, type="fan", fsize=0.1)


## -------------------------------------- ##
## Perform the analyses                   ##
## -------------------------------------- ##


dataset.known.social.trim <- two.d.array(Y.gpa.endo.social.trim)
mass.known.social.trim <- data.frame(known_social_data_trim$Mass)
rownames(mass.known.social.trim) <- known_social_data_trim$species

dat.known.social.trim <- list(
  data = as.matrix(dataset[tree.social.trim_hyp1$tip.label, ] * 1000),
  mass = as.matrix(log(known_social_data_trim$Mass))
)

# fit the model - note: remove the 'error' argument if you don't want to estimate a nuisance parameter
system.time(
  fit.known.social.trim <- mvgls(data ~ mass, data = dat.known.social.trim, tree = simmap.known.social.trim[[1]], model = "BMM", MMSE = TRUE, method = "EmpBayes")
)

class(simmap.known.social.trim) <- "list"
results.known.social.trim <- lapply(simmap.known.social.trim, function(x) mvgls(data ~ mass, data = dat.known.social.trim, tree = x, model = "BMM", method = "EmpBayes", MMSE = TRUE, error = FALSE))

params.known.social.trim <- lapply(results.known.social.trim, function(x) x$param)

combined_params.known.social.trim <- do.call("rbind", params.known.social.trim)

combined_params_long.known.social.trim <- as_tibble(combined_params.known.social.trim) %>%
  pivot_longer(everything(), names_to = "Sociality", values_to = "rate")

social_known.trim <- ggplot(
  combined_params_long.known.social.trim,
  aes(x = rate, col = Sociality, fill = Sociality)
) +
  scale_fill_manual(values = c(
    "Social" = "#1E8449",
    "Solitary" = "#A9DFBF"
  )) +
  scale_colour_manual(values = c(
    "Social" = "#1E8449",
    "Solitary" = "#A9DFBF"
  )) +
  geom_density(alpha = .7) +
  theme_bw()

social_known.trim <- social_known.trim + guides(
  fill = guide_legend(title = "Sociality (known)"),
  col = guide_legend(title = "Sociality (known)")
)
plot(social_known.trim)


##### Locomotion - known data only #####

# Sociality

# remove taxa with unknown sociality
unknown_locomotion_trim <- classifier_trim %>%
  filter(is.na(Locomotion)) %>%
  pull(species)
known_locomotion_data_trim <- classifier_trim %>% filter(!is.na(Locomotion))
tree.locomotion.trim_hyp1 <- drop.tip(tree_hyp1_trimmed, unknown_locomotion)
Y.gpa.endo.locomotion.trim <- Y.gpa.endo[, , which(dimnames(Y.gpa.endo)[[3]] %in% known_locomotion_data_trim$species)]


writeNexus(tree.locomotion.trim_hyp1, file = "locomotion_data_trim_tree_hyp1.nex")
write.csv(known_locomotion_data_trim, file = "known_locomotion_data_trim.csv")
save(Y.gpa.endo.locomotion.trim, file = "Y.gpa.endo.locomotion.trim.R")

locomotion.known.trim <- data.frame(known_locomotion_data_trim$Locomotion)
rownames(locomotion.known.trim) <- known_locomotion_data_trim$species
cat_devo.known.locomotion.trim <- as.factor(locomotion.known.trim[tree.locomotion.trim_hyp1$tip.label, ])
names(cat_devo.known.locomotion.trim) <- tree.locomotion.trim_hyp1$tip.label

locomotion.known.trim <- as.vector(locomotion.known.trim)

simmap.known.locomotion.trim <- make.simmap(tree.locomotion.trim_hyp1, cat_devo.known.locomotion.trim,
  model = "ARD", nsim = 100
)

fit.known.locomotion.trim <- fitMk(tree.locomotion.trim_hyp1, cat_devo.known.locomotion.trim, model = "ARD")

# plotSimmap(tree_simmap, type="fan", fsize=0.1)


## -------------------------------------- ##
## Perform the analyses                   ##
## -------------------------------------- ##


dataset.known.locomotion.trim <- two.d.array(Y.gpa.endo.locomotion.trim)
mass.known.locomotion.trim <- data.frame(known_locomotion_data_trim$Mass)
rownames(mass.known.locomotion.trim) <- known_locomotion_data_trim$species

dat.known.locomotion.trim <- list(
  data = as.matrix(dataset[tree.locomotion.trim_hyp1$tip.label, ] * 1000),
  mass = as.matrix(log(known_locomotion_data_trim$Mass))
)

# fit the model - note: remove the 'error' argument if you don't want to estimate a nuisance parameter
system.time(
  fit.known.locomotion.trim <- mvgls(data ~ mass, data = dat.known.locomotion.trim, tree = simmap.known.locomotion.trim[[1]], model = "BMM", MMSE = TRUE, method = "EmpBayes")
)

class(simmap.known.locomotion.trim) <- "list"
results.known.locomotion.trim <- lapply(simmap.known.locomotion.trim, function(x) mvgls(data ~ mass, data = dat.known.locomotion.trim, tree = x, model = "BMM", method = "EmpBayes", MMSE = TRUE, error = FALSE))

params.known.locomotion.trim <- lapply(results.known.locomotion.trim, function(x) x$param)

combined_params.known.locomotion.trim <- do.call("rbind", params.known.locomotion.trim)

combined_params_long.known.locomotion.trim <- as_tibble(combined_params.known.locomotion.trim) %>%
  pivot_longer(everything(), names_to = "locomotionity", values_to = "rate")

locomotion_known.trim <- ggplot(
  combined_params_long.known.locomotion.trim,
  aes(x = rate, col = locomotionity, fill = locomotionity)
) +
  scale_fill_manual(values = c(
    "Terrestrial" = "#D4E6F1",
    "Arboreal" = "#1A5276",
    "Volant" = "#A9CCE3",
    "Semi-arboreal" = "#2E86C1",
    "Aquatic" = "#1B2631",
    "Fossorial" = "#85929E",
    "Semi-Aquatic" = "#5D6D7E",
    "Semi-fossorial" = "#AEB6BF"
  )) +
  scale_colour_manual(values = c(
    "Terrestrial" = "#D4E6F1",
    "Arboreal" = "#1A5276",
    "Volant" = "#A9CCE3",
    "Semi-arboreal" = "#2E86C1",
    "Aquatic" = "#1B2631",
    "Fossorial" = "#85929E",
    "Semi-Aquatic" = "#5D6D7E",
    "Semi-fossorial" = "#AEB6BF"
  )) +
  geom_density(alpha = .7) +
  theme_bw()

locomotion_known.trim <- locomotion_known.trim + guides(
  fill = guide_legend(title = "locomotionity (known)"),
  col = guide_legend(title = "locomotionity (known)")
)
plot(locomotion_known.trim)
