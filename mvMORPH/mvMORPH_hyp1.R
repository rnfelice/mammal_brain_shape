# Script for comparing between-group rates with mvMORPH
# known data only
# tree_hyp1

library(phytools)
library(readr)
library(snow)
library(geomorph)
library(mvMORPH)
library(tidyverse)

# load phylogeny, trait data and landmark data
tree_hyp1 <- read.nexus('./phylogeny_construction/final_tree_hyp1.nex') # tree hypothesis 1
classifier <- read.csv('./input_data/mammal_data.csv')

# procrustes-aligned landmark data
n_landmarks <- 122
n_dimensions <- 3
raw_data <- read.csv("./input_data/procrutes_aligned_coords.csv", row.names = 1)
Y.gpa.endo <- arrayspecs(raw_data, p = n_landmarks, k = n_dimensions)

dataset <- two.d.array(Y.gpa.endo)


# Trim ecology dataset to known values

# Sociality

# remove taxa with unknown sociality
unknown_social <- classifier %>% filter(is.na(Sociality)) %>% pull(species)
known_social_data<- classifier %>% filter(!is.na(Sociality))
tree.social_hyp1 <- drop.tip(tree_hyp1, unknown_social)
Y.gpa.endo.social <- Y.gpa.endo[,,which(dimnames(Y.gpa.endo)[[3]] %in% known_social_data$species)]

writeNexus(tree.social_hyp1, file = "social_data_tree_hyp1.nex")
write.csv(known_social_data, file = "known_social_data.csv")
save(Y.gpa.endo.social, file = "Y.gpa.endo.social.R")


# Activity

# remove taxa with unknown sociality
unknown_activity <- classifier %>% filter(is.na(Activity)) %>% pull(species)
known_activity_data<- classifier %>% filter(!is.na(Activity))
tree.activity_hyp1 <- drop.tip(tree_hyp1, unknown_activity)
Y.gpa.endo.activity <- Y.gpa.endo[,,which(dimnames(Y.gpa.endo)[[3]] %in% known_activity_data$species)]

writeNexus(tree.activity_hyp1, file = "activity_data_tree_hyp1.nex")
write.csv(known_activity_data, file = "known_activity_data.csv")
save(Y.gpa.endo.activity, file = "Y.gpa.endo.activity.R")


# Locomotion

# remove taxa with unknown sociality
unknown_locomotion <- classifier %>% filter(is.na(Locomotion)) %>% pull(species)
known_locomotion_data<- classifier %>% filter(!is.na(Locomotion))
tree.locomotion_hyp1 <- drop.tip(tree_hyp1, unknown_locomotion)
Y.gpa.endo.locomotion <- Y.gpa.endo[,,which(dimnames(Y.gpa.endo)[[3]] %in% known_locomotion_data$species)]

writeNexus(tree.locomotion_hyp1, file = "locomotion_data_tree.nex")
write.csv(known_locomotion_data, file = "known_locomotion_data.csv")
save(Y.gpa.endo.locomotion, file = "Y.gpa.endo.locomotion.R")


# Diet

# remove taxa with unknown sociality
unknown_diet <- classifier %>% filter(is.na(Diet)) %>% pull(species)
known_diet_data<- classifier %>% filter(!is.na(Diet))
tree.diet_hyp1 <- drop.tip(tree_hyp1, unknown_diet)
Y.gpa.endo.diet <- Y.gpa.endo[,,which(dimnames(Y.gpa.endo)[[3]] %in% known_diet_data$species)]

writeNexus(tree.diet_hyp1, file = "diet_data_tree_hyp1.nex")
write.csv(known_diet_data, file = "known_diet_data.csv")
save(Y.gpa.endo.diet, file = "Y.gpa.endo.diet.R")


##########################
#### mvMORPH analysis ####
##########################

##### Activity #####

# Perform ancestral reconstruction

activity.known <- data.frame(known_activity_data$Activity)
rownames(activity.known) <- known_activity_data$species
cat_devo.known.activity_hyp1 = as.factor(activity.known[tree.activity_hyp1$tip.label,])
names(cat_devo.known.activity_hyp1) = tree.activity_hyp1$tip.label

activity.known <-as.vector(activity.known)

simmap.known.activity_hyp1 <- make.simmap(tree.activity_hyp1, cat_devo.known.activity_hyp1, 
                                     model="ARD", nsim=100)

# plotSimmap(tree_simmap, type="fan", fsize=0.1)


## -------------------------------------- ##
## Perform the analyses                   ##
## -------------------------------------- ##


dat.known.activity_hyp1 <- list(data=as.matrix(dataset[tree.activity_hyp1$tip.label,]*1000), 
                           mass=as.matrix(log(known_activity_data$Mass)))

dataset.known.activity <-two.d.array(Y.gpa.endo.activity)
mass.activity <-data.frame(known_activity_data$Mass)
rownames(mass.activity) = known_activity_data$species


# fit the model - note: remove the 'error' argument if you don't want to estimate a nuisance parameter
system.time(
  fit.known.activity_hyp1 <- mvgls(data~mass, data=dat.known.activity_hyp1, tree=simmap.known.activity_hyp1[[1]], model="BMM", MMSE=TRUE, method = "EmpBayes")
)

class(simmap.known.activity_hyp1)<-'list'
results.known.activity_hyp1 <- lapply(simmap.known.activity_hyp1, function(x) mvgls(data~mass, data=dat.known.activity_hyp1, tree=x, model="BMM", method = "EmpBayes", MMSE=TRUE, error = FALSE))

params.known.activity_hyp1 <- lapply(results.known.activity_hyp1, function(x) x$param)

combined_params.known.activity_hyp1 <- do.call("rbind",params.known.activity_hyp1)

combined_params_long.known.activity_hyp1 <- as_tibble(combined_params.known.activity_hyp1) %>% 
  pivot_longer(everything(), names_to = "Activity", values_to = "rate")

activity_known_hyp1 <- ggplot(combined_params_long.known.activity_hyp1,
                         aes(x=rate, col = Activity, fill = Activity))+
  scale_fill_manual(values = c(  "Mixed"           = "#9B59B6",
                                 "Diurnal"         = "#D7BDE2",
                                 "Nocturnal"       = "#5B2C6F"))+
  scale_colour_manual(values = c("Mixed"           = "#9B59B6",
                                 "Diurnal"         = "#D7BDE2",
                                 "Nocturnal"       = "#5B2C6F"))+
  geom_density(alpha=.7, size = 1)+
  theme_bw()


activity_known_hyp1 <- activity_known_hyp1 + guides(fill=guide_legend(title="Activity (known)"),
                                          col=guide_legend(title="Activity (known)"))+
  theme(legend.position = "inside")+
  theme(legend.position.inside = c(0.8, 0.7))
plot(activity_known_hyp1)



##### Sociality #####

social.known <- data.frame(known_social_data$Sociality)
rownames(social.known) <- known_social_data$species
cat_devo.known.social_hyp1 = as.factor(social.known[tree.social_hyp1$tip.label,])
names(cat_devo.known.social_hyp1) = tree.social_hyp1$tip.label

social.known <-as.vector(social.known)

simmap.known.social_hyp1 <- make.simmap(tree.social_hyp1, cat_devo.known.social_hyp1, 
                                   model="ARD", nsim=100)

fit.social_hyp1 <-fitMk(tree.social_hyp1,cat_devo.known.social_hyp1,model="ARD")

# plotSimmap(tree_simmap, type="fan", fsize=0.1)


## -------------------------------------- ##
## Perform the analyses                   ##
## -------------------------------------- ##

dataset.known.social <-two.d.array(Y.gpa.endo.social)
mass.social <-data.frame(known_social_data$Mass)
rownames(mass.social) = known_social_data$species


dat.known.social_hyp1 <- list(data=as.matrix(dataset[tree.social_hyp1$tip.label,]*1000), 
                         mass=as.matrix(log(known_social_data$Mass)))

# fit the model - note: remove the 'error' argument if you don't want to estimate a nuisance parameter
system.time(
  fit.known.social_hyp1 <- mvgls(data~mass, data=dat.known.social_hyp1, tree=simmap.known.social_hyp1[[1]], model="BMM", MMSE=TRUE, method = "EmpBayes")
)


class(simmap.known.social_hyp1)<-'list'
results.known.social_hyp1 <- lapply(simmap.known.social_hyp1, function(x) mvgls(data~mass, data=dat.known.social_hyp1, tree=x, model="BMM", method = "EmpBayes", MMSE=TRUE, error = FALSE))

params.known.social_hyp1 <- lapply(results.known.social_hyp1, function(x) x$param)

combined_params.known.social_hyp1 <-do.call("rbind",params.known.social_hyp1)

combined_params_long.known.social_hyp1 <- as_tibble(combined_params.known.social_hyp1) %>% 
  pivot_longer(everything(), names_to = "Sociality", values_to = "rate")

social_known_hyp1 <- ggplot(combined_params_long.known.social_hyp1,
                       aes(x=rate, col = Sociality, fill = Sociality))+
  scale_fill_manual(values = c(  "Social"          = "#1E8449",
                                 "Solitary"        = "#A9DFBF"))+
  scale_colour_manual(values = c("Social"          = "#1E8449",
                                 "Solitary"        = "#A9DFBF"))+
  geom_density(alpha=.7, size = 1)+
  theme_bw()

social_known_hyp1 <- social_known_hyp1 + guides(fill=guide_legend(title="Sociality (known)"),
                                      col=guide_legend(title="Sociality (known)"))+
  theme(legend.position = "inside")+
  theme(legend.position.inside = c(0.2, 0.7))
plot(social_known_hyp1)



##### Locomotion #####

locomotion.known <- data.frame(known_locomotion_data$Locomotion)
rownames(locomotion.known) <- known_locomotion_data$species
cat_devo.known.locomotion = as.factor(locomotion.known[tree.locomotion_hyp1$tip.label,])
names(cat_devo.known.locomotion) = tree.locomotion_hyp1$tip.label

locomotion.known <-as.vector(locomotion.known)

simmap.known.locomotion_hyp1 <- make.simmap(tree.locomotion_hyp1, cat_devo.known.locomotion, 
                                       model="ARD", nsim=100)

fit.locomotion_hyp1 <-fitMk(tree.locomotion_hyp1,cat_devo.known.locomotion,model="ARD")

# plotSimmap(tree_simmap, type="fan", fsize=0.1)


## -------------------------------------- ##
## Perform the analyses                   ##
## -------------------------------------- ##

dat.known.locomotion_hyp1 <- list(data=as.matrix(dataset[tree.locomotion_hyp1$tip.label,]*1000), 
                             mass=as.matrix(log(known_locomotion_data$Mass)))

dataset.known.locomotion <-two.d.array(Y.gpa.endo.locomotion)
mass.locomotion <-data.frame(known_locomotion_data$Mass)
rownames(mass.locomotion) = known_locomotion_data$species


# fit the model - note: remove the 'error' argument if you don't want to estimate a nuisance parameter
system.time(
  fit.known.locomotion_hyp1 <- mvgls(data~mass, data=dat.known.locomotion_hyp1, tree=simmap.known.locomotion_hyp1[[1]], model="BMM", MMSE=TRUE, method = "EmpBayes")
)


class(simmap.known.locomotion_hyp1)<-'list'
results.known.locomotion_hyp1 <- lapply(simmap.known.locomotion_hyp1, function(x) mvgls(data~mass, data=dat.known.locomotion_hyp1, tree=x, model="BMM", method = "EmpBayes", MMSE=TRUE, error = FALSE))

params.known.locomotion_hyp1 <- lapply(results.known.locomotion_hyp1, function(x) x$param)

combined_params.known.locomotion_hyp1 <-do.call("rbind",params.known.locomotion_hyp1)

combined_params_long.known.locomotion_hyp1 <- as_tibble(combined_params.known.locomotion_hyp1) %>% 
  pivot_longer(everything(), names_to = "Locomotion", values_to = "rate")

loco_known_hyp1 <- ggplot(combined_params_long.known.locomotion_hyp1,
                     aes(x=rate, col = Locomotion, fill = Locomotion))+
  scale_fill_manual(values = c(  "Terrestrial"     = "#D4E6F1",
                                 "Arboreal"        = "#1A5276",
                                 "Volant"          = "#A9CCE3",
                                 "Semi-arboreal"   = "#2E86C1",
                                 "Aquatic"         = "#1B2631",
                                 "Fossorial"       = "#85929E",
                                 "Semi-Aquatic"    = "#5D6D7E",
                                 "Semi-fossorial"  = "#AEB6BF"))+
  scale_colour_manual(values = c("Terrestrial"     = "#D4E6F1",
                                 "Arboreal"        = "#1A5276",
                                 "Volant"          = "#A9CCE3",
                                 "Semi-arboreal"   = "#2E86C1",
                                 "Aquatic"         = "#1B2631",
                                 "Fossorial"       = "#85929E",
                                 "Semi-Aquatic"    = "#5D6D7E",
                                 "Semi-fossorial"  = "#AEB6BF"))+
  geom_density(alpha=.7, size = 1)+
  theme_bw()

loco_known_hyp1 <- loco_known_hyp1 + guides(fill=guide_legend(title="Locomotion (known)"),
                                  col=guide_legend(title="Locomotion (known)"))+
  theme(legend.position = "inside")+
  theme(legend.position.inside = c(0.8, 0.6))
plot(loco_known_hyp1)



##### Diet #####

diet.known <- data.frame(known_diet_data$Diet)
rownames(diet.known) <- known_diet_data$species
cat_devo.known.diet_hyp1 = as.factor(diet.known[tree.diet_hyp1$tip.label,])
names(cat_devo.known.diet_hyp1) = tree.diet_hyp1$tip.label

diet.known <-as.vector(diet.known)

simmap.known.diet_hyp1 <- make.simmap(tree.diet_hyp1, cat_devo.known.diet_hyp1, 
                                 model="ARD", nsim=100)

# plotSimmap(tree_simmap, type="fan", fsize=0.1)


## -------------------------------------- ##
## Perform the analyses                   ##
## -------------------------------------- ##

dat.known.diet_hyp1 <- list(data=as.matrix(dataset[tree.diet_hyp1$tip.label,]*1000), 
                       mass=as.matrix(log(known_diet_data$Mass)))

dataset.known.diet <-two.d.array(Y.gpa.endo.diet)
mass.diet <-data.frame(known_diet_data$Mass)
rownames(mass.diet) = known_diet_data$species


# fit the model - note: remove the 'error' argument if you don't want to estimate a nuisance parameter
system.time(
  fit.known.diet_hyp1 <- mvgls(data~mass, data=dat.known.diet_hyp1, tree=simmap.known.diet_hyp1[[1]], model="BMM", MMSE=TRUE, method = "EmpBayes")
)


class(simmap.known.diet_hyp1)<-'list'
results.known.diet_hyp1 <- lapply(simmap.known.diet_hyp1, function(x) mvgls(data~mass, data=dat.known.diet_hyp1, tree=x, model="BMM", method = "EmpBayes", MMSE=TRUE, error = FALSE))

params.known.diet_hyp1 <- lapply(results.known.diet_hyp1, function(x) x$param)

combined_params.known.diet_hyp1 <-do.call("rbind",params.known.diet_hyp1)

combined_params_long.known.diet_hyp1 <- as_tibble(combined_params.known.diet_hyp1) %>% 
  pivot_longer(everything(), names_to = "Diet", values_to = "rate")

diet_known_hyp1 <- ggplot(combined_params_long.known.diet_hyp1,
                     aes(x=rate, col = Diet, fill = Diet))+
  scale_fill_manual(values = c(  "Omnivore"        = "#E67E22",
                                 "Carnivore"       = "#C0392B",
                                 "Herbivore"       = "#F9E79F"))+
  scale_colour_manual(values = c("Omnivore"        = "#E67E22",
                                 "Carnivore"       = "#C0392B",
                                 "Herbivore"       = "#F9E79F"))+
  geom_density(alpha=.7, size = 1)+
  theme_bw()

diet_known_hyp1 <- diet_known_hyp1 + 
                   guides(fill=guide_legend(title="Diet (known)"), col=guide_legend(title="Diet (known)"))+
                   theme(legend.position = "inside")+
                   theme(legend.position.inside = c(0.2, 0.7))
plot(diet_known_hyp1)


# write.csv(combined_params.known.diet, file = "~/Dropbox/_UCL/Projects/Team Tetrapod/mvmorph/EmpBay_devo_rates_NO_error.csv")
# saveRDS(results, "~/Dropbox/_UCL/Projects/Team Tetrapod/mvmorph/EmpBay_devo_rates_NO_error_full_results.RDS")

ggarrange(activity_known_hyp1,social_known_hyp1, loco_known_hyp1, diet_known_hyp1, 
          ncol=2, nrow=2, common.legend = F)


col.mvmorph <- values = c(
  # Set 1 – Purples (activity pattern)
  "Mixed"           = "#9B59B6",
  "Diurnal"         = "#D7BDE2",
  "Nocturnal"       = "#5B2C6F",
  
  # Set 2 – Reds/Oranges (diet)
  "Omnivore"        = "#E67E22",
  "Carnivore"       = "#C0392B",
  "Herbivore"       = "#F9E79F",
  
  # Set 3 – Greens (sociality)
  "Social"          = "#1E8449",
  "Solitary"        = "#A9DFBF",
  
  # Set 4 – Blues (habitat)
  "Terrestrial"     = "#D4E6F1",
  "Arboreal"        = "#1A5276",
  "Volant"          = "#A9CCE3",
  "Semi-arboreal"   = "#2E86C1",
  "Aquatic"         = "#1B2631",
  "Fossorial"       = "#85929E",
  "Semi-Aquatic"    = "#5D6D7E",
  "Semi-fossorial"  = "#AEB6BF"
),  na.value = "transparent"