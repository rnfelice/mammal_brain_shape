# this script creates the additional xml code BEAST2 tip dating analysis
# this includes  tip date priors and node age constraints
# at the bottom of this script there are instructions for where to paste the code into the main BEAST xml file which should be created in BEAUTi


library(ape)
library(paleotree)
library(geiger)
library(tidyverse)
library(RRphylo)

# load the trees and tip date data:

# this is the time-calibrated mean tree from Alvarez-Carretero et al
time_tree <- read.nexus(file = "./input_data/mammaltree4705_editednames.nex")

# this is the tree with uncalibrated node dates that contains all
grafted_tree <- read.nexus(file = "./input_data/merged_tree.nex")

# import species data
species_data <- read.csv("./input_data/mammal_data.csv", fileEncoding = "latin1", row.names = 1)

# write down tip occurence ranges:
occurences_in_sample <- species_data[, c("age_max", "age_min")]

nc <- name.check(phy = grafted_tree, data = occurences_in_sample)
nc


# some taxa are in the tree for which we lack shape data.
# we still want them for tip dating
# give them fixed tip ages at present time:
other_extant_taxa <- matrix(0, nrow = length(nc$tree_not_data), ncol = 2)
rownames(other_extant_taxa) <- nc$tree_not_data
colnames(other_extant_taxa) <- c("age_max", "age_min")
fads_and_lads <- rbind(occurences_in_sample, other_extant_taxa)
name.check(phy = grafted_tree, data = fads_and_lads)

# write this full table out in a format that BEAUTi likes. We will use this to
# tell BEAUTI which tips are extant and extinct by
mean_tipdates <- as.matrix(rowMeans(fads_and_lads))
write.table(mean_tipdates,
  file = "./phylogeny_construction/BEAST2/tipdates.dat",
  sep = "\t",
  col.names = FALSE,
  quote = FALSE
)


# before we do anything else, choose tha name that you will give your tree in BEAUTi later.

treename <- "mammal_tree"
##########
# for each FOSSIL tip, make a uniform prior distibution for its age.
# note that the "taxonset id" must be unique (different from the name of the tip label itself)
# so we add the word "tip" to that.
##########

ages <- subset(fads_and_lads, !(age_max == 0 & age_min == 0))

make_tip_date_priors <- function(species, fad, lad) {
  paste0(
    '                <distribution id="', species, '_tip.prior" spec="sa.math.distributions.SAMRCAPrior" tipsonly="true" tree="@Tree.t:', treename, '">\n',
    '                    <taxonset id="', species, '_tip" spec="TaxonSet">\n',
    '                        <taxon idref="', species, '" spec="Taxon"/>\n',
    "                    </taxonset>\n",
    '                    <Uniform id="Uniform.', species, '" lower="', lad, '" name="distr" upper="', fad, '"/>\n',
    "                </distribution>\n"
  )
}


tip_date_text <- mapply(
  make_tip_date_priors,
  rownames(ages),
  ages$age_max,
  ages$age_min,
  SIMPLIFY = FALSE
)

# write this out to a file, then paste into the main BEAST xml file after the node contraints
writeLines(unlist(tip_date_text), "./phylogeny_construction/BEAST2/ages.xml")

##########
# for each tip date prior need an operator that allows the MCMC to randomly sample the distribution
##########

make_tip_operator <- function(taxon) {
  paste0(
    '        <operator id="tipDatesSampler.', taxon,
    '_tip" spec="sa.evolution.operators.SampledNodeDateRandomWalker" taxonset="@',
    taxon,
    '_tip" tree="@Tree.t:', treename, '" weight="1.0" windowSize="1.0"/>'
  )
}

tip_operators <- sapply(rownames(ages), make_tip_operator)


# write out to a file, then paste this into the main BEAST xml file at the end of
# the operators section, just before the logger section
writeLines(tip_operators, "./phylogeny_construction/BEAST2/tip_date_operators.xml")


############
# Node constraints and priors
############


# first find all the most recent common ancestors with known dates
time_tree_splits <- prop.part(time_tree)
time_tree_splits <- time_tree_splits[sapply(time_tree_splits, length) != Ntip(time_tree)]
time_tree_clades <- lapply(time_tree_splits, function(x) time_tree$tip.label[x])

# first corresponding nodes in the undated supertree
grafted_tree_splits <- prop.part(grafted_tree)
grafted_tree_splits <- grafted_tree_splits[sapply(grafted_tree_splits, length) != Ntip(grafted_tree)]
grafted_tree_clades <- lapply(grafted_tree_splits, function(x) grafted_tree$tip.label[x])


# get node ages from Alvarez-Carratero time tree
time_tree$root.time <- max(vcv.phylo(time_tree))
time_tree_dates <- dateNodes(time_tree)

node_date_data <- tibble(
  node_number = double(),
  tiplabs = character()
)

for (i in 1:length(grafted_tree_clades)) {
  this_anc <- getMRCA(phy = grafted_tree, tip = grafted_tree_clades[[i]])
  subclade <- extract.clade(grafted_tree, this_anc)
  subclade_tips <- paste(subclade$tip.label, collapse = " ")
  node_info <- data.frame(node_number = this_anc, tiplabs = subclade_tips)
  node_date_data <- bind_rows(node_date_data, node_info)
}

node_date_data <- node_date_data %>%
  mutate(
    related_clade_in_timetree = NA,
    calibration_date = NA
  )

for (i in 1:length(time_tree_clades)) {
  this_anc <- getMRCA(phy = time_tree, tip = time_tree_clades[[i]])
  this_anc_age <- time_tree_dates[this_anc]
  node_in_grafted_tree <- getMRCA(phy = grafted_tree, tip = time_tree_clades[[i]])
  subclade <- extract.clade(grafted_tree, node_in_grafted_tree)
  subclade_tips <- paste(subclade$tip.label, collapse = " ")
  if (!any(fads_and_lads[subclade$tip.label, "age_max"] > this_anc_age)) {
    this_anc_age <- this_anc_age
  } else {
    this_anc_age <- NA
  }
  node_date_data$calibration_date[which(node_date_data$node_number == node_in_grafted_tree)] <- this_anc_age
  node_date_data$related_clade_in_timetree[which(node_date_data$node_number == node_in_grafted_tree)] <- this_anc
}


make_constraint_xml <- function(tip_list, node_number, calibration_date) {
  # Split taxa on spaces
  taxa <- strsplit(tip_list, "\\s+")[[1]]

  # Build taxon lines
  taxon_lines <- paste0(
    '                        <taxon idref="', taxa, '"/>',
    collapse = "\n"
  )
  # Build node date lines(for nodes with known dates)
  if (!is.na(calibration_date)) {
    node_date_constraint <- paste0(
      '                    <Uniform id="Uniform.', node_number, '_node" lower="', calibration_date * .95, '" name="distr" upper="', calibration_date * 1.05, '"/>\n'
    )
  } else {
    node_date_constraint <- ""
  }
  # Assemble XML
  paste0(
    '                <distribution id="', node_number, '_node.prior" spec="beast.base.evolution.tree.MRCAPrior" monophyletic="true" tree="@Tree.t:', treename, '">\n',
    '                    <taxonset id="', node_number, '_node" spec="TaxonSet">\n',
    taxon_lines, "\n",
    "                    </taxonset>\n",
    node_date_constraint,
    "                </distribution>\n"
  )
}
node_constraints_xml <- mapply(
  make_constraint_xml,
  node_date_data$tiplabs,
  node_date_data$node_number,
  node_date_data$calibration_date,
  SIMPLIFY = FALSE
)
writeLines(unlist(node_constraints_xml), "./phylogeny_construction/BEAST2/constraints.xml")

# make a constraint for the toot of the whole tree!

root_constraints_xml <-
  make_constraint_xml(
    paste0(grafted_tree$tip.label, collapse = " "),
    0,
    202,
  )
writeLines(root_constraints_xml, "./phylogeny_construction/BEAST2/root.xml")


# add loggers for the tip dates and node ages
log_lines <- paste0(
  '        <log idref="',
  rownames(ages),
  '.prior"/>'
)
writeLines(log_lines, "./phylogeny_construction/BEAST2/prior_logs.xml")

node_log_lines <- paste0(
  '        <log idref="',
  node_date_data$node_number,
  '_node.prior"/>'
)
writeLines(node_log_lines, "./phylogeny_construction/BEAST2/node_prior_logs.xml")


# The taxon sets all use the "idref" notation. we need to make sure that all the
# tips are known to BEAST as taxa before hand, we can do that with a plate like this:
# note that unlike above we do not specify monophyly

all_taxa_set <- paste0(
  '<taxonset id="alltips" spec="TaxonSet">\n',
  '  <plate var="n" range="', paste(grafted_tree$tip.label, collapse = ","), '">\n',
  '  <taxon id="$(n)" spec="Taxon"/>\n',
  "  </plate>\n",
  "  </taxonset>\n"
)

writeLines(all_taxa_set, "./phylogeny_construction/BEAST2/all_taxa_set.xml")

# export the dummy trait data for BEAUTi.
trait_data <- as.matrix(rep("?", length(grafted_tree$tip.label)))
rownames(trait_data) <- grafted_tree$tip.label
write.nexus.data(trait_data,
  format = "standard",
  file = "./phylogeny_construction/BEAST2/beauti_morphological_data.nex"
)
# now go edit that file and add the following:
#  CHARSTATELABELS
# 1  dummy / pesent. absent.,
# ;
# after the "FORMAT" line and before the "MATRIX" line
# (BEAUTi only accepts morphological data with char state labels)

##############
# now go to BEAUTi to create the template for the BEAST xml file.
# we did this with Beast 2.7.8
# and the following packages:
# BEAST_CLASSIC 1.6.4
# BEASTLabs 2.0.3
# CCD 1.0.3
# FixedTreeAnalysis 0.0.2 (might not be necessary but we did have it installed)
# GEO_SPHERE 1.4.2 (might not be necessary but we did have it installed)
# MM 1.2.1
# MODEL_SELECTION 1.6.2
# NS 1.2.2
# ORC 1.2.1
# SA 2.1.1
# starbeast3 1.2.1
# timtam 0.4.0

# in beast, choose File->Add Morphological Data and choose your beauti_morphological_data.nex
# it will ask you whether to condition on variable characters, choose yes
# rename the site model, trait model, and tree model to be a bit more readable
# important: make the name of your tree matches what you picked for "treename" above
# we chose "traits" "clock" and "mammal_tree" respectively
# go to the it Tip Dates tab. check "Use tip dates"
# Click Auto-configure, choose "read from file"
# import tipdats.dat. back in the Tip Dates menu, use the dropdown menu to have dates specified as years before the present
# in the Site Model it should have Gamma Site model with Substitution rate estimated and Subst Model = LewisMK
# in the Clock Model tab, choose Optimised Relaxed Clock
# in the priors menu, switch the tree prior to Fossilized Birth Death Model. Click the down arrow to update the options:
# Origin should be an estimate of the root age- for us this is 202Ma based on Alvarez Carretero's mean tree
# make sure estimate is ticked for this.
# Set Rho (the sampling probablity of extant tips) and make sure "estimated" is UNTICKED
# alvarez-carratero et al has 4705 extant taxa. there are 6871 extant taxa currenty recognized.
# so the sampling prob must be about 4705/6871 = 0.685
# going back to the Origin parameter, we will set priors on it: go down to the originFBD.t prior, click the down arrow
# there are lots of ways to parameterize this- we found the best results with an exponential offset prior with an offset of 188, mean = 85, standard deviation = 5
# set the sampling ProportionFBD prior to Beta with Alpha = 1 and Beta =1
# this is a flat prior and follows Matzke and Wright 2016 and the paleotree package
# for diversification rate, https://taming-the-beast.org/tutorials/FBD-tutorial/ recommends an exponential distribution
# we have chosen exponential distribution with mean=2 to be slightly more vague
# turnoverFBD will be left as an uniformative Uniform(0,1) prior
# finally we go to the MCMC tab and tick the "sample from prior" button and save the xml file.
# then open that file and paste in all the chunks of code we made above and execute in BEAST!

# you also need to add this operator with your tree mame in it
#<operator spec='sa.evolution.operators.SampledNodeDateRandomWalker' windowSize="1"  tree="@Tree.t:bears" weight="10">
