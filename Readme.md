# Data and Code from Knapp et al.
This repository contains all data and code for carrying out the analyses in _Title of the paper goes here_. 

## Folders:

📁 BayesTraits: this folder contains the script `BT_input_prep.R`. This script is used to prepare the inputs for BayesTraits evolutionary rate analyses has everal functions: (1) Carry out phylogenetic PCA to reduce the dimensionality of the Procrustes-aligned landmark data (2) create BayesTraits command file (3) organise these input files into folders to analyses the evolutiooanry rates of the whole endocasts and each subregion (cerebellum, neocortex, olfactory lobe and brainstem).

Other contents of this folder is the subfolder 📁 brainsize which does the same task for preparing inputs for the variable-rates phylogenetic regression in BayesTriats to analyse brain size allometry evolution. Finally, the R script `compare_bayes_factors.R` is used to compare the marginal likelihoods of BayesTraits results for model selection. 

---

📁 endocast_pts: this folder contains the raw (i.e., not Procrustes-aligned) landmark and semilandmark data for all specimens

---

📁 phylogeny_construction: this folder includes the scripts to generate the eight dated supertrees used in the analyeses

    📦phylogeny_construction
    ┣ 📂BEAST2
    ┣ 1_initial_tree_construction.R
    ┣ 2_FS_taxa_graft.R
    ┣ 3_marsupial_graft.R
    ┣ 4_fossil_taxa_graft.R
    ┣ 5_make_additional_tip_and_node_priors.R
    ┣ 6_make_alternative_toplogies.R
    ┣ alternate_topologies.csv


The first four numbered scripts build an informal supertree by grafting together three dated phylogenies and several additional tips. We then carry out a tip-dating analysis using BEAST2. First we used BEAUTi to make `BEAST2/initial_beast_input_file.xml`, then the script `5_make_additional_tip_and_node_priors.R` produces several other xml files containing additional priors for the BEAST analysis and writes them into the 📁 BEAST2 folder. We then manually edited BEAST xml file, adding those additional priors, to make `BEAST2/updated_beast_input_file.xml`. That xml is executed in BEAST2, producing a posterior distribution of trees. We used TreeAnnotator to produce a single summary tree (maximum clade credibility with mean node heights)