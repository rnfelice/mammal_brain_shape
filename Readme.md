# Data and Code from Knapp et al.
This repository contains all data and code for carrying out the analyses in _Mosaic evolution of brain shape and size in mammals_. 

## Overview:

📁 BayesTraits: Code for carrying out BayesTraits analyses and processing results

📁 endocat_pts: Raw 3D landmark configurations

📁 input_data: Morphometric, ecological, phylogenetic, and body size data for analyses

📁 mvMORPH: Code for evolutionary rate analysis with mvMORPH

📁 phylogeny_construction: code for building dated supertrees with BEAST 2

📁 ply_ASCII: 

📁 scripts: R scripts for comparative analyses including convergent evolution tests

## Detailed Folder Contents

📁 BayesTraits: this folder contains the script `BT_input_prep.R`. This script is used to prepare the inputs for BayesTraits evolutionary rate analyses has everal functions: (1) Carry out phylogenetic PCA to reduce the dimensionality of the Procrustes-aligned landmark data (2) create BayesTraits command file (3) organise these input files into folders to analyses the evolutiooanry rates of the whole endocasts and each subregion (cerebellum, neocortex, olfactory lobe and brainstem). 

Other contents of this folder:

📁 brainsize: which does the same task for preparing inputs for the variable-rates phylogenetic regression in BayesTriats to analyse brain size allometry evolution. 
📁 postprocessing: contains the processed output of the BayesTraits models with the highest marginal likelihood. The raw output files from BayesTraits were post-processed using the variable-rates postprocessor here https://www.evolution.reading.ac.uk/VarRatesWebPP/

Finally, several R scripts provide the analysis and plotting code for the results of the BayesTraits analyses:
`check_convergence.R` uses Gelman-Rubin test diagnostic to confirm good convergence of MCMC chains
`compare_alternative_topologies.R` Makes the supplemental tree and rate-through-time figures for all 8 alternative phylogenetic hypotheses 
`compare_bayes_factors.R` is used to compare the marginal likelihoods of BayesTraits results for model selection. Makes figures S19-S31
`plot_figure_2.R` Makes main text figure 2
`plot_rates_per_region.R` makes main text figure 3 and supplmental figure S16

---

📂input_data: raw data for the analyses

    📦input_data
    ┣ 📂plotting_links
    ┣ 📜Proc.endo.mammals.rda
    ┣ 📜centroid_sizes.csv
    ┣ 📜mammal_data.csv
    ┣ 📜mammaltree4705_editednames.nex
    ┣ 📜procrutes_aligned_coords.csv
    ┗ 📜raw_landmark_coords.csv

📂plotting_links contains files for making connections between landmark points in plots (e.g., Figs S12-S15). `Proc.endo.mammals.rda` is the raw output of Procrustes alignment containing coordinates, centroid size, mean shape, etc. `centroid_sizes.csv` contains the centroid size of landmark configurations. `mammal_data.csv` contains the taxonomic, ecological, life history, and body size data for all taxa. Importantly, this file also contains links to the original source files for the 3D scan and/or 3D model data for each specimen in the dataset. `mammaltree4705_editednames.nex` is the phylogeny from Alvarez Carretero et al with genus names modified to be capitalized. `procrutes_aligned_coords.csv` contains the Procrustes-aligned landmark configurations whereas `raw_landmark_coords.csv` is the landmark data with semilandmarks having been slid to minimize bending energy but not yet Procrustes aligned.

---

📁 mvMORPH: These scripts are used to compare rates of evolution across ecological/behavioral groups using the mvMORPH R pacakge and create figures 4 and S17-18. `mvMORPH_hyp1.R` analyses per-group rates with phylogenetic hypothesis 1. `mvMORPH_hyp4.R` carries out a sensitivity analysis with an alternative phylogenetic hypotheis with Pantodonta, Brontotheriidae, and Creodonta in alternative positions in the tree. `mvMORPH_no_cetaceans_anthropoids.R` carries out a sensitivity analysis with phylogenetic hypotheis 1 with but the potentially influential outlier groups Cetacea (whales) and Anthropods (humans and close relatives) removed. 


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


The first four numbered scripts build an informal supertree by grafting together three dated phylogenies and several additional tips. We then carry out a tip-dating analysis using BEAST2. First we used BEAUTi to make `BEAST2/initial_beast_input_file.xml`, then the script `5_make_additional_tip_and_node_priors.R` produces several other xml files containing additional priors for the BEAST analysis and writes them into the 📁 BEAST2 folder. We then manually edited BEAST xml file, adding those additional priors, to make `BEAST2/updated_beast_input_file.xml`. That xml is executed in BEAST2, producing a posterior distribution of trees. We used TreeAnnotator to produce a single summary tree (maximum clade credibility with mean node heights). Finally, `6_make_alternative_toplogies.R` uses the `move.lineage` function in the RRPhylo R package to create alternative phylogenetic hypotheses by moving clades with uncertain phylogenetic affinities to alternative parts of the tree. the folder 📁 input_trees contains the published time trees used to assemble our supertree topology. The results are the 8 numbered nexus files `final_tree_hyp1.nex` through `final_tree_hyp8.nex` which are used for all comparative analyses and evolutionary models. 

---

📂ply_ASCII: contains the 3D mesh file for *Abrocoma cinerea*. Used for plotting. This model corresponds to this Morphosource record: https://n2t.net/ark:/87602/m4/M141118 

---
📁 scripts: Contains R scripts for main analyses

    📦scripts
    ┣ 📂utilities
    ┃ ┣ 📜BayesTraitsModelComparison.R
    ┃ ┣ 📜BayesTraitsPlottingFunctions.R
    ┃ ┗ 📜Functions.R
    ┣ 📜1_Procrustes_alignment_and_PCA.R
    ┣ 📜2_allometry.R
    ┣ 📜3_PhySig_and_rates.R
    ┣ 📜4_convergent_evolution_script.R
    ┗ 📜5_brain_region_analysis.R

📂utilities contains utility functions for reading and manipulating landmark data, model comparison including calculating and plotting BayesFactors, and functions for reading BayesTraits results into R.`1_Procrustes_alignment_and_PCA.R` Carries out PCA and make morphospace plots. `2_allometry.R` fits phylogenetic regressions to test allometric effects. `3_PhySig_and_rates.R` calculates phylogenetic signal evolutionary rates (sigma mult). `4_convergent_evolution_script.R` tests convergent evolution (Fig S3). `5_brain_region_analysis.R` makes per-region morphospace plots (Fig S12-S15)

