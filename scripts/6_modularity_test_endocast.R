library(geomorph)

# load the data
# aligned specimens
n_landmarks <- 122
n_dimensions <- 3
raw_data <- read.csv("./input_data/procrutes_aligned_coords.csv", row.names = 1)
Y.gpa.endo <- arrayspecs(raw_data, p = n_landmarks, k = n_dimensions)



# import point by point module definition file
module.df <- as.data.frame(read.csv("./input_data/lmk_module_guide.csv"))



# designate module guide to use
module.guide <-as.vector(module.df$X6)

# run CR analysis for each module hypothesis

MT.2 <-modularity.test(Y.gpa.endo, module.df$X2)
MT.3 <-modularity.test(Y.gpa.endo, module.df$X3)
MT.4 <-modularity.test(Y.gpa.endo, module.df$X4)
MT.5 <-modularity.test(Y.gpa.endo, module.df$X5)
MT.6 <-modularity.test(Y.gpa.endo, module.df$X6)


# compare CR test
MT.compare <- compare.CR(MT.2,MT.3,MT.4,MT.5,MT.6,
                               CR.null = TRUE)
summary(MT.compare)


#################################

#### Phylogenetic modularity ####

#################################

Phy.MT.2 <-phylo.modularity(Y.gpa.endo, module.df$X2, phy = tree_hyp1)
Phy.MT.3 <-phylo.modularity(Y.gpa.endo, module.df$X3, phy = tree_hyp1)
Phy.MT.4 <-phylo.modularity(Y.gpa.endo, module.df$X4, phy = tree_hyp1)
Phy.MT.5 <-phylo.modularity(Y.gpa.endo, module.df$X5, phy = tree_hyp1)
Phy.MT.6 <-phylo.modularity(Y.gpa.endo, module.df$X6, phy = tree_hyp1)

Phy.MT.compare <- compare.CR(Phy.MT.2,Phy.MT.3,Phy.MT.4,Phy.MT.5,Phy.MT.6,
                         CR.null = TRUE)
summary(Phy.MT.compare)

