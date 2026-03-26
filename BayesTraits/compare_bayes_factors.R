source("./scripts/BayesTraitsModelComparison.R")
stone_paths <- list.files(path="/BT_results/whole_brain/hyp_1/pPCscores/", #put your file path here, our stones files were stored elsewhere on another drive rather than in the project folder
                          pattern="a_run.txt.Stones.txt",
                          recursive = T,
                          full.names = T)

stone_paths # confirm if correct


#Calculate marginal likelihoods

# make sure the order of the labels matches the order of the path names

marginal_likelihoods <- getStones(stone_paths, labels = c("Single\nRate\nBM",
                                                          "Single\nRate\ndelta",
                                                          "Single\nRate\nkappa",
                                                          "Single\nRate\nlambda",
                                                          "Single\nRate\nOU",
                                                          "Variable\nRate\nBM",
                                                          "Variable\nRate\ndelta",
                                                          "Variable\nRate\nkappa",
                                                          "Variable\nRate\nlambda",
                                                          "Variable\nRate\nOU"))

plot(marginal_likelihoods)
