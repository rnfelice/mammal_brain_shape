####BayesTraits Plotting functions
#Modified from https://github.com/IanGBrennan/Tiliquini/
# process one or more BayesTraits PPPostProcess file(s)
require(Rmisc)
require(phytools)
process_PPP <- function(res.files, phy, col.palette, save_summary_trees = F, relative.rates=F){
  # warn users (me) that they will now be using relative rates
  #warning("now converting to relative rates (via mean scalars) \n
  #  this will affect downstream visualizations")
  # make some empty objects for us to store results
  all.trees <- NULL; all.names <- NULL; all.sig <- NULL; 
  all.res <- NULL; scalar.trees <- NULL; rate.trees <- NULL
  m.scalar.trees <- NULL;
  for (p in 1:length(res.files)){
    # get the name of the trait
    #trait.name <- strsplit(res.files[[p]], "_")[[1]][1]
    # read in the file
    BTres <- read.delim(res.files[[p]])
    # loop through and establish the ape-style node numbers from BayesTraits
    node.no <- NULL
    for (k in 1:length(BTres$Taxa.List)){
      curr.taxa <- strsplit(BTres$Taxa.List[[k]], ",")[[1]]
      if(length(curr.taxa)>1){node.no <- append(node.no, getMRCA(phy=phy, tip=curr.taxa))}
      if(length(curr.taxa)==1){node.no <- append(node.no, which(phy$tip.label==curr.taxa))}
    }
    # add a column for the new node numbers
    BTres$Node.No <- node.no
    # reorder the dataframe by new node numbers
    BTres <- BTres[order(BTres$Node.No),]
    # remove the root edge (Node.ID==0)
    BTres <- BTres[-which(BTres$Node.ID==0),]
    # get the edge/branch number for each
    BTres$edge <- unlist(sapply(BTres$Node.No, function(x) which(phy$edge[,2]==x)))
    #edge.vec <- append(edge.vec, 0, after=Ntip(phy)) # have to do this because the root is NA in above
    #BTres$edge <- edge.vec
    
    # establish the parent/child nodes
    BTres$child.node <- BTres$Node.No
    BTres$parent.node <- unlist(sapply(BTres$child.node, function(x) getParent(phy, x)))
    BTres$Node.No[1:Ntip(phy)] <- phy$tip.label
    
    # establish the start and end times of each branch/edge
    BTres$timestart <- sapply(BTres$parent.node, function(x) nodeheight(phy, x))
    BTres$timestop <-  sapply(BTres$child.node,  function(x) nodeheight(phy, x))
    BTres$timestart <- BTres$timestart - max(nodeHeights(phy))
    BTres$timestop <-  BTres$timestop  - max(nodeHeights(phy))
    BTres$timestop <- round(BTres$timestop, 5)
    
    # establish the start and end rates of each branch/edge
    # get the trait values at the parent and child nodes (start and end of edge)
    root.rate <- function(y){
      if(y==Ntip(phy)+1){0}
      else{BTres$Mean.SigV[which(BTres$child.node==y)]}
    }
    BTres$ratestart <- unlist(sapply(BTres$parent.node, root.rate))
    BTres$ratestop  <- unlist(sapply(BTres$child.node,  function(x) BTres$Mean.SigV[which(BTres$child.node==x)]))
    
    # establish the start and end rates of each branch/edge
    # get the trait values at the parent and child nodes (start and end of edge)
    root.scalar <- function(y){
      if(y==Ntip(phy)+1){0}
      else{BTres$Mean.Scalar[which(BTres$child.node==y)]}
    }    
    BTres$scalarstart <- unlist(sapply(BTres$parent.node, root.scalar))
    BTres$scalarstop  <- unlist(sapply(BTres$child.node,  function(x) BTres$Mean.Scalar[which(BTres$child.node==x)]))
    BTres$edge.scalar <- BTres$Mean.Scalar
    
    # scale the rates by the fastest rate (for ease of plotting)
    BTres$edge.rate <- BTres$Mean.SigV
    BTres$rounded.rates <- round((BTres$Mean.SigV - min(BTres$Mean.SigV))/diff(range(BTres$Mean.SigV)) * 99) + 1
    
    # make a color ramp that suits the data
    if(col.palette %in% c("magma", "inferno", "plasma", "viridis","cividis", "rocket", "mako")){
      new.cols <- viridis(n=100, option=col.palette, direction = -1)
    }else{
      col.ramp <- colorRampPalette(brewer.pal(9, col.palette)[2:9])
      new.cols <- (col.ramp(100))
      if(col.palette %in% c("Spectral","RdYlGn","RdYlBu","RdGy","RdBu","PuOr","PRGn","PiYG","BrBG")){new.cols <- rev(new.cols)}
    }
    
    #all.edges$edge.color <- sapply(all.edges$rounded.rates, function(x) new.cols[x]) # this is definitely wrong
    BTres$edge.color <- new.cols[BTres$rounded.rates]
    BTres$palette <- col.palette
    
    # identify branches with "positive phenotypic selection" sensu Baker et al. 2016
    BTres$Mean.deltaVB <- BTres$Mean.DeltaV/BTres$Mean.DeltaB
    BTres$Mode.deltaVB <- BTres$Mode.DeltaV/BTres$Mode.DeltaB
    
    # save the results df
    all.res[[p]] <- dplyr::select(BTres, Original.BL, 
                                  Mean.Scalar, Median.Scalar,
                                  Mean.SigV, Median.SigV, 
                                  Mean.deltaVB, Mode.deltaVB, Pct.time.scaled,
                                  Node.No, edge,
                                  child.node, parent.node,
                                  timestart, timestop,
                                  ratestart, ratestop,
                                  scalarstart, scalarstop, 
                                  edge.rate, edge.scalar, rounded.rates, 
                                  edge.color, palette)
    # make a vector of sigma and scale the tree by the sigma rate
    sig.bt <- BTres$Mean.SigV; names(sig.bt) <- BTres$child.node
    all.sig[[p]] <- sig.bt
    bt.tree <- phy
    bt.tree$edge.length <- sapply(1:nrow(bt.tree$edge), function(x) bt.tree$edge.length[x] * sig.bt[which(names(sig.bt)==bt.tree$edge[x,2])])
    plot(bt.tree, cex=0.3); axisPhylo()
    if (save_summary_trees){
    write.tree(bt.tree, file=paste0("BayesTraits_VarRates_",".tre"))
    }
    all.trees[[p]] <- bt.tree
    
    # make a vector of sigma and change the tree edges to match the rates
    rt.tree <- phy
    rt.tree$edge.length <- sapply(1:nrow(rt.tree$edge), function(x) sig.bt[which(names(sig.bt)==rt.tree$edge[x,2])])
    rate.trees[[p]] <- rt.tree
    
    # make a vector of median scalar and scale the tree by the median scalar
    scalar.bt <- BTres$Median.Scalar; names(scalar.bt) <- BTres$child.node
    sc.tree <- phy
    sc.tree$edge.length <- sapply(1:nrow(sc.tree$edge), function(x) sc.tree$edge.length[x] * scalar.bt[which(names(scalar.bt)==sc.tree$edge[x,2])])
    #plot(bt.tree, cex=0.3); axisPhylo()
    #write.tree(bt.tree, file=paste0("BayesTraits_VarRates_",trait.name,".tre"))
    scalar.trees[[p]] <- sc.tree
    
    # make a vector of mean scalar and scale the tree by the mean scalar
    m.scalar.bt <- BTres$Mean.Scalar; names(m.scalar.bt) <- BTres$child.node
    m.sc.tree <- phy
    m.sc.tree$edge.length <- sapply(1:nrow(m.sc.tree$edge), function(x) m.sc.tree$edge.length[x] * m.scalar.bt[which(names(m.scalar.bt)==m.sc.tree$edge[x,2])])
    m.scalar.trees[[p]] <- m.sc.tree
    
    # make a vector of 
    
  }
  names(all.trees) <- all.names; names(all.sig) <- all.names; 
  names(all.res) <- all.names; names(scalar.trees) <- all.names; 
  names(rate.trees) <- all.names; names(m.scalar.trees) <- all.names;
  class(all.trees) <- "multiPhylo"; class(scalar.trees) <- "multiPhylo"; class(m.scalar.trees) <- "multiPhylo"
  write.nexus(all.trees, file="BayesTraits_VarRates_SigV.trees")
  write.nexus(scalar.trees, file="BayesTraits_VarRates_MedianScalar.trees")
  write.nexus(m.scalar.trees, file="BayesTraits_VarRates_MeanScalar.trees")
  
  return(list(all.trees = all.trees, all.sig = all.sig, all.res = all.res, 
              scalar.trees = scalar.trees, mean.scalar.trees = m.scalar.trees, rate.trees = rate.trees))
}
# e.g.: testo <- process_PPP(res.files = in.files, phy=egernia.tree)

process_FPP <- function(res.files, phy, trait.name=NULL,save_summary_trees=F){
  # make some empty objects for us to store results
  all.trees <- NULL; all.names <- NULL; all.scalar <- NULL; 
  all.res <- NULL; scalar.trees <- NULL; rate.trees <- NULL
  m.scalar.trees <- NULL;
  for (p in 1:length(res.files)){
    # get the name of the trait
    if(is.null(trait.name)){
    trait.name <- strsplit(res.files[[p]], "_")[[1]][1]
    }else{
    trait.name<-trait.name}
    # read in the file
    BTres <- read.delim(res.files[[p]], sep=",")
    colnames(BTres)[c(16,17,18,29,30,31)] <- c("Beta_P_less_0", "Beta_P_equal_0", "Beta_P_greater_0",
                                               "Scalar_P_less_1", "Scalar_P_equal_1", "Scalar_P_greater_1")
    # loop through and establish the ape-style node numbers from BayesTraits
    node.no <- NULL
    for (k in 1:nrow(BTres)){
      curr.taxa <- unlist(BTres[k,44:ncol(BTres)]); names(curr.taxa) <- NULL
      curr.taxa <- curr.taxa[curr.taxa != ""]
      if(length(curr.taxa)>1){node.no <- append(node.no, getMRCA(phy=phy, tip=curr.taxa))}
      if(length(curr.taxa)==1){node.no <- append(node.no, which(phy$tip.label==curr.taxa))}
    }
    # add a column for the new node numbers
    BTres$Node.No <- node.no
    # reorder the dataframe by new node numbers
    BTres <- BTres[order(BTres$Node.No),]
    # remove the root edge (Node.ID==0)
    BTres <- BTres[-which(BTres$ID==0),]
    # get the edge/branch number for each
    BTres$edge <- unlist(sapply(BTres$Node.No, function(x) which(phy$edge[,2]==x)))
    #edge.vec <- append(edge.vec, 0, after=Ntip(phy)) # have to do this because the root is NA in above
    #BTres$edge <- edge.vec
    
    # establish the parent/child nodes
    BTres$child.node <- BTres$Node.No
    BTres$parent.node <- unlist(sapply(BTres$child.node, function(x) getParent(phy, x)))
    BTres$Node.No[1:Ntip(phy)] <- phy$tip.label
    
    # establish the start and end times of each branch/edge
    BTres$timestart <- sapply(BTres$parent.node, function(x) nodeheight(phy, x))
    BTres$timestop <-  sapply(BTres$child.node,  function(x) nodeheight(phy, x))
    BTres$timestart <- BTres$timestart - max(nodeHeights(phy))
    BTres$timestop <-  BTres$timestop  - max(nodeHeights(phy))
    BTres$timestop <- round(BTres$timestop, 5)
    
    # establish the start and end rates of each branch/edge
    # get the trait values at the parent and child nodes (start and end of edge)
    #root.rate <- function(y){
    #  if(y==Ntip(phy)+1){0}
    #  else{BTres$Mean.SigV[which(BTres$child.node==y)]}
    #}
    #BTres$ratestart <- unlist(sapply(BTres$parent.node, root.rate))
    #BTres$ratestop  <- unlist(sapply(BTres$child.node,  function(x) BTres$Mean.SigV[which(BTres$child.node==x)]))
    
    # scale the rates by the fastest rate (for ease of plotting)
    #BTres$edge.rate <- BTres$Mean.SigV
    #BTres$rounded.rates <- round((BTres$Mean.SigV - min(BTres$Mean.SigV))/diff(range(BTres$Mean.SigV)) * 99) + 1
    
    # make a color ramp that suits the data
    #if(col.palette %in% c("magma", "inferno", "plasma", "viridis","cividis", "rocket", "mako")){
    #  new.cols <- viridis(n=100, option=col.palette)
    #}else{
    #  col.ramp <- colorRampPalette(brewer.pal(9, col.palette))
    #  new.cols <- (col.ramp(100))
    #  if(col.palette %in% c("Spectral","RdYlGn","RdYlBu","RdGy","RdBu","PuOr","PRGn","PiYG","BrBG")){new.cols <- rev(new.cols)}
    #}
    
    #all.edges$edge.color <- sapply(all.edges$rounded.rates, function(x) new.cols[x]) # this is definitely wrong
    #BTres$edge.color <- new.cols[BTres$rounded.rates]
    #BTres$palette <- col.palette
    
    # identify branches with significant 
    
    # save the results df
    all.res[[p]] <- dplyr::select(BTres,
                                  ID, Branch.Length,
                                  Mean..Beta...BL..NZ,
                                  Beta_P_less_0, Beta_P_equal_0, Beta_P_greater_0,
                                  Mean.Non.1.Scalar, P..Beta...BL.,
                                  Scalar_P_less_1, Scalar_P_equal_1, Scalar_P_greater_1,
                                  "Z..Beta...BL.","Sig..Beta...BL.",
                                  P.Scalar,
                                  Prop.beta.followed.by.a.variance.scalar,
                                  Prop.Variance.scalar.followed.by.a.beta,
                                  Z.Scalar,Sig.Scalar,
                                  Node.No, edge,
                                  child.node, parent.node,
                                  timestart, timestop)
    
    
    # make a vector of scalar values and scale the tree
    #sc.bt <- BTres$Mean.Non.1.Scalar; names(sc.bt) <- BTres$child.node
    #all.scalar[[p]] <- sc.bt
    #sc.tree <- phy
    #bt.tree$edge.length <- sapply(1:nrow(bt.tree$edge), function(x) bt.tree$edge.length[x] * sig.bt[which(names(sig.bt)==bt.tree$edge[x,2])])
    #plot(bt.tree, cex=0.3); axisPhylo()
    #write.tree(bt.tree, file=paste0("BayesTraits_VarRates_",trait.name,".tre"))
    all.trees[[p]] <- phy
    all.names[[p]] <- trait.name
    
    # make a vector of scalar scalar and scale the tree
    scalar.bt <- BTres$Mean.Non.1.Scalar; names(scalar.bt) <- BTres$child.node
    sc.tree <- phy
    sc.tree$edge.length <- sapply(1:nrow(sc.tree$edge), function(x) sc.tree$edge.length[x] * scalar.bt[which(names(scalar.bt)==sc.tree$edge[x,2])])
    scalar.trees[[p]] <- sc.tree
    
    # make a vector of sigma and change the tree edges to match the rates
    #rt.tree <- phy
    #rt.tree$edge.length <- sapply(1:nrow(rt.tree$edge), function(x) sig.bt[which(names(sig.bt)==rt.tree$edge[x,2])])
    #rate.trees[[p]] <- rt.tree
    
    
    
    # make a vector of mean scalar and scale the tree by the mean scalar
    m.scalar.bt <- BTres$Mean.Scalar; names(m.scalar.bt) <- BTres$child.node
    m.sc.tree <- phy
    m.sc.tree$edge.length <- sapply(1:nrow(m.sc.tree$edge), function(x) m.sc.tree$edge.length[x] * m.scalar.bt[which(names(m.scalar.bt)==m.sc.tree$edge[x,2])])
    m.scalar.trees[[p]] <- m.sc.tree
    
    # make a vector of 
    
  }

  names(all.res) <- all.names; names(scalar.trees) <- all.names; 
  names(m.scalar.trees) <- all.names;
  class(all.trees) <- "multiPhylo"; class(scalar.trees) <- "multiPhylo"; class(m.scalar.trees) <- "multiPhylo"
  if(save_summary_trees){
  write.nexus(all.trees, file="BayesTraits_VarRates_SigV.trees")
  write.nexus(scalar.trees, file="BayesTraits_VarRates_MedianScalar.trees")
  write.nexus(m.scalar.trees, file="BayesTraits_VarRates_MeanScalar.trees")
  }
  return(list(all.trees = all.trees, all.res = all.res, 
              scalar.trees = scalar.trees, mean.scalar.trees = m.scalar.trees, rate.trees = rate.trees))
}
# e.g.: testo <- process_PPP(res.files = in.files, phy=egernia.tree)


process_combined_FPP <- function(res.file, phy, trait.name=NULL,save_summary_trees=F){
  # make some empty objects for us to store results
  all.trees <- NULL; all.names <- NULL; all.scalar <- NULL; 
  all.res <- NULL; scalar.trees <- NULL; rate.trees <- NULL
  m.scalar.trees <- NULL;

      trait.name<-trait.name
    # read in the file
    BTres <- read.delim(res.file, sep=",")
    lookup <- c("P...0"  = "Beta_P_less_0",
                                               "P....0" = "Beta_P_equal_0",
                                               "P...0.1" = "Beta_P_greater_0",
                                               "P...1" = "Scalar_P_less_1",
                                               "P....1" = "Scalar_P_equal_1",
                                               "P...1.1"  ="Scalar_P_greater_1")
    
    for (i in names(lookup)) {
      colnames(BTres)[colnames(BTres) == i] <- lookup[[i]]
    }
    # loop through and establish the ape-style node numbers from BayesTraits
    node.no <- NULL
    first_taxon_column<-grep("^Taxa",colnames(BTres))[1]
    for (k in 1:nrow(BTres)){
      curr.taxa <- unlist(BTres[k,first_taxon_column:ncol(BTres)]); names(curr.taxa) <- NULL
      curr.taxa <- curr.taxa[curr.taxa != ""]
      if(length(curr.taxa)>1){node.no <- append(node.no, getMRCA(phy=phy, tip=curr.taxa))}
      if(length(curr.taxa)==1){node.no <- append(node.no, which(phy$tip.label==curr.taxa))}
    }
    # add a column for the new node numbers
    BTres$Node.No <- node.no
    # reorder the dataframe by new node numbers
    BTres <- BTres[order(BTres$Node.No),]
    # remove the root edge (Node.ID==0)
    BTres <- BTres[-which(BTres$ID==0),]
    # get the edge/branch number for each
    BTres$edge <- unlist(sapply(BTres$Node.No, function(x) which(phy$edge[,2]==x)))
    #edge.vec <- append(edge.vec, 0, after=Ntip(phy)) # have to do this because the root is NA in above
    #BTres$edge <- edge.vec
    
    # establish the parent/child nodes
    BTres$child.node <- BTres$Node.No
    BTres$parent.node <- unlist(sapply(BTres$child.node, function(x) getParent(phy, x)))
    BTres$Node.No[1:Ntip(phy)] <- phy$tip.label
    
    # establish the start and end times of each branch/edge
    BTres$timestart <- sapply(BTres$parent.node, function(x) nodeheight(phy, x))
    BTres$timestop <-  sapply(BTres$child.node,  function(x) nodeheight(phy, x))
    BTres$timestart <- BTres$timestart - max(nodeHeights(phy))
    BTres$timestop <-  BTres$timestop  - max(nodeHeights(phy))
    BTres$timestop <- round(BTres$timestop, 5)
    
    # establish the start and end rates of each branch/edge
    # get the trait values at the parent and child nodes (start and end of edge)
    #root.rate <- function(y){
    #  if(y==Ntip(phy)+1){0}
    #  else{BTres$Mean.SigV[which(BTres$child.node==y)]}
    #}
    #BTres$ratestart <- unlist(sapply(BTres$parent.node, root.rate))
    #BTres$ratestop  <- unlist(sapply(BTres$child.node,  function(x) BTres$Mean.SigV[which(BTres$child.node==x)]))
    
    # scale the rates by the fastest rate (for ease of plotting)
    #BTres$edge.rate <- BTres$Mean.SigV
    #BTres$rounded.rates <- round((BTres$Mean.SigV - min(BTres$Mean.SigV))/diff(range(BTres$Mean.SigV)) * 99) + 1
    
    # make a color ramp that suits the data
    #if(col.palette %in% c("magma", "inferno", "plasma", "viridis","cividis", "rocket", "mako")){
    #  new.cols <- viridis(n=100, option=col.palette)
    #}else{
    #  col.ramp <- colorRampPalette(brewer.pal(9, col.palette))
    #  new.cols <- (col.ramp(100))
    #  if(col.palette %in% c("Spectral","RdYlGn","RdYlBu","RdGy","RdBu","PuOr","PRGn","PiYG","BrBG")){new.cols <- rev(new.cols)}
    #}
    
    #all.edges$edge.color <- sapply(all.edges$rounded.rates, function(x) new.cols[x]) # this is definitely wrong
    #BTres$edge.color <- new.cols[BTres$rounded.rates]
    #BTres$palette <- col.palette
    
    # identify branches with significant 
    
    # save the results df
    all.res <- dplyr::select(BTres,
                                  ID, Branch.Length,
                                  Mean..Beta...BL..NZ,
                                  Beta_P_less_0, Beta_P_equal_0, Beta_P_greater_0,
                                  Mean.Non.1.Scalar, P..Beta...BL.,
                                  Scalar_P_less_1, Scalar_P_equal_1, Scalar_P_greater_1,
                                  "Z..Beta...BL.","Sig..Beta...BL.",
                                  P.Scalar,
                                  Prop.beta.followed.by.a.variance.scalar,
                                  Prop.Variance.scalar.followed.by.a.beta,
                                  Z.Scalar,Sig.Scalar,
                                  Node.No, edge,
                                  child.node, parent.node,
                                  timestart, timestop)
    
    
    # make a vector of scalar values and scale the tree
    #sc.bt <- BTres$Mean.Non.1.Scalar; names(sc.bt) <- BTres$child.node
    #all.scalar[[p]] <- sc.bt
    #sc.tree <- phy
    #bt.tree$edge.length <- sapply(1:nrow(bt.tree$edge), function(x) bt.tree$edge.length[x] * sig.bt[which(names(sig.bt)==bt.tree$edge[x,2])])
    #plot(bt.tree, cex=0.3); axisPhylo()
    #write.tree(bt.tree, file=paste0("BayesTraits_VarRates_",trait.name,".tre"))
    all.trees <- phy
    all.names <- trait.name
    
    # make a vector of scalar scalar and scale the tree
    scalar.bt <- BTres$Mean.Non.1.Scalar; names(scalar.bt) <- BTres$child.node
    sc.tree <- phy
    sc.tree$edge.length <- sapply(1:nrow(sc.tree$edge), function(x) sc.tree$edge.length[x] * scalar.bt[which(names(scalar.bt)==sc.tree$edge[x,2])])
    scalar.trees <- sc.tree
    
    # make a vector of sigma and change the tree edges to match the rates
    #rt.tree <- phy
    #rt.tree$edge.length <- sapply(1:nrow(rt.tree$edge), function(x) sig.bt[which(names(sig.bt)==rt.tree$edge[x,2])])
    #rate.trees[[p]] <- rt.tree
    
    
    
    # make a vector of mean scalar and scale the tree by the mean scalar
    m.scalar.bt <- BTres$Mean.Non.1.Scalar; names(m.scalar.bt) <- BTres$child.node
    m.sc.tree <- phy
    m.sc.tree$edge.length <- sapply(1:nrow(m.sc.tree$edge), function(x) m.sc.tree$edge.length[x] * m.scalar.bt[which(names(m.scalar.bt)==m.sc.tree$edge[x,2])])
    m.scalar.trees <- m.sc.tree
    
    # make a vector of 
    
  

  #class(all.trees) <- "multiPhylo"; class(scalar.trees) <- "multiPhylo"; class(m.scalar.trees) <- "multiPhylo"
  if(save_summary_trees){
    write.nexus(all.trees, file="BayesTraits_VarRates_SigV.trees")
    write.nexus(scalar.trees, file="BayesTraits_VarRates_MedianScalar.trees")
    write.nexus(m.scalar.trees, file="BayesTraits_VarRates_MeanScalar.trees")
  }
  return(list(all.trees = all.trees, all.res = all.res, 
              scalar.trees = scalar.trees, mean.scalar.trees = m.scalar.trees, rate.trees = rate.trees))
}

# plot a phylogeny with branches colored by rate
# coloring options are a bit tricky, so stick to 'log.rates=F relative.rates=T'
plot.VarRates.tree <- function(BT, phy, col.palette = "Blues", legend = F,
                               tree.type = c("phylogram", "fan"),
                               log.rates = F, relative.rates = F,
                               trait=NULL, outline=F,
                               pos.selection = F, shift.rates = F,
                               annotations=F, tip.labs=T,plot.axis=T){
  
  #ntraits <- length(RR)
  #opt.layout <- n2mfrow(ntraits)
  #if(!(ntraits %% 2) == 0){mat.len <- ntraits + 1} else{mat.len <- ntraits}
  
  #all.node.rates <- NULL
  
  # set the plot layout
  #layout(matrix(nrow = opt.layout[1], ncol = opt.layout[2], 1:mat.len))
  
  # Set the layout
  if(legend==T){layout(
    matrix(c(1,2,1,3), ncol=2, byrow=T), 
    widths=c(4,1), 
    heights=c(1,1))
  }
  
  # extract the rates and edge/node information
  #if(log.rates==T){BT$Mean.SigV <- log(BT$Mean.SigV)}
  all.edges <- data.frame(edge.rate = BT$Mean.SigV,
                          edge = BT$edge,
                          child.node = BT$child.node,
                          scalar = BT$Mean.Scalar) # could also use 'scalar = BT$edge.scalar'
  if(log.rates==T && relative.rates==F){all.edges$raw.rate <- all.edges$edge.rate; all.edges$edge.rate <- log(all.edges$edge.rate)}
  #if(log.rates==T && relative.rates==T){all.edges$edge.rate <- mean(all.edges$edge.rate) / all.edges$edge.rate}
  if(log.rates==T && relative.rates==T){all.edges$raw.rate <- all.edges$edge.rate / mean(all.edges$edge.rate); all.edges$edge.rate <- log(all.edges$edge.rate)}
  #if(log.rates==F && relative.rates==T){all.edges$edge.rate <- all.edges$edge.rate / mean(all.edges$edge.rate)} # this was the old dumb way of getting the relative rates!
  if(log.rates==F && relative.rates==T){all.edges$edge.rate <- all.edges$scalar}  # this is the correct way of getting relative rates because it's already available as the scalar!
  #all.edges <- all.edges[-which(all.edges$child.node==Ntip(phy)+1),]
  
  # scale the rates by the fastest rate (for ease of plotting)
  all.edges$rounded.rates <- round((all.edges$edge.rate - min(all.edges$edge.rate))/diff(range(all.edges$edge.rate)) * 99) + 1
  
  # sort it by the original order!
  all.edges <- all.edges[order(match(all.edges$child.node, phy$edge[,2])),]  
  
  # make a color ramp that suits the data
  if(col.palette %in% c("magma", "inferno", "plasma", "viridis","cividis", "rocket", "mako","turbo")){
    new.cols <- viridis(n=100, option=col.palette, direction = 1)
    #if(col.palette %in% c("mako", "plasma", "inferno")){new.cols <- rev(new.cols)}
  }else if(col.palette == "GyYlRd"){
    colfunc <- colorRampPalette(c("#c7c7c7", "#ffee2e", "#e81d13"))
    new.cols <- colfunc(100)
  }else if(col.palette == "GnBu"){
    colfunc <- colorRampPalette(c("#C7E9B4","#7FCDBB","#41B6C4","#2C7FB8","#253494","#081D58"))
    new.cols <- colfunc(100)
  }else{
    col.ramp <- colorRampPalette(brewer.pal(9, col.palette)[2:9])
    new.cols <- (col.ramp(100))
    if(col.palette %in% c("Spectral","RdYlGn","RdYlBu","RdGy","RdBu","PuOr","PRGn","PiYG","BrBG")){new.cols <- rev(new.cols)}
  }
  
  
  #all.edges$edge.color <- sapply(all.edges$rounded.rates, function(x) new.cols[x]) # this is definitely wrong
  all.edges$edge.color <- new.cols[all.edges$rounded.rates]
  
  # plot the phylogeny with colored branches
  if(outline==T){plot.phylo(phy, edge.width=5, type=tree.type, cex=0.3, show.tip.label = tip.labs); par(new=T)}
  plot.phylo(phy, edge.color = unlist(all.edges$edge.color), edge.width=3, type=tree.type, cex=0.3, open.angle=5, show.tip.label = tip.labs)
  if(plot.axis){axisPhylo()}
  title(trait)
  
  # plot the location of edges undergoing putative shifts in evolutionary rate
  if(shift.rates==T){
    sr <- dplyr::filter(BT, Pct.time.scaled >= 50 & Mean.Scalar >= 2)
    if(nrow(sr) > 0){
      for (j in 1:nrow(sr)){
        edgelabels(text="", edge=noquote(sr$edge[j]), frame="circle", bg="white", cex=0.4, col=0.5)
      }      
    }
    if(annotations==T){title(xlab = paste(nrow(sr), "instances of rate shifts"))}
  }
  
  # plot the location of edges undergoing shifts in evolutionary rate
  if(shift.rates==T){
    sr <- dplyr::filter(BT, Pct.time.scaled >= 70 & Pct.time.scaled < 95 & Mean.Scalar >= 2)
    if(nrow(sr) > 0){
      for (j in 1:nrow(sr)){
        edgelabels(text="", edge=noquote(sr$edge[j]), frame="circle", bg="grey", cex=0.3, col=0.5)
      }      
    }
    #if(annotations==T){title(xlab = paste(nrow(sr), "instances of rate shifts"))}
  }
  
  # plot the location of edges undergoing positive selection sensu Baker et al. 2016
  if(pos.selection==T){
    #g2 <- dplyr::filter(BT, Mean.deltaVB >= 2 & Pct.time.scaled >= 95)
    g2 <- dplyr::filter(BT, Mean.Scalar >= 2 & Pct.time.scaled >= 95)
    
    if(nrow(g2) > 0){
      for (j in 1:nrow(g2)){
        edgelabels(text="", edge=noquote(g2$edge[j]), frame="circle", bg="black", cex=0.2, col=0.5)
      }   
    }
    if(annotations==T){title(sub = paste(nrow(g2), "instances of positive selection"))}
  }
  
  # if we chose to plot a legend, do that now
  if(legend==T && log.rates==F && relative.rates==F){plot(0,type='n',axes=FALSE,ann=FALSE);
    color.bar(new.cols, 
              min=round(min(all.edges$edge.rate),5),
              max=round(max(all.edges$edge.rate),5))}
  if(legend==T && log.rates==F && relative.rates==T){plot(0,type='n',axes=FALSE,ann=FALSE);
    color.bar(new.cols, 
              min=round(min(all.edges$edge.rate),5),
              max=round(max(all.edges$edge.rate),5))}
  if(legend==T && log.rates==T && relative.rates==F){plot(0,type='n',axes=FALSE,ann=FALSE);
    color.bar(new.cols, 
              min=round(min(all.edges$raw.rate),5),
              max=round(max(all.edges$raw.rate),5))}
  if(legend==T && log.rates==T && relative.rates==T){plot(0,type='n',axes=FALSE,ann=FALSE);
    color.bar(new.cols, 
              min=round(min(log(all.edges$raw.rate)),2),
              max=round(max(log(all.edges$raw.rate)),2))}
}


# this function is just to plot the scale bar, which is dumb, but it works, so hey.
color.bar <- function(lut, min=0, max=100, nticks=11, ticks=seq(min, max, len=nticks), title='') {
  scale = (length(lut)-1)/(max-min)
  
  #dev.new(width=1.75, height=5)
  plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
  axis(2, ticks, las=1)
  for (i in 1:(length(lut)-1)) {
    y = (i-1)/scale + min
    rect(0,y,10,y+1/scale, col=lut[i], border=NA)
  }
}



# get the rates at timeslices from a processed dataframe
rate.at.time.df <- function(timeslices, obj, plot=F, relative.rates=c("F","mean","median")){
  # define the timeslices for extracting trait values, and apply it to your tree
  ts <- c(seq(from=min(c(obj$timestart,obj$timestop)), to=max(c(obj$timestart,obj$timestop)), by=timeslices), max(c(obj$timestart,obj$timestop)))
  
  # for each timepoint reconstruct the trait value as a function of the distance between the parent and child node
  Yt <- list()
  for (i in 1:(length(ts)-1)) {      
    # which edges overlap timeslice i
    curr.edges <- filter(obj, timestart <= ts[[i]] & timestop >= ts[[i]])
    
    # extract the rates for the branches of interest
    curr.rates <- curr.edges$edge.rate
    
    # apply the current slice time
    slice <- rep(ts[i],times = length(curr.rates))
    
    # add all the trait values at timeslice i to the list
    Yt[[i]] <- cbind(curr.rates, slice)
  }
  # remove the fractional period just before the tips (combine the penultimate and ultimate windows)
  ts2 <- ts[-(length(ts)-1)]
  
  # now make a data frame of the trait values at given times
  shuffled.ages <- NULL
  for (k in 1:length(ts2)){
    curr.slice <- unlist(Yt[[k]])
    curr.slice.df <- data.frame(curr.slice, "time" = ts2[k])
    shuffled.ages <- rbind(shuffled.ages, curr.slice.df)
  }
  
  # extract just the trait and timeslices
  rate.time <- shuffled.ages[,c("curr.rates","time")]
  
  # reorder the columns so 'time' comes first
  rate.time <- relocate(rate.time, time)
  colnames(rate.time) <- c("time", "rate")
  
  # rescale rates relative to the mean if requested
  if(relative.rates=="F"){rate.time$rate <- rate.time$rate}
  if(relative.rates=="mean"){rate.time$rate <- rate.time$rate/mean(rate.time$rate)}
  if(relative.rates=="median"){rate.time$rate <- rate.time$rate/median(rate.time$rate)}
  
  # plot the values if you're interested, it should look like the horizontal branches of the tree
  if(plot==T){plot(rate.time$rate ~ rate.time$time, xlab="Time", ylab="Rate", pch=16)}
  
  # give up the object
  return(rate.time)
}


extract.stat <- function(rate.time.obj, stat=c("mean", "median", "scale"), range=c("confidence","quantile"),
                         plot=c("average", "corrected", "sideXside", FALSE)){
  
  # create a data frame of the individual time slices
  time = unique(rate.time.obj$time)
  time.mean <- data.frame(time)
  
  # establish the minimum and maximum values at each slice
  if(stat=="mean"){rate.mean <- aggregate(rate.time.obj[,2], list(rate.time.obj$time), mean)}
  if(stat=="median"){rate.mean <- aggregate(rate.time.obj[,2], list(rate.time.obj$time), median)}
  if(stat=="scale"){
    rate.mean <- aggregate(rate.time.obj[,2], list(rate.time.obj$time), mean)
    rate.mean[,2] <- ((rate.mean[,2] - min(rate.mean[,2]))/diff(range(rate.mean[,2])) * 99) + 1
  }
  rate.sd <- aggregate(rate.time.obj[,2], list(rate.time.obj$time), sd)
  
  
  # if you want to estimate quantiles on the rates:
  if(range=="quantile"){
    rate.qt5 <- aggregate(rate.time.obj[,2], list(rate.time.obj$time), FUN='quantile', probs=0.025)
    rate.qt95 <- aggregate(rate.time.obj[,2], list(rate.time.obj$time), FUN='quantile', probs=0.975)
    rate.mean <- cbind(rate.mean, rate.sd[,2], rate.qt5[,2], rate.qt95[,2])    
  }
  
  # if you prefer to estimate confidence intervals on the rates:
  if(range=="confidence"){
    rate.all <- aggregate(rate.time.obj[,2], list(rate.time.obj$time), FUN = function(x) Rmisc::CI(x, ci = 0.95))
    rate.mean <- data.frame(time = time.mean[,1], rate = rate.mean[,2], 
                            sd = rate.sd[,2], "5%" = rate.all$x[,3], "95%" = rate.all$x[,1])   
  }
  # regardless, rename the columns
  colnames(rate.mean) <- c("time", "rate", "sd", "5%", "95%")
  
  # establish the number of species living in each time slice
  spp.slice <- aggregate(rate.time.obj[,2], list(rate.time.obj$time), length)
  
  # add the richness information to the data frame
  rate.mean$richness <- spp.slice[,2]
  
  # correct rate by the number of species living at that time
  rate.mean$rate.rich <- rate.mean$rate/rate.mean$richness
  
  # get the name of the current variable
  var.name <- colnames(rate.time.obj)[2]
  
  # plot the results if you'd like
  if(plot=="average"){plot(rate.mean$rate ~ rate.mean$time, 
                           xlab="Time", ylab=paste(stat, "Evolutionary Rate"), type="l")}
  if(plot=="corrected"){plot(rate.mean$rate.rich ~ rate.mean$time, 
                             xlab="Time", ylab=paste(stat, "Evolutionary Rate"), type="l")}
  if(plot=="sideXside"){layout(matrix(1:2,ncol=2))
    plot(rate.mean$rate ~ rate.mean$time, 
         xlab="Time", ylab=paste(stat, "Evolutionary Rate"), type="l")
    plot(rate.mean$rate.rich ~ rate.mean$time, 
         xlab="Time", ylab=paste(stat, "Evolutionary Rate"), type="l")
    layout(matrix(1))}
  if(plot==F){}
  
  return(rate.mean)
  
}

# example:
#extract.mean(rate.time.obj = testo[,c(1,5)], stat = "mean", plot = "sideXside")

# reorder any dataframe or matrix following the phylogenetic order
order.by.phylo <- function(jumbled.order, phy){
  # make sure the matrix order matches the tip + node order
  
  # if the data is in a vector (single trait)
  if(is.vector(jumbled.order)){
    phylo.order <- jumbled.order[order(match(names(jumbled.order), 
                                             c(phy$tip.label, (length(phy$tip.label)+1):(Ntip(phy)+Nnode(phy)))))]
  }
  # if the data is in a matrix or data frame (multiple traits)
  if(!is.vector(jumbled.order)){
    phylo.order <- jumbled.order[order(match(rownames(jumbled.order), 
                                             c(phy$tip.label, (length(phy$tip.label)+1):(Ntip(phy)+Nnode(phy))))),]
  }
  return(phylo.order)
}

# example:
# order.by.phylo(rate.matrix, egernia.tree)
