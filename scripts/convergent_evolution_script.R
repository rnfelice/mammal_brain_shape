

# 3. Convergence testing

library(convevol)

# 3.1 Designate groups for convergence testing
conv.dims <-PCA.3D.endo$x[,c(1:3)]

burrowers <- c("Notoryctes_typhlops","Amblysomus_hottentotus","Talpa_europaea",
               "Chlamyphorus_truncatus")

aquatic <-c("Enhydra_lutris","Neomonachus_tropicalis","Inia_geoffrensis")


# run Ct convergence check, using cumulative 95% of PCA values
conv.burrow <-convSigCt(tree_hyp1, PCA.3D.endo$x[,c(1:which(cumsum(phypc$d/sum(phypc$d))>0.95)[1])], 
                        focaltaxa = burrowers, groups = NULL, nsim = 1)


conv.burrow <- convratsig(tree_hyp1, PCA.3D.endo$x[,c(1:10)], 
                          convtips = burrowers, nsim = 100)

conv.burrow <-convSigCt(final_tree, PCA.3D.endo$x[,c(1:3)], 
                        focaltaxa = burrowers, nsim = 100)

plotCt(output = conv.burrow, phy = tree_hyp1, focaltaxa = burrowers)

conv.burrow.2 <-convSigCt(final_tree, PCA.3D.endo$x[,c(1:10)], 
                        focaltaxa = c("Notoryctes_typhlops","Amblysomus_hottentotus"), nsim = 100)

plotCt(output = conv.burrow.2, phy = final_tree, focaltaxa = c("Notoryctes_typhlops","Amblysomus_hottentotus"))



conv.aquatic <- convSigCt(tree_hyp1, PCA.3D.endo$x[,c(1:3)], 
                          focaltaxa = aquatic, nsim = 100)

plotCt(output = conv.aquatic, phy = tree_hyp1, focaltaxa = aquatic)



# 3.3 Convergence values for each morphotype
cluster <-gdf.eco.curves$cluster
# convtips <-tree.names
# convtips <-convtips[match(tree.names,names(cluster))]
# levels(gdf.eco.curves$cluster)

# convergence values
ljaw.conv.morph1 <-convrat(mandPhy, PCA.3D.nopatch$x[,c(1:35)], convtips = dimnames(Y.gpa.nopatch[,,morph1])[[3]])
ljaw.conv.morph2 <-convrat(mandPhy, PCA.3D.nopatch$x[,c(1:35)], convtips = dimnames(Y.gpa.nopatch[,,morph2])[[3]])
ljaw.conv.morph3 <-convrat(mandPhy, PCA.3D.nopatch$x[,c(1:35)], convtips = dimnames(Y.gpa.nopatch[,,morph3])[[3]])
ljaw.conv.morph4 <-convrat(mandPhy, PCA.3D.nopatch$x[,c(1:35)], convtips = dimnames(Y.gpa.nopatch[,,morph4])[[3]])
ljaw.conv.morph5 <-convrat(mandPhy, PCA.3D.nopatch$x[,c(1:35)], convtips = dimnames(Y.gpa.nopatch[,,morph5])[[3]])
ljaw.conv.morph6 <-convrat(mandPhy, PCA.3D.nopatch$x[,c(1:35)], convtips = dimnames(Y.gpa.nopatch[,,morph6])[[3]])



ljaw.conv.morph1 <- convSig(phy = mandPhy, traits = PCA.3D.nopatch, focaltaxa = dimnames(PCA.3D.nopatch[morph1,])[[1]], nsim = 100)





# Burrowers

burrowers <- c("Notoryctes_typhlops","Amblysomus_hottentotus","Talpa_europaea",
               "Chlamyphorus_truncatus")

# 1.1.3 Import mesh file for template specimen
Notoryctes.mesh <- vcgImport( paste( "./ply_ASCII/Notoryctes_typhlops_endocast.ply",
                                   sep = "" ) )
# plot mesh
shade3d( Notoryctes.mesh, col = "cornflowerblue", specular = "black" )


# 1.1.3 Import mesh file for template specimen
Amblysomus.mesh <- vcgImport( paste( "./ply_ASCII/Amblysomus_hottentotus_endocast.ply",
                                     sep = "" ) )
# plot mesh
shade3d( Amblysomus.mesh, col = "cornflowerblue", specular = "black" )


# 1.1.3 Import mesh file for template specimen
Chlamyphorus.mesh <- vcgImport( paste( "./ply_ASCII/Chlamyphorus_truncatus_endocast.ply",
                                     sep = "" ) )
# plot mesh
shade3d( Chlamyphorus.mesh, col = "cornflowerblue", specular = "black" )


# 1.1.3 Import mesh file for template specimen
Talpa.mesh <- vcgImport( paste( "./ply_ASCII/Talpa_europaea_endocast.ply",
                                       sep = "" ) )
# plot mesh
shade3d( Talpa.mesh, col = "cornflowerblue", specular = "black" )


# Swimmers


aquatic <-c("Enhydra_lutris","Neomonachus_tropicalis","Inia_geoffrensis")

# 1.1.3 Import mesh file for template specimen
Enhydra.mesh <- vcgImport( paste( "./ply_ASCII/Enhydra_lutris_endocast.ply",
                                     sep = "" ) )
# plot mesh
shade3d( Enhydra.mesh, col = "cornflowerblue", specular = "black" )


# 1.1.3 Import mesh file for template specimen
Neomonachus.mesh <- vcgImport( paste( "./ply_ASCII/Neomonachus_tropicalis_endocast.ply",
                                     sep = "" ) )
# plot mesh
shade3d( Neomonachus.mesh, col = "cornflowerblue", specular = "black" )


# 1.1.3 Import mesh file for template specimen
Inia.mesh <- vcgImport( paste( "./ply_ASCII/Inia_geoffrensis_endocast.ply",
                                       sep = "" ) )
# plot mesh
shade3d( Inia.mesh, col = "cornflowerblue", specular = "black" )



conv <-c(burrowers,aquatic)

col.conv.1 <-c("#f2464e", "#e1d69c","#1e798a","#efc35a",
               "#1e798a", "#1e798a","#1e798a")

names( col.conv.1 )   <- conv 

col.conv   <- col.conv.1[ match( conv, names( col.conv.1) ) ]


# plot PCA with convergent taxa
PhM.conv.1 <-ggplot()+ #PCs 1 and 2
  geom_segment(data=edgecoords.1, aes(x=x.x,xend=x.y, y=y.x, yend=y.y),col = "grey", size=0.5) +
  # geom_convexhull(data = PCs, aes(x=Comp1, y = Comp2, colour = superorder, fill = superorder), alpha = 0.2)+#,show.legend = FALSE)+
  geom_point(data = PCs, aes(x=Comp1, y = Comp2, fill = col.conv), fill = "grey", colour = "grey",shape = 21,  size=3)+#,show.legend = FALSE)+
  labs(x = xlab, y= ylab)+
  # scale_fill_manual(values = col.gp)+
  # scale_colour_manual(values = col.gp)+
  # geom_text(aes(label = shape.names), show.legend = FALSE)+
  theme_bw()
plot(PhM.conv.1)
