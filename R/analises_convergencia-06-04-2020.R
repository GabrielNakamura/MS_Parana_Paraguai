####loading packages####
library(vegan)
library(ape)
library(betapart)
library(SYNCSA)
library(ade4)
library(geiger)
library(phytools)

#######Taxonomic beta diversity#######
#between basins
beta.taxonomic<-beta.pair(comm.presau,index.family="sorensen") #taxonomic beta diversity
total.taxonomic<- beta.taxonomic$beta.sor
turn.taxonomic<-beta.taxonomic$beta.sim #turnover component of taxonomic beta diversity
nest.taxonomic<-beta.taxonomic$beta.sne #nestedness component of taxonomic beta diversity

#plotting total beta-div
pcoa.total.taxonomic<- cmdscale(sqrt(total.taxonomic), eig = T) #PCoA in taxonomic turnover matrix
scores.total.taxonomic<- vegan::scores(pcoa.total.taxonomic)
quartz()
plot(scores.total.taxonomic[,1],scores.total.taxonomic[,2],type="n",
     main="PCoA taxonomic beta-div")
points(scores.total.taxonomic[1:10,1],scores.total.taxonomic[1:10,2],pch=19,cex=1.8)
points(scores.total.taxonomic[11:20,1],scores.total.taxonomic[11:20,2],pch=2,cex=1.8)
text(scores.total.taxonomic[,1],scores.total.taxonomic[,2], labels=rownames(scores.total.taxonomic))
abline(v=0,h=0,lty=2)


#plotting turnover
pcoa.turn.taxonomic<- cmdscale(sqrt(turn.taxonomic), eig = T) #PCoA in taxonomic turnover matrix
scores.turn.taxonomic<- vegan::scores(pcoa.turn.taxonomic)
quartz()
plot(scores.turn.taxonomic[,1],scores.turn.taxonomic[,2],type="n",
     main="PCoA taxonomic turnover")
points(scores.turn.taxonomic[1:10,1],scores.turn.taxonomic[1:10,2],pch=19,cex=1.8)
points(scores.turn.taxonomic[11:20,1],scores.turn.taxonomic[11:20,2],pch=2,cex=1.8)
text(scores.turn.taxonomic[,1],scores.turn.taxonomic[,2], labels=rownames(scores.turn.taxonomic))
abline(v=0,h=0,lty=2)

#plotting nestedness
pcoa.nest.taxonomic<- cmdscale(sqrt(nest.taxonomic), eig = T)
scores.nest.taxonomic<- vegan::scores(pcoa.nest.taxonomic)
quartz()
plot(scores.nest.taxonomic[,1],scores.nest.taxonomic[,2],type="n",
     main="PCoA taxonomic nestedness")
points(scores.nest.taxonomic[1:10,1],scores.nest.taxonomic[1:10,2],pch=19,cex=1.8)
points(scores.nest.taxonomic[11:20,1],scores.nest.taxonomic[11:20,2],pch=2,cex=1.8)
text(scores.nest.taxonomic[,1],scores.nest.taxonomic[,2], labels=rownames(scores.nest.taxonomic))
abline(v=0,h=0,lty=2)

#######Phylogenetic beta diversity#######
beta.phylogenetic<-phylo.beta.pair(comm.presau,tree = phylo,index.family="sorensen") #phylogenetic beta diversity
total.phylogenetic<- beta.phylogenetic$phylo.beta.sor
turn.phylogenetic<-beta.phylogenetic$phylo.beta.sim #turnover component of phylogenetic beta diversity
nest.phylogenetic<-beta.phylogenetic$phylo.beta.sne #nestedness component of phylogenetic beta diversity

#phylogenetic total beta div
pcoa.total.phylogenetic<- cmdscale(sqrt(total.phylogenetic), eig = T) #PCoA in phylogenetic turnover matrix
scores.total.phylogenetic<- vegan::scores(pcoa.total.phylogenetic)
quartz()
plot(scores.total.phylogenetic[,1],scores.total.phylogenetic[,2],type="n",
     main="PCoA phylogenetic beta-div")
points(scores.total.phylogenetic[1:10,1],scores.total.phylogenetic[1:10,2],pch=19,cex=1.8)
points(scores.total.phylogenetic[11:20,1],scores.total.phylogenetic[11:20,2],pch=2,cex=1.8)
text(scores.total.phylogenetic[,1], (scores.total.phylogenetic[,2] + 0.02), labels=rownames(scores.total.phylogenetic))
abline(v=0,h=0,lty=2)

#phylogenetic turnover
pcoa.turn.phylogenetic<- cmdscale(sqrt(turn.phylogenetic), eig = T) #PCoA in phylogenetic turnover matrix
scores.turn.phylogenetic<- vegan::scores(pcoa.turn.phylogenetic)
quartz()
plot(scores.turn.phylogenetic[,1],scores.turn.phylogenetic[,2],type="n",
     main="PCoA phylogenetic turnover")
points(scores.turn.phylogenetic[1:10,1],scores.turn.phylogenetic[1:10,2],pch=19,cex=1.8)
points(scores.turn.phylogenetic[11:20,1],scores.turn.phylogenetic[11:20,2],pch=2,cex=1.8)
text(scores.turn.phylogenetic[,1], (scores.turn.phylogenetic[,2] + 0.02), labels=rownames(scores.turn.phylogenetic))
abline(v=0,h=0,lty=2)

#phylogenetic nestedness
pcoa.nest.phylogenetic<-cmdscale(sqrt(nest.phylogenetic), eig = T)
scores.nest.phylogenetic<-vegan::scores(pcoa.nest.phylogenetic)
quartz()
plot(scores.nest.phylogenetic[,1],scores.nest.phylogenetic[,2],type="n",
     main="PCoA phylogenetic nestedness")
points(scores.nest.phylogenetic[1:10,1],scores.nest.phylogenetic[1:10,2],pch=19,cex=1.8)
points(scores.nest.phylogenetic[11:20,1],scores.nest.phylogenetic[11:20,2],pch=2,cex=1.8)
text(scores.nest.phylogenetic[,1],scores.nest.phylogenetic[,2], labels=rownames(scores.nest.phylogenetic))
abline(v=0,h=0,lty=2)

##################################
#adonis with environmental factors
#################################
env_std<- scale(x = env[, -10], center = T, scale = T)
env_std<- data.frame(cbind(env_std, bacia= as.factor(env$Bacia)))

#testing phylogenetic composition
total_phylogenetic_org_std<- as.matrix(as.dist(total.phylogenetic, diag = T, upper = T))[match(rownames(env_std), rownames(as.matrix(as.dist(total.phylogenetic, diag = T, upper = T)))),
                                                                                     match(rownames(env_std), colnames(as.matrix(as.dist(total.phylogenetic, diag = T, upper = T))))]


adonis_PhyloTotal_std<- adonis2(total_phylogenetic_org_std~Turbidez + pH + Condutividade + O2.porc + Altitude + Temp + Profundidade + Velocidade + bacia, 
                            data= env_std, permutations = 999, by = "margin")

turn_phylogenetic_org_std<- as.matrix(as.dist(turn.phylogenetic, diag = T, upper = T))[match(rownames(env_std), rownames(as.matrix(as.dist(turn.phylogenetic, diag = T, upper = T)))),
                                                                                   match(rownames(env_std), colnames(as.matrix(as.dist(turn.phylogenetic, diag = T, upper = T))))]

adonis_PhyloTurn_std<- adonis2(turn_phylogenetic_org_std~Turbidez + pH + Condutividade + O2.porc + Altitude + Temp + Profundidade + Velocidade + bacia, 
                      data= env_std, permutations = 999, by = "margin")


nest_phylogenetic_org_std<- as.matrix(as.dist(nest.phylogenetic, diag = T, upper = T))[match(rownames(env_std), rownames(as.matrix(as.dist(nest.phylogenetic, diag = T, upper = T)))),
                                                                                   match(rownames(env_std), colnames(as.matrix(as.dist(nest.phylogenetic, diag = T, upper = T))))]

adonis_PhyloNest_std<- adonis2(nest_phylogenetic_org_std~Turbidez + pH + Condutividade + O2.porc + Altitude + Temp + Profundidade + Velocidade + bacia, 
        data= env_std, permutations = 999, by = "margin")


#testing taxonomic composition
total_taxonomic_org_std<- as.matrix(as.dist(total.taxonomic, diag = T, upper = T))[match(rownames(env_std), rownames(as.matrix(as.dist(total.taxonomic, diag = T, upper = T)))),
                                                                                 match(rownames(env_std), colnames(as.matrix(as.dist(total.taxonomic, diag = T, upper = T))))]

adonis_TaxTotal_std<- adonis2(total_taxonomic_org_std~Turbidez + pH + Condutividade + O2.porc + Altitude + Temp + Profundidade + Velocidade + bacia, 
                             data= env_std, permutations = 999, by = "margin")


turn_taxonomic_org_std<- as.matrix(as.dist(turn.taxonomic, diag = T, upper = T))[match(rownames(env_std), rownames(as.matrix(as.dist(turn.taxonomic, diag = T, upper = T)))),
                                                                             match(rownames(env_std), colnames(as.matrix(as.dist(turn.taxonomic, diag = T, upper = T))))]

adonis_TaxTurn_std<- adonis2(turn_taxonomic_org_std~Turbidez + pH + Condutividade + O2.porc + Altitude + Temp + Profundidade + Velocidade + bacia, 
        data= env_std, permutations = 999, by = "margin")

nest_taxonomic_org_std<- as.matrix(as.dist(nest.taxonomic, diag = T, upper = T))[match(rownames(env_std), rownames(as.matrix(as.dist(nest.taxonomic, diag = T, upper = T)))),
                                                                             match(rownames(env_std), colnames(as.matrix(as.dist(nest.taxonomic, diag = T, upper = T))))]
adonis_TaxNest_std<- adonis2(nest_taxonomic_org_std~Turbidez + pH + Condutividade + O2.porc + Altitude + Temp + Profundidade + Velocidade + bacia, 
        data= env_std, permutations = 999, by = "margin")

