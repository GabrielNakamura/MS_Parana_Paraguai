###MS Taxonomic and phylogenetic turnover in Parana and Paraguay streams

####packages####
library(vegan)
library(ape)
library(betapart)
library(SYNCSA)
library(ade4)
library(geiger)
library(phytools)

#####data####
comm<- read.table(here::here("data", "processed", "comm_Parana-Paraguai.txt"), header= TRUE)
phylo<- read.tree(here::here("data", "processed", "tree_update_20-03-20.new"))
env<- read.table(here::here("data", "processed", "amb.txt"), header= T)
comm.presau<- ifelse(comm[,3:ncol(comm)]>=1,1,0) #incidence matrix
rownames(comm.presau)<- comm[,"pontos"] #naming the rows with the names of streams

###saving proccessed files####
write.table(here::here("data", "processed", "comm.txt"))
write.table(here::here("data", "processed", "env.txt"))

#######Taxonomic beta diversity#######
#between basins
beta.taxonomic<-beta.pair(comm.presau,index.family="sorensen") #taxonomic beta diversity
total.taxonomic<- beta.taxonomic$beta.sor
turn.taxonomic<-beta.taxonomic$beta.sim #turnover component of taxonomic beta diversity
nest.taxonomic<-beta.taxonomic$beta.sne #nestedness component of taxonomic beta diversity


#######Phylogenetic beta diversity#######
beta.phylogenetic<-phylo.beta.pair(comm.presau,tree = phylo,index.family="sorensen") #phylogenetic beta diversity
total.phylogenetic<- beta.phylogenetic$phylo.beta.sor
turn.phylogenetic<-beta.phylogenetic$phylo.beta.sim #turnover component of phylogenetic beta diversity
nest.phylogenetic<-beta.phylogenetic$phylo.beta.sne #nestedness component of phylogenetic beta diversity

#phylogenetic total beta div
pcoa.total.phylogenetic<- cmdscale(sqrt(total.phylogenetic), eig = T) #PCoA in phylogenetic turnover matrix
scores.total.phylogenetic<- scores(pcoa.total.phylogenetic)


##################################
#adonis with environmental factors
#################################
env_std<- scale(x = env[, -10], center = T, scale = T)
env_std<- data.frame(cbind(env_std, bacia= as.factor(env$Bacia)))

#testing phylogenetic composition
total_phylogenetic_org_std<- as.matrix(as.dist(total.phylogenetic, diag = T, upper = T))[match(rownames(env_std), rownames(as.matrix(as.dist(total.phylogenetic, diag = T, upper = T)))),
                                                                                     match(rownames(env_std), colnames(as.matrix(as.dist(total.phylogenetic, diag = T, upper = T))))]

adonis_PhyloTotal_std<- adonis2(total_phylogenetic_org_std~Turbidez + pH + Condutividade + O2.porc + Altitude + Temp + Profundidade + Velocidade + bacia, 
                            data= env_std, permutations = 999, by = "margin") #adonis with standardized variables


#testing phylogenetic turnover
turn_phylogenetic_org_std<- as.matrix(as.dist(turn.phylogenetic, diag = T, upper = T))[match(rownames(env_std), rownames(as.matrix(as.dist(turn.phylogenetic, diag = T, upper = T)))),
                                                                                   match(rownames(env_std), colnames(as.matrix(as.dist(turn.phylogenetic, diag = T, upper = T))))]


adonis_PhyloTurn_std<- adonis2(turn_phylogenetic_org_std~Turbidez + pH + Condutividade + O2.porc + Altitude + Temp + Profundidade + Velocidade + bacia, 
                      data= env_std, permutations = 999, by = "margin")


#testing phylogenetic nestedness
nest_phylogenetic_org_std<- as.matrix(as.dist(nest.phylogenetic, diag = T, upper = T))[match(rownames(env_std), rownames(as.matrix(as.dist(nest.phylogenetic, diag = T, upper = T)))),
                                                                                   match(rownames(env_std), colnames(as.matrix(as.dist(nest.phylogenetic, diag = T, upper = T))))]

adonis_PhyloNest_std<- adonis2(nest_phylogenetic_org_std~Turbidez + pH + Condutividade + O2.porc + Altitude + Temp + Profundidade + Velocidade + bacia, 
        data= env_std, permutations = 999, by = "margin")


####testing taxonomic composition#####
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


