#############################################################################
#############Análises convergência funcional Paraná - Paraguai###############
#############################################################################

#pacotes necessários
library(vegan)
library(ape)
library(betapart)
library(SYNCSA)

#conjunto de dados necessários
traits.Parana<-read.table("clipboard",header=TRUE)
traits.Paraguai<-read.table("clipboard",header=TRUE)
comu.Parana<-read.table("clipboard",header=TRUE)
comu.Paraguai<-read.table("clipboard",header=TRUE)
comu.todas<-read.table("clipboard",header=TRUE)
comu.todas.pad<-decostand(comu.todas,"total",MARGIN=1) #padronizando pelo total
Parana.Paraguai<-ifelse(read.table("clipboard",header=TRUE)>=1,1,0)
phylo_para_parag<- read.tree("phylo_parana_paraguai.txt")
windows()
plot(phylo_para_parag, cex= 0.7)

#################################################################
#diversidade Beta para dimensão taxonômica
#################################################################
#Para ordenação entre as comunidades e comparação das configurações

#1 - utilizar a função beta.pair() com todas as comunidades juntas em uma 
#matriz

presau.comu.todas<-ifelse(comu.todas>=1,1,0)
beta.total<-beta.pair(presau.comu.todas,index.family="sorensen")

turn.comu.todas<-beta.total$beta.sim #componente de turnover
sqrt.turn.comu.todas<-sqrt(turn.comu.todas) #raiz quadrada para diminuir propriedades não euclidianas presentes nesta matriz

mean(turn.comu.todas)
mean(sqrt.turn.comu.todas) #componente de turnover com os dados transformados na raiz quadrada
sd(turn.comu.todas)
sd(sqrt.turn.comu.todas)

nest.comu.todas<-beta.total$beta.sne #componente de aninhamento
sqrt.nest.comu.todas<-sqrt(nest.comu.todas)
mean(nest.comu.todas)
mean(sqrt.nest.comu.todas)
sd(nest.comu.todas)
sd(sqrt.nest.comu.todas)

mean(beta.total$beta.sor)
mean(sqrt(beta.total$beta.sor))
sd(beta.total$beta.sor)
sd(sqrt(beta.total$beta.sor))

      #calcular para bacias separadas
comu.Paraguai<-ifelse(comu.Paraguai>=1,1,0)
comu.Parana<-ifelse(comu.Parana>=1,1,0)
beta.tax.Parana<-beta.pair(comu.Parana,index.family="sorensen")

mean(beta.tax.Parana$beta.sor) #div beta taxonomic total
mean(sqrt(beta.tax.Parana$beta.sor)) #corrigido pela raiz quadrada
sd(beta.tax.Parana$beta.sor) #desvio div beta taxonomic total
sd(sqrt(beta.tax.Parana$beta.sor))

turn.tax.Parana<-beta.tax.Parana$beta.sim #componente de turnover do Paraná
mean(turn.tax.Parana) #média de turnover do Paraná
mean(sqrt(turn.tax.Parana)) #média dos valores transformados pela raiz quadrada
sd(turn.tax.Parana) #desvio do Paraná
sd(sqrt(turn.tax.Parana))

nest.tax.Parana<-beta.tax.Parana$beta.sne
mean(nest.tax.Parana)
mean(sqrt(nest.tax.Parana))
sd(nest.tax.Parana)
sd(sqrt(nest.tax.Parana))

total.tax.Parana<-beta.tax.Parana$beta.sor
mean(total.tax.Parana)
sd(total.tax.Parana)

mean(turn.tax.Parana/total.tax.Parana)
mean(nest.tax.Parana/total.tax.Parana)
sd(turn.tax.Parana/total.tax.Parana)

        
beta.tax.Paraguai<-beta.pair(comu.Paraguai,index.family="sorensen")
turn.tax.Paraguai<-beta.tax.Paraguai$beta.sim #componente de turnover do Paraguai
mean(turn.tax.Paraguai) #média de turnover do Paraguai
sd(turn.tax.Paraguai) #desvio do Paraguai

nest.tax.Paraguai<-beta.tax.Paraguai$beta.sne
mean(nest.tax.Paraguai)
sd(nest.tax.Paraguai)

total.tax.Paraguai<-beta.tax.Paraguai$beta.sor
mean(total.tax.Paraguai) #media total div beta taxonomic paraguai
sd(total.tax.Paraguai) #desvio total div beta taxonomic paraguai

mean(turn.tax.Paraguai/total.tax.Paraguai)
sd(turn.tax.Paraguai/total.tax.Paraguai)
mean(nest.tax.Paraguai/total.tax.Paraguai)

              #usando somente duas comunidades Paraná e Paraguai
beta.tax.Parana.Paraguai<-beta.pair(Parana.Paraguai,index.family="sorensen")
(beta.tax.Parana.Paraguai$beta.sim)/(beta.tax.Parana.Paraguai$beta.sor) #porcentagem
                                              #da diversidade beta decorrente de turnover
(beta.tax.Parana.Paraguai$beta.sne)/(beta.tax.Parana.Paraguai$beta.sor) #porcentagem da
                                              #diversidade beta decorrente de aninhamento

NMDS.turn<-metaMDS(turn.comu.todas,autotransform=FALSE,k=3) #NMDS utilizando as 
        #distâncias baseadas em turnover
scores.turn<-scores(NMDS.turn)
windows()
plot(scores.turn[,1],scores.turn[,2],type="n",main="turnover",
     xlim=c(-0.6,0.6)) #gráfico baseado no Turnover
abline(v=0,h=0,lty=2)
points(scores.turn[1:10,1],scores.turn[1:10,2],pch=19)
points(scores.turn[11:20,1],scores.turn[11:20,2],pch=2)
text(scores.turn[,1],scores.turn[,2],labels=rownames(scores.turn))

#usando PCoA para gerar os gráficos com os dados transformados na raiz quadrada
pcoa.sqrt.tax.turn<-cmdscale(sqrt.turn.comu.todas)
scores.pcoa.sqrt.tax.turn<-scores(pcoa.sqrt.tax.turn)
windows()
plot(scores.pcoa.sqrt.tax.turn[,1],scores.pcoa.sqrt.tax.turn[,2],type="n",
     main="PCoA taxonomic turnover")
points(scores.pcoa.sqrt.tax.turn[1:10,1],scores.pcoa.sqrt.tax.turn[1:10,2],pch=19,cex=1.8)
points(scores.pcoa.sqrt.tax.turn[11:20,1],scores.pcoa.sqrt.tax.turn[11:20,2],pch=2,cex=1.8)
text(scores.pcoa.sqrt.tax.turn[,1],scores.pcoa.sqrt.tax.turn[,2], labels=rownames(scores.pcoa.sqrt.tax.turn))
abline(v=0,h=0,lty=2)

NMDS.nest<-metaMDS(nest.comu.todas,autotransform=FALSE,k=3) #NMDS utilizando
        #as distâncias baseadas em turnover
scores.nest<-scores(NMDS.nest)
windows()
plot(scores.nest[,1],scores.nest[,2],type="n",main="nestedness",
     xlim=c(-0.5,0.4)) #gráfico baseado no aninhamento
abline(v=0,h=0,lty=2)
points(scores.nest[1:10,1],scores.nest[1:10,2],pch=19)
points(scores.nest[11:20,1],scores.nest[11:20,2],pch=2)
text(scores.nest[,1],scores.nest[,2],labels=rownames(scores.nest))

pcoa.sqrt.tax.nest<-cmdscale(sqrt.nest.comu.todas)
scores.pcoa.sqrt.tax.nest<-scores(pcoa.sqrt.tax.nest)
windows()
plot(scores.pcoa.sqrt.tax.nest[,1],scores.pcoa.sqrt.tax.nest[,2],type="n",
     main="PCoA taxonomic nestedness")
points(scores.pcoa.sqrt.tax.nest[1:10,1],scores.pcoa.sqrt.tax.nest[1:10,2],pch=19,cex=1.8)
points(scores.pcoa.sqrt.tax.nest[11:20,1],scores.pcoa.sqrt.tax.nest[11:20,2],pch=2,cex=1.8)
text(scores.pcoa.sqrt.tax.nest[,1],scores.pcoa.sqrt.tax.nest[,2], labels=rownames(scores.pcoa.sqrt.tax.nest))
abline(v=0,h=0,lty=2)

NMDS.tax<-metaMDS(beta.total$beta.sor,autotransform=FALSE,k=3)
scores.tax.total<-scores(NMDS.tax)  
?plot

#2- comparar através de uma análise de procrustes as configurações provenientes
#da ordenação a partir da div beta taxonômica e da div beta funcional
?protest
summary(procrustes(scores.turn,scores.nest))
protest(scores.turn[,1:2],scores.nest[,1:2],scale=TRUE) #taxonomic
                                        #turnover and nestdness
plot(protest(scores.turn[-11,1:2],scores.func.turn[,1:2],scale=TRUE)) #taxonomic
                              #turnover e functional turnover
protest(scores.func.nest[,1:2],scores.func.turn[,1:2],scale=TRUE)
r<-sqrt(1-0.9853)
text()
perm

#3- calcular a diversidade beta multi-sites
multi.beta.tax<-beta.multi(presau.comu.todas,index.family="sorensen") #para todos os riachos
multi.beta.tax.Parana<-beta.multi(presau.comu.todas[1:10,],
                                  index.family="sorensen")
multi.beta.tax.Paraguai<-beta.multi(presau.comu.todas[11:20,],index.family=
                                      "sorensen")
?beta.pair
###########################Fim div beta para dim tax#######################


#################################
#div Beta para dimensão funcional
#################################

#tabelas necessárias
traits.todos<-read.table("clipboard",header=TRUE)
names(traits.todos)


#1 - Calcular distância de Gower entre os atributos
      #primeiro separar os atributos
traits.todos
diet.traits<-traits.todos[,c(4,6)] #traços quantitativos de dieta
habitat.traits<-prep.binary(traits.todos[,1:3],col.blocks=3) #traços binarios de habitat
reprod.traits<-as.data.frame(traits.todos[,5]) #traços quantitativos de reprodução

      #criando a matriz de distancia
ktab1<-ktab.list.df(list(diet.traits,habitat.traits,reprod.traits))

distraits<-dist.ktab(ktab1,c("Q","B","Q"),"scaledBYrange") #distância entre
                        #os as espécies com base nos traços


#2 - calcular o espaço multidimensional baseado nos eixos de uma 
#PCoA
PCoA.traits<-cmdscale(distraits,k=3,eig=TRUE)
summary(PCoA.traits)
(PCoA.traits$eig[1]+PCoA.traits$eig[2]+PCoA.traits$eig[3])/
  sum(PCoA.traits$eig) #quantidade de
                        #variação explicada pelos três primeiros eixos

#3 - utilizar a função functional.beta.pair() para calcular a beta diversidade
#e seus componente
rowSums(presau.comu.todas) #verificar quantas espécies tem cada um dos riachos
       
  #remover o córrego apa das análises e utilizar apenas duas dimensões funcionais
presau.comu1<-presau.comu.todas[-11,] #excluindo córrego Apa por ter uma espécie apenas
  #calcular a beta diversidade funcional
beta.functional.total<-functional.beta.pair(presau.comu1,
                                            PCoA.traits$points[,1:2],
                                            index.family="sorensen")

beta.functional.nest<-beta.functional.total$funct.beta.sne #componente de
                                                #aninhamento
sqrt.beta.functional.nest<-sqrt(beta.functional.nest) #transformação da raiz quadrada para diminuir as propriedades não-
                                                          #euclidianas presentes nesta matrix

beta.functional.turn<-beta.functional.total$funct.beta.sim #componente de 
                                                  #turnover
sqrt.beta.functional.turn<-sqrt(beta.functional.turn) #transformação da raiz quadrada para diminuir as propriedades não-
                                                      #euclidianas presentes nesta matrix


mean(beta.functional.turn)
sd(beta.functional.turn)
mean(beta.functional.nest)
sd(beta.functional.nest)
mean(beta.functional.total$funct.beta.sor)
sd(beta.functional.total$funct.beta.sor)

    #Diversidade funcional beta apenas para Paraná
beta.func.Parana<-functional.beta.pair(comu.Parana,
                     PCoA.traits$points[colnames(comu.Parana),c(1:2)],
                     index.family="sorensen")
turn.func.Parana<-beta.func.Parana$funct.beta.sim
nest.func.Parana<-beta.func.Parana$funct.beta.sne
total.func.Parana<-beta.func.Parana$funct.beta.sor
mean(turn.func.Parana)
sd(turn.func.Parana)
mean(nest.func.Parana)
sd(nest.func.Parana)
mean(total.func.Parana)
sd(total.func.Parana)
mean(turn.func.Parana/total.func.Parana)
mean(nest.func.Parana/total.func.Parana)

      #Diversidade funcional beta apenas para Paraguai
beta.func.Paraguai<-functional.beta.pair(comu.Paraguai[-1,],
                                       PCoA.traits$points[colnames(comu.Paraguai),c(1:2)],
                                       index.family="sorensen")
turn.func.Paraguai<-beta.func.Paraguai$funct.beta.sim
nest.func.Paraguai<-beta.func.Paraguai$funct.beta.sne
total.func.Paraguai<-beta.func.Paraguai$funct.beta.sor
mean(turn.func.Paraguai) #média de turnover
sd(turn.func.Paraguai)
mean(nest.func.Paraguai) #média de aninhamento
sd(nest.func.Paraguai)
mean(total.func.Paraguai) #média total
sd(total.func.Paraguai)
mean(turn.func.Paraguai/total.func.Paraguai) #representatividade de turnover
                                             #em %
mean(nest.func.Paraguai/total.func.Paraguai) #representatividade de aninhamento
                                             #em %


    #calculando levando em conta apenas Paraná e Paraguai
beta.func.Parana.Paraguai<-functional.beta.pair(Parana.Paraguai,
                                                PCoA.traits$points[,1:2],
                                                index.family="sorensen")
par(bty="l")
stripchart(c(mean(nest.func.Paraguai),mean(turn.func.Paraguai))~
             c("Nestedness","Turnvover"),pch=19,vertical=TRUE,at=c(1.3,1.7),
           ylim=c(-0.1,1))
?stripchart
arrows(x0=c(1.3,1.7),y0=c((mean(nest.func.Paraguai)-sd(nest.func.Paraguai))
                         ,(mean(turn.func.Paraguai)-sd(turn.func.Paraguai))),
      x1=c(1.3,1.7),y1=c((mean(nest.func.Paraguai)+sd(nest.func.Paraguai))
                         ,(mean(turn.func.Paraguai)+sd(turn.func.Paraguai))),
      angle=90,length=0.04,code=3)
?arrows
#4 - ordenar as comunidades através de um NMDS

NMDS.func.turn<-metaMDS(beta.functional.turn,autotransform=FALSE,k=3) #NMDS utilizando as 
#distâncias baseadas em turnover
scores.func.turn<-scores(NMDS.func.turn)
windows()
plot(scores.func.turn[,1],scores.func.turn[,2],type="n",
     main="functional turnover",xlim=c(-0.4,0.7)) #gráfico baseado no Turnover
abline(v=0,h=0,lty=2)
points(scores.func.turn[1:10,1],scores.func.turn[1:10,2],pch=19)
points(scores.func.turn[11:19,1],scores.func.turn[11:19,2],pch=2)
text(scores.func.turn[,1],scores.func.turn[,2],
     labels=rownames(scores.func.turn))

#########################################################

NMDS.sqrt.func.turn<-metaMDS(sqrt.beta.functional.turn,autotransform=FALSE,k=3) #NMDS utilizando as 
#distâncias baseadas em turnover
scores.sqrt.func.turn<-scores(NMDS.sqrt.func.turn)
windows()
plot(scores.sqrt.func.turn[,1],scores.sqrt.func.turn[,2],type="n",
     main="functional turnover",xlim=c(-0.4,0.7)) #gráfico baseado no Turnover transformado sqrt()
abline(v=0,h=0,lty=2)
points(scores.sqrt.func.turn[1:10,1],scores.sqrt.func.turn[1:10,2],pch=19)
points(scores.sqrt.func.turn[11:19,1],scores.sqrt.func.turn[11:19,2],pch=2)
text(scores.sqrt.func.turn[,1],scores.sqrt.func.turn[,2],
     labels=rownames(scores.sqrt.func.turn))

#com PCoA para turnover funcional

pcoa.sqrt.funct.turn<-cmdscale(sqrt.beta.functional.turn)
?cmdscale
scores.pcoa.sqrt.funct<-scores(pcoa.sqrt.funct.turn)
windows()
plot(scores.pcoa.sqrt.funct[,1],scores.pcoa.sqrt.funct[,2],type="n", 
     main="PCoA functional turnover")
points(scores.pcoa.sqrt.funct[1:10,1],scores.pcoa.sqrt.funct[1:10,2],pch=19,cex=1.8)
points(scores.pcoa.sqrt.funct[11:19,1],scores.pcoa.sqrt.funct[11:19,2],pch=2,cex=1.8)
text(scores.pcoa.sqrt.funct[,1],scores.pcoa.sqrt.funct[,2], labels=rownames(scores.pcoa.sqrt.funct))
abline(v=0,h=0,lty=2)

#########################################################


NMDS.func.nest<-metaMDS(beta.functional.nest,autotransform=FALSE,k=3) #NMDS utilizando
#as distâncias baseadas em aninhamento
scores.func.nest<-scores(NMDS.func.nest)
windows()
plot(scores.func.nest[,1],scores.func.nest[,2],type="n",main="nestedness",
     xlim=c(-0.6,0.6)) #gráfico baseado no aninhamento
abline(v=0,h=0,lty=2)
points(scores.func.nest[1:10,1],scores.func.nest[1:10,2],pch=19)
points(scores.func.nest[11:19,1],scores.func.nest[11:19,2],pch=2)
text(scores.func.nest[,1],scores.func.nest[,2],
     labels=rownames(scores.func.nest))

#com PCoA para aninhamento funcional
pcoa.sqrt.funct.nest<-cmdscale(sqrt.beta.functional.nest)
scores.pcoa.sqrt.funct.nest<-scores(pcoa.sqrt.funct.nest)
windows()
plot(scores.pcoa.sqrt.funct.nest[,1],scores.pcoa.sqrt.funct.nest[,2],type="n",
     main="PCoA functional nestedness")
points(scores.pcoa.sqrt.funct.nest[1:10,1],scores.pcoa.sqrt.funct.nest[1:10,2],pch=19,cex=1.8)
points(scores.pcoa.sqrt.funct.nest[11:19,1],scores.pcoa.sqrt.funct.nest[11:19,2],pch=2,cex=1.8)
text(scores.pcoa.sqrt.funct.nest[,1],scores.pcoa.sqrt.funct.nest[,2], labels=rownames(scores.pcoa.sqrt.funct.nest))
abline(v=0,h=0,lty=2) #com PCoA para aninhamento funcional





################Fim div beta para dim funcional#############################

######################################
#div Beta para a dimensão filogenética
######################################

phylo_para_parag
phylo_para_parag$tip.label==colnames(presau.comu1)
colnames(presau.comu1)<-colnames(presau.comu1)[match(phylo_para_parag$tip.label,colnames(presau.comu1))]
match(phylo_para_parag$tip.label, colnames(presau.comu1))
phylo_sobra_parana<- drop.tip(phy = phylo_para_parag, tip = colnames(comu.Parana))
phylo_parana<- drop.tip(phy = phylo_para_parag, tip = phylo_sobra_parana$tip.label)
colnames(comu.Parana)<-colnames(comu.Parana)[match(phylo_parana$tip.label, colnames(comu.Parana))]
comu.Paraguai.sApa<-comu.Paraguai[-1,]
phylo_sobra_paraguai<- drop.tip(phy = phylo_para_parag, tip = colnames(comu.Paraguai.sApa))
phylo_paraguai<- drop.tip(phy = phylo_para_parag, tip = phylo_sobra_paraguai$tip.label)
length(phylo_paraguai$tip.label)
length(colnames(comu.Paraguai.sApa)
colnames(comu.Paraguai.sApa)<-colnames(comu.Paraguai.sApa)[match(phylo_paraguai$tip.label, colnames(comu.Paraguai.sApa))]


beta.phylo.nest<-beta.phylo.total$phylo.beta.sne #componente de aninhamento
sqrt.beta.phylo.nest<-sqrt(beta.phylo.nest) #transformação da raiz quadrada para diminuir as propriedades não-euclidianas presentes nesta matrix

beta.phylo.turn<-beta.phylo.total$phylo.beta.sim #componente de turnover
sqrt.beta.phylo.turn<-sqrt(beta.phylo.turn) #transformação da raiz quadrada para diminuir as propriedades não-euclidianas presentes nesta matrix


mean(beta.phylo.turn)
sd(beta.phylo.turn)
mean(beta.phylo.nest)
sd(beta.phylo.nest)
mean(beta.phylo.total$phylo.beta.sor)
sd(beta.phylo.total$phylo.beta.sor)
(mean(beta.phylo.turn)/mean(beta.phylo.total$phylo.beta.sor)) #prop de turnover
(mean(beta.phylo.nest)/mean(beta.phylo.total$phylo.beta.sor)) #prop de aninhamento
mean(beta.phylo.total)
#Diversidade filogenetica beta apenas para Paraná
beta.phylo.Parana<-phylo.beta.pair(comu.Parana,
                                       phylo_parana,
                                       index.family="sorensen")
turn.phylo.Parana<-beta.phylo.Parana$phylo.beta.sim
nest.phylo.Parana<-beta.phylo.Parana$phylo.beta.sne
total.phylo.Parana<-beta.phylo.Parana$phylo.beta.sor
mean(turn.phylo.Parana)
sd(turn.phylo.Parana)
mean(nest.phylo.Parana)
sd(nest.phylo.Parana)
mean(total.phylo.Parana)
sd(total.phylo.Parana)
mean(turn.phylo.Parana/total.phylo.Parana)
mean(nest.phylo.Parana/total.phylo.Parana)

#Diversidade filogenetica beta apenas para Paraguai
beta.phylo.Paraguai<-phylo.beta.pair(comu.Paraguai.sApa,
                                         phylo_paraguai,
                                         index.family="sorensen")
turn.phylo.Paraguai<-beta.phylo.Paraguai$phylo.beta.sim
nest.phylo.Paraguai<-beta.phylo.Paraguai$phylo.beta.sne
total.phylo.Paraguai<-beta.phylo.Paraguai$phylo.beta.sor
mean(turn.phylo.Paraguai) #média de turnover
sd(turn.phylo.Paraguai)
mean(nest.phylo.Paraguai) #média de aninhamento
sd(nest.phylo.Paraguai)
mean(total.phylo.Paraguai) #média total
sd(total.phylo.Paraguai)
mean(turn.phylo.Paraguai/total.phylo.Paraguai) #representatividade de turnover em %
mean(nest.phylo.Paraguai/total.phylo.Paraguai) #representatividade de aninhamento em %
library(vegan)

#4 - ordenar as comunidades através de um NMDS
NMDS.phylo.turn<-metaMDS(beta.phylo.turn,autotransform=FALSE,k=3) #NMDS utilizando as 
#distâncias baseadas em turnover
scores.phylo.turn<-scores(NMDS.phylo.turn)
windows()
plot(scores.phylo.turn[,1],scores.phylo.turn[,2],type="n",
     main="phylogenetic turnover",xlim=c(-0.4,0.7)) #gráfico baseado no Turnover
abline(v=0,h=0,lty=2)
points(scores.phylo.turn[1:10,1],scores.phylo.turn[1:10,2],pch=19)
points(scores.phylo.turn[11:19,1],scores.phylo.turn[11:19,2],pch=2)
text(scores.phylo.turn[,1],scores.phylo.turn[,2],
     labels=rownames(scores.phylo.turn))

#usando a matriz P para ordenar as comunidades 

matrixP.para.parag<-matrix.p(presau.comu1, dist.spp = cophenetic(phylo_para_parag))
windows()
biplot(prcomp(matrixP.para.parag$matrix.P, method = "bray"), center = FALSE, scale. = FALSE)
str(matrixP.para.parag)]
windows()
biplot(cmdscale(d = vegdist(matrixP.para.parag$matrix.P)))
pcoa.phylo.para.parag<-cmdscale(d = vegdist(matrixP.para.parag$matrix.P))
loads.phylo.para.parag<- wascores(x = pcoa.phylo.para.parag, w = matrixP.para.parag$matrix.P)
windows()
plot(pcoa.phylo.para.parag[,1], pcoa.phylo.para.parag[,2], type= "n")
points(pcoa.phylo.para.parag[,1], pcoa.phylo.para.parag[,2])
text(pcoa.phylo.para.parag[,1], pcoa.phylo.para.parag[,2], labels= rownames(pcoa.phylo.para.parag))
text(loads.phylo.para.parag[,1], loads.phylo.para.parag[,2], labels= rownames(loads.phylo.para.parag), cex= 0.7)

#pcoa com a matriz de turnover 
Pcoa_phyloTurn<-cmdscale(beta.phylo.turn, k=2)
windows()
plot(Pcoa_phyloTurn[,1], Pcoa_phyloTurn[,2], type= "n")
points(Pcoa_phyloTurn[1:10,1], Pcoa_phyloTurn[1:10,2], pch=19)
points(Pcoa_phyloTurn[11:19,1], Pcoa_phyloTurn[11:19,2], pch=2)
text(Pcoa_phyloTurn[,1], Pcoa_phyloTurn[,2], labels= rownames(Pcoa_phyloTurn))

#pcoa para matriz de aninhamento
Pcoa_phyloNest<-cmdscale(beta.phylo.nest, k=2)
windows()
plot(Pcoa_phyloNest[,1], Pcoa_phyloNest[,2], type= "n")
points(Pcoa_phyloNest[1:10,1], Pcoa_phyloNest[1:10,2], pch=19)
points(Pcoa_phyloNest[11:19,1], Pcoa_phyloNest[11:19,2], pch=2)
text(Pcoa_phyloNest[,1], Pcoa_phyloNest[,2], labels= rownames(Pcoa_phyloNest))
library(ape)
plot(phylo_para_parag)
#########################################################

NMDS.sqrt.func.turn<-metaMDS(sqrt.beta.functional.turn,autotransform=FALSE,k=3) #NMDS utilizando as 
#distâncias baseadas em turnover
scores.sqrt.func.turn<-scores(NMDS.sqrt.func.turn)
windows()
plot(scores.sqrt.func.turn[,1],scores.sqrt.func.turn[,2],type="n",
     main="functional turnover",xlim=c(-0.4,0.7)) #gráfico baseado no Turnover transformado sqrt()
abline(v=0,h=0,lty=2)
points(scores.sqrt.func.turn[1:10,1],scores.sqrt.func.turn[1:10,2],pch=19)
points(scores.sqrt.func.turn[11:19,1],scores.sqrt.func.turn[11:19,2],pch=2)
text(scores.sqrt.func.turn[,1],scores.sqrt.func.turn[,2],
     labels=rownames(scores.sqrt.func.turn))

#com PCoA para turnover funcional

pcoa.sqrt.funct.turn<-cmdscale(sqrt.beta.functional.turn)
?cmdscale
scores.pcoa.sqrt.funct<-scores(pcoa.sqrt.funct.turn)
windows()
plot(scores.pcoa.sqrt.funct[,1],scores.pcoa.sqrt.funct[,2],type="n", 
     main="PCoA functional turnover")
points(scores.pcoa.sqrt.funct[1:10,1],scores.pcoa.sqrt.funct[1:10,2],pch=19,cex=1.8)
points(scores.pcoa.sqrt.funct[11:19,1],scores.pcoa.sqrt.funct[11:19,2],pch=2,cex=1.8)
text(scores.pcoa.sqrt.funct[,1],scores.pcoa.sqrt.funct[,2], labels=rownames(scores.pcoa.sqrt.funct))
abline(v=0,h=0,lty=2)

#########################################################


NMDS.func.nest<-metaMDS(beta.functional.nest,autotransform=FALSE,k=3) #NMDS utilizando
#as distâncias baseadas em aninhamento
scores.func.nest<-scores(NMDS.func.nest)
windows()
plot(scores.func.nest[,1],scores.func.nest[,2],type="n",main="nestedness",
     xlim=c(-0.6,0.6)) #gráfico baseado no aninhamento
abline(v=0,h=0,lty=2)
points(scores.func.nest[1:10,1],scores.func.nest[1:10,2],pch=19)
points(scores.func.nest[11:19,1],scores.func.nest[11:19,2],pch=2)
text(scores.func.nest[,1],scores.func.nest[,2],
     labels=rownames(scores.func.nest))

#com PCoA para aninhamento funcional
pcoa.sqrt.funct.nest<-cmdscale(sqrt.beta.functional.nest)
scores.pcoa.sqrt.funct.nest<-scores(pcoa.sqrt.funct.nest)
windows()
plot(scores.pcoa.sqrt.funct.nest[,1],scores.pcoa.sqrt.funct.nest[,2],type="n",
     main="PCoA functional nestedness")
points(scores.pcoa.sqrt.funct.nest[1:10,1],scores.pcoa.sqrt.funct.nest[1:10,2],pch=19,cex=1.8)
points(scores.pcoa.sqrt.funct.nest[11:19,1],scores.pcoa.sqrt.funct.nest[11:19,2],pch=2,cex=1.8)
text(scores.pcoa.sqrt.funct.nest[,1],scores.pcoa.sqrt.funct.nest[,2], labels=rownames(scores.pcoa.sqrt.funct.nest))
abline(v=0,h=0,lty=2) #com PCoA para aninhamento funcional

###########################################Fim Beta Diversidade Filogenética########################################################

################################################################
#para testar os valores de diversidade beta
################################################################
#função para aleatorização da matriz e calculo de div.beta 
library(picante)
library(betapart)
library(vegan)


rnd.beta.pair.tax<-function(x){
  rnd.comu<-randomizeMatrix(x,null.model="frequency")
  beta.rnd<-beta.pair(rnd.comu,index.family="sorensen")
  return(list(turnover=beta.rnd$beta.sim, nest=beta.rnd$beta.sne))
  
}
test.replicate<-replicate(999,rnd.beta.pair.tax(Parana.Paraguai))

hist(as.numeric(test.replicate[1,]))
abline(v=c(0.35,0.45),lty=2)

rnd.beta.pair.funct<-function(Parana.Paraguai){
  rnd.comu<-randomizeMatrix(Parana.Paraguai,null.model="frequency")
  beta.funct.rand<-functional.beta.pair(rnd.comu,PCoA.traits$points[,1:2],
                                        index.family="sorensen")
  return(list(turnover=beta.funct.rnd$beta.sim, nest=beta.funct.rnd$beta.sne))
}
test.replicate.funct<-replicate(999,rnd.beta.pair.funct(Parana.Paraguai))


windows()
hist(as.numeric(test.replicate[2,]))
is.list()

rnd.resul<-numeric(0)
for(i in 1:ncol(test.replicate)){
    rnd.resul[i]<-mean(test.replicate[2,i]$nest)
}
?quantile
?unlist
as.matrix(test.replicate)
quantile(unlist(test.replicate[1,]),probs=c(0.025,0.75))
############################################################### 

#################################################################
#Para testar convergência funcional, calcular a div beta apenas
#com espécies de ocorrência endêmica
#################################################################

comu.endemico<-ifelse(read.table("clipboard",header=TRUE)>=1,1,0)

traits.endemico<-read.table("clipboard",header=TRUE)

beta.pair(comu.endemico,index.family="sorensen")

functional.beta.pair(comu.endemico,
                     PCoA.traits$points[colnames(comu.endemico),1:2])

organize_endemic<-organize.syncsa(comm=comu.endemico,dist.spp=cophenetic(phylo_para_parag))
str(organize_endemic)                     
phylo.beta.pair(organize_endemic$community, phylo_para_parag) #pode ser visto como uma medida de redundancia filogenetica na diversidade beta
###########################################################
#testando as relações entre os valores de diversidade beta 
###########################################################

#usar um teste de mantel para isso

Mantel.turnfunc.turntax<-mantel(as.dist(as.matrix(turn.comu.todas)[-11,-11]),beta.functional.turn)
                #teste de mantel entre turnover funcional e taxonomico
windows()
plot(beta.functional.turn,as.dist(as.matrix(turn.comu.todas)[-11,-11]),
     pch=19,xlab="beta functional turnover",ylab="beta taxonomic turnover",cex=1.3)
                #gráfico teste de mantel turnover funcional e taxonomico


Mantel.nestfunc.nesttax<-mantel(as.dist(as.matrix(nest.comu.todas)[-11,-11]),
                                beta.functional.nest) #teste de mantel entre 
                                                  #aninhamento funcional e taxonomico
windows()
plot(beta.functional.nest,as.dist(as.matrix(nest.comu.todas)[-11,-11]),
     pch=19,xlab="beta functional nestedness",
     ylab="beta functional nestedness") #grafico teste mantel aninhamento
                                  #funcional e taxonomico


Mantel.nestfunc.turnfunc<-mantel(beta.functional.turn,beta.functional.nest)
windows()
plot(beta.functional.turn,beta.functional.nest,pch=19,
     xlab="beta functional turnover",
     ylab="beta functional nestedness") #grafico teste mantel turnover X aninhamento
                                      #funcional

Mantel.nesttax.turntax<-mantel(nest.comu.todas,turn.comu.todas) #teste de mantel
                                      #entre nestedness X turnover taxonomico
windows()
plot(turn.comu.todas,nest.comu.todas,pch=19,
     xlab="beta taxonomic turnover",
     ylab="beta taxonomic nestedness") #grafico teste mantel turnover X aninhamento
#funcional


Mantel.nesttax.functurn<-mantel(as.dist(as.matrix(nest.comu.todas)[-11,-11]),
                                beta.functional.turn)
windows()
plot(beta.functional.turn,as.dist(as.matrix(nest.comu.todas)[-11,-11]),
     pch=19,xlab="beta functional turnover",
     ylab="beta taxonomic nestedness") #grafico teste mantel turnover funcional
              #e aninhamento taxonomico

Mantel.turntax.funcnest<-mantel(as.dist(as.matrix(turn.comu.todas)[-11,-11]),
                                beta.functional.nest)
windows()
plot(beta.functional.nest,as.dist(as.matrix(turn.comu.todas)[-11,-11]),
     pch=19,xlab="beta functional nestedness",
     ylab="beta taxonomic turnover") #grafico teste mantel aninhamento funcional
                                #e turnover taxonomico

windows()
plot(beta.functional.nest,as.dist(as.matrix(nest.comu.todas)[-11,-11]),
     pch=19,xlab="beta functional nestedness",
     ylab="beta taxonomic nestedness") #grafico teste mantel aninhamento funcional
#e aninhamento taxonomico

#plotar todas as relações em um mesmo quadro com os respectivos 
#valores de p



#gráfico teste de mantel turnover funcional e taxonomico

#####################################################################################
#investigando a relação de fatores ambientais com beta diversidade - BDEHR Hypothesis
#####################################################################################

#matriz de fatores ambientais

amb.todos<-read.table("clipboard",header=TRUE)
colnames(amb.todos)

amb.todos<-decostand(amb.todos,"standardize") #padronizando para média 0
                                              # e mesma variância

dist.amb.jac<-vegdist(amb.todos,"jaccard")
dist.amb<-dist(amb.todos,"euclidean")
PCA.amb<-prcomp(amb.todos,scale.=FALSE)
summary(PCA.amb)

PCA.amb.summary<-scores(PCA.amb)[,1:5]  #fatores ambientais sintéticos

dist.PCA.amb<-dist(PCA.amb.summary[,1:5],"euclidean")

Mantel.functotal.env<-mantel(beta.functional.total$funct.beta.sor,
    as.dist(as.matrix(dist.PCA.amb)[rownames(as.matrix(beta.functional.total$funct.beta.sor)),
                                    colnames(as.matrix(beta.functional.total$funct.beta.sor))]))
windows()
plot(as.dist(as.matrix(dist.PCA.amb)[rownames(as.matrix(beta.functional.turn)),
                                     colnames(as.matrix(beta.functional.turn))]
     ),beta.functional.total$funct.beta.sor,
     xlab="environmental distance",ylab="functional beta diversity",
     pch=19) #plotando distância ambiental contra diversidade beta funcional total


Mantel.functurn.env<-mantel(beta.functional.turn,as.dist(as.matrix(dist.PCA.amb)[rownames(as.matrix(beta.functional.turn)),
                                          colnames(as.matrix(beta.functional.turn))]))
windows()
plot(as.dist(as.matrix(dist.PCA.amb)[rownames(as.matrix(beta.functional.turn)),colnames(as.matrix(beta.functional.turn))]),
     beta.functional.turn,xlab="environmental distance",ylab="beta functional turnover", pch=19)

Mantel.functnest.env<-mantel(beta.functional.nest,as.dist(as.matrix(dist.PCA.amb)[rownames(as.matrix(beta.functional.nest)),
             colnames(as.matrix(beta.functional.nest))]))
windows()
plot(as.dist(as.matrix(dist.PCA.amb)[rownames(as.matrix(beta.functional.nest)),colnames(as.matrix(beta.functional.nest))]
     ),beta.functional.nest, xlab="environmental distance",ylab="beta functional nestedness",pch=19)

Mantel.taxtotal.env<-mantel(beta.total$beta.sor,
                     as.dist(as.matrix(dist.PCA.amb)[rownames(as.matrix(beta.total$beta.sor)),
                                                            colnames(as.matrix(beta.total$beta.sor))]))
windows()
plot(as.dist(as.matrix(dist.PCA.amb)[rownames(as.matrix(nest.comu.todas)),
                                     colnames(as.matrix(nest.comu.todas))]
     ),beta.total$beta.sor,
     xlab="environmental distance",ylab="taxonomic beta diversity",
     pch=19)

Mantel.taxnest.env<-mantel(nest.comu.todas,
                           as.dist(as.matrix(dist.PCA.amb)[rownames(as.matrix(nest.comu.todas)),
                           colnames(as.matrix(nest.comu.todas))]))
windows()
plot(as.dist(as.matrix(dist.PCA.amb)[rownames(as.matrix(nest.comu.todas)), 
                                 colnames(as.matrix(nest.comu.todas))]
     ), nest.comu.todas,xlab="environmental distance",ylab="beta taxonomic nestedness",
     pch=19)

Mantel.taxturn.env<-mantel(turn.comu.todas,as.dist(as.matrix(dist.PCA.amb)[rownames(as.matrix(turn.comu.todas)),
                            colnames(as.matrix(turn.comu.todas))]))
windows()
plot(as.dist(as.matrix(dist.amb)[rownames(as.matrix(turn.comu.todas)),
                                 colnames(as.matrix(turn.comu.todas))]
     ),turn.comu.todas,xlab="environmental distance",ylab="beta taxonomic turnvover",
     pch=19)

#investigando a relação dos fatores ambientais para cada Bacia

#Paraná Basin

mantel.taxturn.env.Parana<-mantel(turn.tax.Parana,as.dist(as.matrix(dist.amb)[rownames(as.matrix(turn.tax.Parana)),
                                                                              colnames(as.matrix(turn.tax.Parana))]))

mantel.taxnest.env.Parana<-mantel(nest.tax.Parana,as.dist(as.matrix(dist.amb)[rownames(as.matrix(nest.tax.Parana)),
                                                                              colnames(as.matrix(nest.tax.Parana))]))

mantel.turnfun.env.Parana<-mantel(turn.func.Parana,as.dist(as.matrix(dist.amb)[rownames(as.matrix(turn.func.Parana)),
                                                                              colnames(as.matrix(turn.func.Parana))]))

mantel.nestfun.env.Parana<-mantel(nest.func.Parana,as.dist(as.matrix(dist.amb)[rownames(as.matrix(nest.func.Parana)),
                                                                               colnames(as.matrix(nest.func.Parana))]))

#Paraguai Basin
mantel.taxturn.env.Paraguai<-mantel(turn.tax.Paraguai,as.dist(as.matrix(dist.amb)[rownames(as.matrix(turn.tax.Paraguai)),
                                                                              colnames(as.matrix(turn.tax.Paraguai))]))

mantel.taxnest.env.Paraguai<-mantel(nest.tax.Paraguai,as.dist(as.matrix(dist.amb)[rownames(as.matrix(nest.tax.Paraguai)),
                                                                                  colnames(as.matrix(nest.tax.Paraguai))]))

mantel.functurn.env.Paraguai<-mantel(turn.func.Paraguai,as.dist(as.matrix(dist.amb)[rownames(as.matrix(turn.func.Paraguai)),
                                                                                  colnames(as.matrix(turn.func.Paraguai))]))

mantel.nestturn.env.Paraguai<-mantel(nest.tax.Paraguai,as.dist(as.matrix(dist.amb)[rownames(as.matrix(nest.tax.Paraguai)),
                                                                                  colnames(as.matrix(nest.tax.Paraguai))]))

#apenas na bacia do Paraguai foram encontradas relações significativas entre heterogeineidade ambiental e 
#beta diversidade.

#######################################################################
#acessando quanto cada variável contribui na explicação da  variação
#da diversidade beta
#######################################################################

?envfit #usando a função envfit para verificar quais fatores ambientais
  #melhor explicam os padrões de diversidade beta funcional e taxonomic
library(vegan)
library(betapart)

fit.amb.tax.turn<-envfit(pcoa.sqrt.tax.turn[rownames(amb.todos),],amb.todos,
       permutations=999) #taxonomic turnover
fit.amb.tax.turn$vectors

fit.amb.tax.nest<-envfit(pcoa.sqrt.tax.nest,amb.todos,permutations=999) #taxonomic nest
fit.amb.tax.nest$vectors

#funcional
fit.amb.func.turn<-envfit(pcoa.sqrt.funct.turn,amb.todos[rownames(pcoa.sqrt.funct.turn),]
       ,permutations=999) #functional turnover
fit.amb.func.turn$vectors


fit.amb.func.nest<-envfit(pcoa.sqrt.funct.nest,amb.todos[rownames(pcoa.sqrt.funct.nest),]
       ,permutations=999) #functional nestedness
fit.amb.func.nest$vectors

#comparando os valores de r2 entre as diferentes matrizes
r2.env<-read.table("clipboard",header=TRUE)
cor(r2.env,method="pearson")
is.matrix(r2.env)
is.data.frame(r2.env)
windows()
plot(r2.env,pch=19,cex=1.25)
windows()
plot(r2.env[,1],r2.env[,2],pch=19,cex=1.25,xlab="functional turnover",
     ylab="functional nestednesss")
windows()
plot(r2.env[,1],r2.env[,3],pch=19,cex=1.25,xlab="functional turnover",
     ylab="taxonomic turnover")
windows()
plot(r2.env[,1],r2.env[,3],pch=19,cex=1.25,xlab="functional turnover",
     ylab="taxonomic turnover")



windows()
plot(r2.env,pch=19,cex=1.25)

##################################################################
#Mantel Parcial com os diferentes aspectos de diversidade Beta
# e distância ambiental
##################################################################
?mantel.partial
beta.functional.nest
beta.functional.total
nrow(as.matrix(beta.functional.turn))
str(beta.functional.nest)
nrow(as.matrix(beta.functional.total$funct.beta.sim))

mant.partial<-mantel.partial(beta.functional.nest,
                             as.dist(as.matrix(dist.PCA.amb)[rownames(as.matrix(beta.functional.nest)),
                            colnames(as.matrix(beta.functional.nest))]),
                            beta.functional.total$funct.beta.sim)

mantel(beta.functional.nest,
       as.dist(as.matrix(dist.PCA.amb)[rownames(as.matrix(beta.functional.nest)),
                                       colnames(as.matrix(beta.functional.nest))]))




