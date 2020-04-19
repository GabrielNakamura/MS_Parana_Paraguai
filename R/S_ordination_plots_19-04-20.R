####plotting total taxonomic beta-div####
pcoa.total.taxonomic<- cmdscale(sqrt(total.taxonomic), eig = T) #PCoA in taxonomic turnover matrix
scores.total.taxonomic<- scores(pcoa.total.taxonomic)
quartz()
plot(scores.total.taxonomic[,1],scores.total.taxonomic[,2],type="n",
     main="PCoA taxonomic beta-div")
points(scores.total.taxonomic[1:10,1],scores.total.taxonomic[1:10,2],pch=19,cex=1.8)
points(scores.total.taxonomic[11:20,1],scores.total.taxonomic[11:20,2],pch=2,cex=1.8)
text(scores.total.taxonomic[,1],scores.total.taxonomic[,2], labels=rownames(scores.total.taxonomic))
abline(v=0,h=0,lty=2)


#plotting turnover
pcoa.turn.taxonomic<- cmdscale(sqrt(turn.taxonomic), eig = T) #PCoA in taxonomic turnover matrix
scores.turn.taxonomic<- scores(pcoa.turn.taxonomic)
quartz()
plot(scores.turn.taxonomic[,1],scores.turn.taxonomic[,2],type="n",
     main="PCoA taxonomic turnover")
points(scores.turn.taxonomic[1:10,1],scores.turn.taxonomic[1:10,2],pch=19,cex=1.8)
points(scores.turn.taxonomic[11:20,1],scores.turn.taxonomic[11:20,2],pch=2,cex=1.8)
text(scores.turn.taxonomic[,1],scores.turn.taxonomic[,2], labels=rownames(scores.turn.taxonomic))
abline(v=0,h=0,lty=2)



#plotting nestedness
pcoa.nest.taxonomic<- cmdscale(sqrt(nest.taxonomic), eig = T)
scores.nest.taxonomic<-scores(pcoa.nest.taxonomic)
quartz()
plot(scores.nest.taxonomic[,1],scores.nest.taxonomic[,2],type="n",
     main="PCoA taxonomic nestedness")
points(scores.nest.taxonomic[1:10,1],scores.nest.taxonomic[1:10,2],pch=19,cex=1.8)
points(scores.nest.taxonomic[11:20,1],scores.nest.taxonomic[11:20,2],pch=2,cex=1.8)
text(scores.nest.taxonomic[,1],scores.nest.taxonomic[,2], labels=rownames(scores.nest.taxonomic))
abline(v=0,h=0,lty=2)

#plot with betadisper function
quartz()
plot(betadisper(d = total.taxonomic, group = comm$bacia, sqrt.dist = T), cex= 1.3, main= "Total taxonomic diversity", 
     xlab = "PCoA 1 (21.8 %)", ylab= "PCoA 2 (12 %)")
permutest(betadisper(d = total.taxonomic, group = comm$bacia, sqrt.dist = T))

quartz()
plot(betadisper(d = turn.taxonomic, group = comm$bacia, sqrt.dist = T), cex= 1.3, main= "Taxonomic turnover", 
     xlab = "PCoA 1 (35.4 %)", ylab= "PCoA 2 (15.9 %)")

quartz()
plot(betadisper(d = nest.taxonomic, group = comm$bacia, sqrt.dist = T), cex= 1.3, main= "Taxonomic nestedness", 
     xlab = "PCoA 1 (66 %)", ylab= "PCoA 2 (49 %)")
permutest(betadisper(d = nest.taxonomic, group = comm$bacia, sqrt.dist = T))


####ploting phylogenetic beta diversity####
quartz()
plot(scores.total.phylogenetic[,1],scores.total.phylogenetic[,2],type="n",
     main="PCoA phylogenetic beta-div")
points(scores.total.phylogenetic[1:10,1],scores.total.phylogenetic[1:10,2],pch=19,cex=1.8)
points(scores.total.phylogenetic[11:20,1],scores.total.phylogenetic[11:20,2],pch=2,cex=1.8)
text(scores.total.phylogenetic[,1], (scores.total.phylogenetic[,2] + 0.02), labels=rownames(scores.total.phylogenetic))
abline(v=0,h=0,lty=2)

#phylogenetic turnover
pcoa.turn.phylogenetic<- cmdscale(sqrt(turn.phylogenetic), eig = T) #PCoA in phylogenetic turnover matrix
scores.turn.phylogenetic<- scores(pcoa.turn.phylogenetic)
quartz()
plot(scores.turn.phylogenetic[,1],scores.turn.phylogenetic[,2],type="n",
     main="PCoA phylogenetic turnover")
points(scores.turn.phylogenetic[1:10,1],scores.turn.phylogenetic[1:10,2],pch=19,cex=1.8)
points(scores.turn.phylogenetic[11:20,1],scores.turn.phylogenetic[11:20,2],pch=2,cex=1.8)
text(scores.turn.phylogenetic[,1], (scores.turn.phylogenetic[,2] + 0.02), labels=rownames(scores.turn.phylogenetic))
abline(v=0,h=0,lty=2)

#phylogenetic nestedness
pcoa.nest.phylogenetic<-cmdscale(sqrt(nest.phylogenetic), eig = T)
scores.nest.phylogenetic<-scores(pcoa.nest.phylogenetic)
quartz()
plot(scores.nest.phylogenetic[,1],scores.nest.phylogenetic[,2],type="n",
     main="PCoA phylogenetic nestedness")
points(scores.nest.phylogenetic[1:10,1],scores.nest.phylogenetic[1:10,2],pch=19,cex=1.8)
points(scores.nest.phylogenetic[11:20,1],scores.nest.phylogenetic[11:20,2],pch=2,cex=1.8)
text(scores.nest.phylogenetic[,1],scores.nest.phylogenetic[,2], labels=rownames(scores.nest.phylogenetic))
abline(v=0,h=0,lty=2)

#betadisper
quartz()
plot(betadisper(d = total.phylogenetic, group = comm$bacia, sqrt.dist = T), cex= 1.3, main= "Total phylogenetic", 
     xlab = "PCoA 1 (21.8 %)", ylab= "PCoA 2 (14.8 %)")
permutest(betadisper(d = nest.taxonomic, group = comm$bacia, sqrt.dist = T))

quartz()
plot(betadisper(d = turn.phylogenetic, group = comm$bacia, sqrt.dist = T), cex= 1.3, main= "Phylogenetic turnover", 
     xlab = "PCoA 1 (28.6 %)", ylab= "PCoA 2 (18.3 %)")

quartz()
plot(betadisper(d = nest.phylogenetic, group = comm$bacia, sqrt.dist = T), cex= 1.3, main= "Phylogenetic nestedness", 
     xlab = "PCoA 1 (50 %)", ylab= "PCoA 2 (23 %)")

