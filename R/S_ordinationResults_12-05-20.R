#plot with betadisper function
quartz()
layout.show(layout(mat = matrix(c(1, 2, 3,
                      4, 5, 6),
nrow = 2, ncol= 3, byrow = T)
))
plot(betadisper(d = total.taxonomic, group = comm$bacia, sqrt.dist = T), cex= 1.5, main= "Total", 
     xlab = "PCoA 1 (20.9 %)", ylab= "PCoA 2 (13 %)")
plot(betadisper(d = turn.taxonomic, group = comm$bacia, sqrt.dist = T), cex= 1.5, main= "Turnover", 
     xlab = "PCoA 1 (29.1 %)", ylab= "PCoA 2 (17.0 %)")
plot(betadisper(d = nest.taxonomic, group = comm$bacia, sqrt.dist = T), cex= 1.5, main= "Nestedness", 
     xlab = "PCoA 1 (64 %)", ylab= "PCoA 2 (34 %)")
plot(betadisper(d = total.phylogenetic, group = comm$bacia, sqrt.dist = T), cex= 1.5, main= "Total", 
     xlab = "PCoA 1 (22.0 %)", ylab= "PCoA 2 (14.6 %)")
plot(betadisper(d = turn.phylogenetic, group = comm$bacia, sqrt.dist = T), cex= 1.5, main= "Turnover", 
     xlab = "PCoA 1 (25.4 %)", ylab= "PCoA 2 (18.8 %)")
plot(betadisper(d = nest.phylogenetic, group = comm$bacia, sqrt.dist = T), cex= 1.5, main= "Nestedness", 
     xlab = "PCoA 1 (50 %)", ylab= "PCoA 2 (23 %)")

(pcoa.total.taxonomic$eig)/sum(pcoa.total.taxonomic$eig)
