#plot with betadisper function
quartz()
layout.show(layout(mat = matrix(c(1, 2, 3,
                      4, 5, 6),
nrow = 2, ncol= 3, byrow = T)
))
plot(betadisper(d = total.taxonomic, group = comm$bacia, sqrt.dist = T), cex= 1.3, main= "Total", 
     xlab = "PCoA 1 (21.8 %)", ylab= "PCoA 2 (12 %)")
plot(betadisper(d = turn.taxonomic, group = comm$bacia, sqrt.dist = T), cex= 1.3, main= "Turnover", 
     xlab = "PCoA 1 (35.4 %)", ylab= "PCoA 2 (15.9 %)")
plot(betadisper(d = nest.taxonomic, group = comm$bacia, sqrt.dist = T), cex= 1.3, main= "Nestedness", 
     xlab = "PCoA 1 (66 %)", ylab= "PCoA 2 (49 %)")
plot(betadisper(d = total.phylogenetic, group = comm$bacia, sqrt.dist = T), cex= 1.3, main= "Total", 
     xlab = "PCoA 1 (21.8 %)", ylab= "PCoA 2 (14.8 %)")
plot(betadisper(d = turn.phylogenetic, group = comm$bacia, sqrt.dist = T), cex= 1.3, main= "Turnover", 
     xlab = "PCoA 1 (28.6 %)", ylab= "PCoA 2 (18.3 %)")
plot(betadisper(d = nest.phylogenetic, group = comm$bacia, sqrt.dist = T), cex= 1.3, main= "Nestedness", 
     xlab = "PCoA 1 (50 %)", ylab= "PCoA 2 (23 %)")
