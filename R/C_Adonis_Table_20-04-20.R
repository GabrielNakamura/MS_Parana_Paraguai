####writing Adonis results for Table 1 and 2#####

#Adonis table for taxonomic beta diversity
Table_adonis_res_taxonomic<- matrix(unlist(lapply(list(adonis_TaxTotal_std, adonis_TaxTurn_std, adonis_TaxNest_std), 
                                                  function(x) {
                                                    tables_tmp<- cbind(x$Df, x$R2, x$`Pr(>F)`)
                                                    tables_tmp
                                                  }
)
),
nrow= nrow(adonis_TaxTotal_std), ncol= 3*3, byrow= F, 
dimnames = list(rownames(adonis_TaxTotal_std), 
                rep(c("Df", "R2", "P value"), 3)
)
)

write.csv(Table_adonis_res_taxonomic, file = here::here("output", "tables", "Table_01_adonis_taxonomic.csv")) #saving adonis result for taxonomic beta diversity

#Adonis table for phylogenetic beta diversity
Table_adonis_res_phylogenetic<- matrix(unlist(lapply(list(adonis_PhyloTotal_std, adonis_PhyloTurn_std, adonis_PhyloNest_std), 
                                                     function(x) {
                                                       tables_tmp<- cbind(x$Df, x$R2, x$`Pr(>F)`)
                                                       tables_tmp
                                                     }
)
),
nrow= nrow(adonis_PhyloTotal_std), ncol= 3*3, byrow= F, 
dimnames = list(rownames(adonis_PhyloTotal_std), 
                rep(c("Df", "R2", "P value"), 3)
)
)

write.csv(Table_adonis_res_phylogenetic, file = here::here("output", "tables", "Table_01_adonis_phylogenetic.csv")) #saving adonis result for taxonomic beta diversity
