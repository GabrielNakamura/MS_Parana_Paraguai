####loading function for edit data#####
source(here::here("R", "user_opt_printCat.R"))
source(here::here("R", "function_treedata_modif.R"))
source(here::here("R", "function_phylomatchFish_29-04-20.R"))
source(here::here("R", "function_tab_function_19-05-20.R"))

###community data####
comm<- read.table(here::here("data", "processed", "comm_Parana_Paraguai_16-05-20.txt"), header= TRUE)
comm.presau<- ifelse(comm[,3:ncol(comm)]>=1,1,0) #incidence matrix
rownames(comm.presau)<- comm$pontos

###phylogenetic tree####
taxon_data<- tab_function(comm.presau)
phylo<- phyloMatch(data = taxon_data)
write.tree(phylo, file = here::here("data", "processed", "data_phylo_19-05-20.txt"))

####environmental variables####
env<- read.table(here::here("data", "processed", "amb.txt"), header= T)
