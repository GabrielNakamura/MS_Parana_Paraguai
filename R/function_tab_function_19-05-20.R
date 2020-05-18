tab_function<- function(comm){
  library(rfishbase)
  data("fishbase")
  list_genus<- fishbase[match(sub("_.*", "", colnames(comm)), fishbase$Genus), c("Genus", "Family", "Order")]
  list_local<- data.frame(s= colnames(comm), f= list_genus$Family, o= list_genus$Order)
  not_found<- list_local[which(rowSums(is.na(list_local)) > 0), ]
  
  ####manual insertion##
  print_cat_Family<- function(not_found_fishtree){
    cat("tell the Family of ", not_found_fishtree)
    cat("\n")
  }
  print_cat_Order<- function(not_found_fishtree){
    cat("tell the Order of ", not_found_fishtree)
    cat("\n")
  }
  for(i in 1:length(not_found$s)){
    spp_family<- readline(prompt = print_cat_Family(not_found_fishtree = not_found$s[i]))
    spp_order<-  readline(prompt = print_cat_Order(not_found_fishtree = not_found$s[i]))
    list_local[which(rowSums(is.na(list_local)) > 0)[i], "f"]<- spp_family
    list_local[which(rowSums(is.na(list_local)) > 0)[i], "o"]<- spp_order
  }
  sub_Perciformes<- which(list_local$o == "Perciformes")
  list_local[sub_Perciformes, "o"]<- "Cichliformes"
  list_local
} 
