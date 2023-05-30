countdata_new2 <- countdata_new %>% group_by(gene_name) %>% summarise_each(funs(sum)) 
