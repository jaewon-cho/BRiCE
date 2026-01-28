#obj: seurat object
#meta information
#orig.ident
#clonotype
#individual
#condition: unperturb, perturb
#usage: clonal_seurat_process(user_seurat_object, "user_file_name")


library(Seurat)
library(SingleCellExperiment)
library(zellkonverter)


clonal_seurat_process<-function(obj, save_name){
individual_list <- unique(obj$individual)

obj$clonality <- 1
obj$expand <- "no"

filter_cell <- c()
for (j in c(1:length(individual_list))){
	ind <- individual_list[j]
	print(ind)
	a1 <- obj[,obj$individual == ind & obj$condition == "unperturb"]
	p <- table(a1$clonotype[a1$condition == "unperturb"])
q <-which(p>1)
normal_exp <- names(q)
tmp_norm <- names(a1$orig.ident)[a1$clonotype %in% normal_exp]
filter_cell<- c(filter_cell, tmp_norm)


a1 <- obj[,obj$individual == ind & obj$condition == "perturb"]
a1 <- a1[,!a1$clonotype %in% normal_exp]
p1 <- table(a1$clonotype)
temp <- a1$clonality
names(temp) <- a1$clonotype
temp[names(p1)] <- p1
names(temp) <- names(a1$orig.ident)

m <- match(names(temp), names(obj$orig.ident))
obj$clonality[m] <- temp
	
}

obj1<- obj[,!names(obj$orig.ident) %in% filter_cell]
obj1$expand[obj1$clonality > 1]  <- "yes"
saveRDS(obj1, paste(save_name, ".rds", sep = ""))
a2<- as.SingleCellExperiment(obj1)
writeH5AD(a2, paste(save_name, ".h5ad", sep = ""))

}

