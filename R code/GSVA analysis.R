library(GSVA)



#########################################immune cell infiltration score####################################
immuneGeneSet = read.csv("./data/s_data/TILs_genes.txt", sep = "\t", header = F)
rownames(immuneGeneSet) = immuneGeneSet[,1]
immuneGeneSet = immuneGeneSet[,-1]

cells_name = rownames(immuneGeneSet)
l = list()
for (name in cells_name){
  x = as.character(immuneGeneSet[name,])
  x = x[nchar(x)!=0]
  x = as.character(x)
  l[[name]] = x
}

##########Gide dataset#########
expr = read.table("./data/h_data/Gide/Gide_norm3.txt", header = T, sep = "\t")
rownames(expr) = expr[, 1]
expr = expr[, -1]
expr_matrix = as.matrix(expr)
ssgsea = gsva(expr_matrix, l, method = "ssgsea", kcdf="Gaussian", abs.ranking=T)
write.table(ssgsea, file = "./data/h_data/Gide/Gide_TILs.txt", quote = F, sep = "\t")



##########Riaz dataset#########
expr = read.table("./data/h_data/Riaz/Riaz_norm3.txt", header = T, sep = "\t")
rownames(expr) = expr[, 1]
expr = expr[, -1]
expr_matrix = as.matrix(expr)
ssgsea = gsva(expr_matrix, l, method = "ssgsea", kcdf="Gaussian", abs.ranking=T)
write.table(ssgsea, file = "./data/h_data/Riaz/Riaz_TILs.txt", quote = F, sep = "\t")


##########Kim dataset#########
expr = read.table("./data/s_data/Kim/expression_mRNA.norm3.txt", header = T, sep = "\t")
rownames(expr) = expr[, 1]
expr = expr[, -1]
expr_matrix = as.matrix(expr)
ssgsea = gsva(expr_matrix, l, method = "ssgsea", kcdf="Gaussian", abs.ranking=T)
write.table(ssgsea, file = "./data/h_data/Kim/Kim_TILs.txt", quote = F, sep = "\t")

##########计GSE91061 dataset#########
expr = read.table("./data/h_data/GSE91061/GSE91061_tpm.txt", header = T, sep = "\t")
rownames(expr) = expr[, 1]
expr = expr[, -1]
expr_matrix = as.matrix(expr)
ssgsea = gsva(expr_matrix, l, method = "ssgsea", kcdf="Gaussian", abs.ranking=T)
write.table(ssgsea, file = "./data/h_data/GSE91061/GSE91061_TILs.txt", quote = F, sep = "\t")



##########GSE78220 dataset#########
expr = read.table("./data/h_data/GSE78220/GSE78220_tpm.txt", header = T, sep = "\t")
rownames(expr) = expr[, 1]
expr = expr[, -1]
expr_matrix = as.matrix(expr)
ssgsea = gsva(expr_matrix, l, method = "ssgsea", kcdf="Gaussian", abs.ranking=T)
write.table(ssgsea, file = "./data/h_data/GSE78220/GSE78220_TILs.txt", quote = F, sep = "\t")



##########Liu dataset#########
expr = read.table("./data/h_data/Liu/Liu_tpm.txt", header = T, sep = "\t")
rownames(expr) = expr[, 1]
expr = expr[, -1]
expr_matrix = as.matrix(expr)
ssgsea = gsva(expr_matrix, l, method = "ssgsea", kcdf="Gaussian", abs.ranking=T)
write.table(ssgsea, file = "./data/h_data/Liu/Liu_TILs.txt", quote = F, sep = "\t")



############################################pathways###################################

gs_data = readLines("./data/s_data/c2.all.v7.2.symbols.gmt")
geneSets = list()
for (line in gs_data){
  temp_list = strsplit(line, "\t")
  temp_list[[1]] = temp_list[[1]][-2]
  pathway_name = temp_list[[1]][1]
  temp_list[[1]] = temp_list[[1]][-1]
  geneSets[[pathway_name]] = temp_list[[1]]
} 


############GSE91061
expr = read.csv("./data/h_data/GSE91061/GSE91061_tpm.txt", sep = "\t")
rownames(expr) = expr[, 1]
expr = expr[, -1]
gsva_matrix<- gsva(
  expr = data.matrix(expr),
  gset.idx.list = geneSets,
  method='ssgsea',
  kcdf='Gaussian',
  abs.ranking=TRUE)
write.table(gsva_matrix, file = "./data/h_data/GSE91061_ssgsea.txt", quote = FALSE, sep="\t")


############GSE78220
expr = read.csv("./data/h_data/GSE78220/GSE78220_tpm.txt", sep = "\t")
rownames(expr) = expr[, 1]
expr = expr[, -1]
gsva_matrix<- gsva(
  expr = data.matrix(expr),
  gset.idx.list = geneSets,
  method='ssgsea',
  kcdf='Gaussian',
  abs.ranking=TRUE)
write.table(gsva_matrix, file = "./data/h_data/GSE78220_ssgsea.txt", quote = FALSE, sep="\t")

############GSE166449
expr = read.csv("../data/h_data/GSE166449/GSE166449_tpm.txt", sep = "\t")
rownames(expr) = expr[, 1]
expr = expr[, -1]
gsva_matrix<- gsva(
  expr = data.matrix(expr),
  gset.idx.list = geneSets,
  method='ssgsea',
  kcdf='Gaussian',
  abs.ranking=TRUE)
write.table(gsva_matrix, file = "./data/h_data/GSE166449_ssgsea.txt", quote = FALSE, sep="\t")

############Liu
expr = read.csv("./data/h_data/Liu2/Liu_norm3.txt", sep = "\t")
rownames(expr) = expr[, 1]
expr = expr[, -1]
gsva_matrix<- gsva(
  expr = data.matrix(expr),
  gset.idx.list = geneSets,
  method='ssgsea',
  kcdf='Gaussian',
  abs.ranking=TRUE)
write.table(gsva_matrix, file = "./data/h_data/Liu2/Liu_ssgsea.txt", quote = FALSE, sep="\t")


############Gide
expr = read.csv("./data/h_data/Gide/Gide_norm3.txt", sep = "\t")
rownames(expr) = expr[, 1]
expr = expr[, -1]
gsva_matrix<- gsva(
  expr = data.matrix(expr),
  gset.idx.list = geneSets,
  method='ssgsea',
  kcdf='Gaussian',
  abs.ranking=TRUE)
write.table(gsva_matrix, file = "./data/h_data/Gide/Gide_ssgsea.txt", quote = FALSE, sep="\t")



############Riaz
expr = read.csv("./data/h_data/Riaz/Riaz_norm3.txt", sep = "\t")
rownames(expr) = expr[, 1]
expr = expr[, -1]
gsva_matrix<- gsva(
  expr = data.matrix(expr),
  gset.idx.list = geneSets,
  method='ssgsea',
  kcdf='Gaussian',
  abs.ranking=TRUE)
write.table(gsva_matrix, file = "./data/h_data/Riaz/Riaz_ssgsea.txt", quote = FALSE, sep="\t")


############Kim
expr = read.csv("./data/s_data/Kim/expression_mRNA.norm3.txt", sep = "\t")
rownames(expr) = expr[, 1]
expr = expr[, -1]
gsva_matrix<- gsva(
  expr = data.matrix(expr),
  gset.idx.list = geneSets,
  method='ssgsea',
  kcdf='Gaussian',
  abs.ranking=TRUE)
write.table(gsva_matrix, file = "./data/h_data/Kim/Kim_ssgsea.txt", quote = FALSE, sep="\t")
