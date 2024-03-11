library(edgeR)


##########GSE115821
counts = read.csv("./data/s_data/GSE115821/GSE115821_MGH_counts.csv.gz", sep = ",")
gene_names = counts[,1]
tpm <- data.frame(gene_id = gene_names)
for (i in 7:length(colnames(counts))){
  col = counts[[i]]
  len <- counts[[6]]
  rate <- col/len
  N <- sum(col) 
  TPMi <- (rate*1e6)/(sum(rate))
  print(sum(rate))
  TPMi <- pmax(TPMi,0) %>% as.data.frame()
  colnames(TPMi) <- colnames(counts)[i]
  tpm <- cbind(tpm,TPMi)
}
write.table(tpm, file = "./data/h_data/GSE115821/GSE115821_tpm.txt", quote = F, sep = "\t", row.names = F)
tpm = tpm[,-1]
tpm = log2(tpm + 1)
tpm = cbind(data.frame(gene_id = gene_names), tpm)
write.table(tpm, file = "./data/h_data/GSE115821/GSE115821_log2tpm.txt", quote = F, sep = "\t", row.names = F)


###########GSE78220
# FPKM2TPM
counts = read.csv("./data/s_data/GSE78220/GSE78220_fpkm.txt", sep = "\t")
gene_names = counts$Gene
tpm <- data.frame(gene_id = gene_names)
for (i in 2:length(colnames(counts))){
  TPMi <- FPKMtoTPM(counts[[i]]) %>% as.data.frame()
  colnames(TPMi) <- colnames(counts)[i]
  tpm <- cbind(tpm,TPMi)
}
write.table(tpm, file = "./data/h_data/GSE78220/GSE78220_tpm.txt", quote = F, sep = "\t", row.names = F)
tpm = tpm[,-1]
tpm = log2(tpm + 1)
tpm = cbind(data.frame(gene_id = gene_names), tpm)
write.table(tpm, file = "./data/h_data/GSE78220/GSE78220_log2tpm.txt", quote = F, sep = "\t", row.names = F)

FPKMtoTPM <- function(x) {
  return(exp(log(x) - log(sum(x)) + log(1e6)))
}

# patient labels
labels = read.csv("./data/s_data/GSE78220/GSE78220_patient.txt", sep = "\t")
flags = vector()
for (i in 1:length(rownames(labels))){
  if (labels[i,]$Response == "R"){
    flags = append(flags, 1)
  }
  else{
    flags = append(flags, 0)
  }
}
labels = cbind(labels, data.frame(flag = flags))
write.table(labels, file = "./data/h_data/GSE78220/GSE78220_patient.txt", quote = F, sep = "\t", row.names = F)



###############GSE91061
# FPKM2TPM
counts = read.csv("./data/s_data/GSE91061/GSE91061_fpkm.txt", sep = "\t")
# delete invalid gene
invalid_row_index = which(counts[, 1] == "-")
counts = counts[-invalid_row_index, ]
gene_names = counts[,1]
tpm <- data.frame(gene_id =gene_names)
for (i in 2:length(colnames(counts))){
  TPMi <- FPKMtoTPM(counts[[i]]) %>% as.data.frame()
  colnames(TPMi) <- colnames(counts)[i]
  tpm <- cbind(tpm,TPMi)
}
write.table(tpm, file = "./data/h_data/GSE91061/GSE91061_tpm.txt", quote = F, sep = "\t", row.names = F)
tpm = tpm[,-1]
tpm = log2(tpm + 1)
tpm = cbind(data.frame(gene_id = gene_names), tpm)
write.table(tpm, file = "./data/h_data/GSE91061/GSE91061_log2tpm.txt", quote = F, sep = "\t", row.names = F)


###############GSE126044
counts = read.csv("./data/s_data/GSE126044/GSE126044_counts.txt", sep = "\t")
gene_names = counts[,1]
tpm <- data.frame(gene_id = gene_names)
for (i in 2:length(colnames(counts))){
  col = counts[[i]]
  len <- counts[[6]]
  rate <- col/len
  N <- sum(col) 
  TPMi <- (rate*1e6)/(sum(rate))
  print(sum(rate))
  TPMi <- pmax(TPMi,0) %>% as.data.frame()
  colnames(TPMi) <- colnames(counts)[i]
  tpm <- cbind(tpm,TPMi)
}
write.table(tpm, file = "./data/h_data/GSE126044/GSE126044_tpm.txt", quote = F, sep = "\t", row.names = F)
tpm = tpm[,-1]
tpm = log2(tpm + 1)
tpm = cbind(data.frame(gene_id = gene_names), tpm)
write.table(tpm, file = "./data/h_data/GSE126044/GSE126044_log2tpm.txt", quote = F, sep = "\t", row.names = F)



###############GSE136961
tpm = read.csv("./data/s_data/GSE166449/GSE166449_tpm.txt", sep = "\t")
gene_names = tpm[,1]
tpm = tpm[,-1]
tpm = log2(tpm + 1)
tpm = cbind(data.frame(gene_id = gene_names), tpm)
write.table(tpm, file = "./data/h_data/GSE166449/GSE166449_log2tpm.txt", quote = F, sep = "\t", row.names = F)



###############GSE136961
counts = read.csv("./data/s_data/GSE136961/GSE136961_TPM.tsv", sep = "\t")
write.table(counts, file = "./data/h_data/GSE136961/GSE136961_tpm.txt", quote = F, sep = "\t", row.names = F)
