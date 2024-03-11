library(Seurat)

GSE179994_counts = readRDS("./data/s_data/GSE179994_all.Tcell.rawCounts.rds")
R.pre = c( 'P10.pre', 'P19.pre' , 'P29.pre', 'P30.pre', 'P33.pre', 'P35.pre')
R.post = c('P1.post.1','P1.post.3', 'P10.post.1', 'P13.post.1' ,'P19.post.1',
           'P29.post.1', 'P30.post.1', 'P33.post.1', 'P35.post.1')
R = append(R.pre, R.post)

NR.pre = c('P1.pre', 'P13.pre', 'P36.pre', 'P37.pre', 'P38.pre')
NR.post = c('P1.post.2', 'P13.post.2', 'P36.post.1', 'P37.post.1', 'P38.post.1')
NR = append(NR.pre, NR.post)
cell_meta = read.csv("./data/s_data/GSE179994_Tcell.metadata.tsv.gz", sep = "\t")
seurat_obj <- CreateSeuratObject(counts = GSE179994_counts, project = 'NSCLC',min.cells = 3 ,min.features = 200)
seurat_obj[["patient"]] = cell_meta$sample

# normalization
seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
seurat_obj <- ScaleData(seurat_obj, features = (rownames(seurat_obj)))
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 3000)
variable_genes = VariableFeatures(seurat_obj)
cells_id_R = vector()
cells_id_NR = vector()
for (i in 1:length(rownames(cell_meta)))
{
  if (cell_meta[i,]$sample %in% R)
  {
    cells_id_R = append(cells_id_R, cell_meta[i,]$cellid)
  }
  if (cell_meta[i,]$sample %in% NR)
  {
    cells_id_NR = append(cells_id_NR, cell_meta[i,]$cellid)
  }
}
seuobj_R = subset(seurat_obj, cells = cells_id_R)
seuobj_NR = subset(seurat_obj, cells = cells_id_NR)
exprMatrix_R = seuobj_R@assays$RNA@counts
exprMatrix_NR = seuobj_NR@assays$RNA@counts
temp = as.data.frame(exprMatrix_R)
df = temp[variable_genes, ]

if (!file.exists(paste0("../scGRN-L0_output/Data1-Qian et al. Cell Research 2020/weightList_PPCOR_", 
                        cellType, 
                        ".Rdata"))) {
  library("ppcor")
  message(paste0("---------- Now running PPCOR ----------"))
  ptm <- proc.time()
  inputExpr <- df
  geneNames <- rownames(inputExpr)
  rownames(inputExpr) <- c(geneNames)
  pcorResults <- pcor(x = t(as.matrix(inputExpr)), method = "spearman")
  DF <- data.frame(
    Gene1 = geneNames[c(row(pcorResults$estimate))], Gene2 = geneNames[c(col(pcorResults$estimate))],
    corVal = c(pcorResults$estimate), pValue = c(pcorResults$p.value)
  )
  outDF <- DF[order(DF$corVal, decreasing = TRUE), ]
  weightList <- outDF[-(1:8), -4]
  time <- proc.time() - ptm
  runningTime_PPCOR <- time[3]
  save(weightList, file = paste0("../scGRN-L0_output/Data1-Qian et al. Cell Research 2020/weightList_PPCOR_", 
                                 cellType, 
                                 ".Rdata")) ;rm(weightList)
} else {
  runningTime_PPCOR <- "0"
}

if (!file.exists(paste0("../scGRN-L0_output/Data1-Qian et al. Cell Research 2020/weightList_LEAP_", 
                        cellType, 
                        ".Rdata"))) {
  library("LEAP")
  message(paste0("---------- Now running LEAP ----------"))
  ptm <- proc.time()
  inputExpr <- matrixMC
  geneNames <- rownames(inputExpr)
  rownames(inputExpr) <- c()
  MAC_results <- MAC_counter(
    data = inputExpr,
    MAC_cutoff = 0,
    file_name = "temp",
    lag_matrix = FALSE,
    symmetric = FALSE
  )
  regulator <- geneNames[MAC_results[, "Row gene index"]]
  target <- geneNames[MAC_results[, "Column gene index"]]
  weight <- MAC_results[, "Correlation"]
  weightList <- data.frame(regulator, target, weight)
  time <- proc.time() - ptm
  runningTime_LEAP <- time[3]
  save(weightList, file = paste0("../scGRN-L0_output/Data1-Qian et al. Cell Research 2020/weightList_LEAP_", 
                                 cellType, 
                                 ".Rdata")) ;rm(weightList)
  file.remove("MAC_temp.csv")
}else {
  runningTime_LEAP <- "0"
}


R_regulate_network = outDF %>% subset(pValue<0.05) %>% subset(pValue > 0) %>% subset(abs(corVal) > 0.4)
write.csv(R_regulate_network, "./data/h_data/R_regulate_network.txt", row.names = F, quote = F)


GRN_PPCOR = read.table("./data/h_data/GRN_PPCOR.txt", sep = "\t", header = F)
