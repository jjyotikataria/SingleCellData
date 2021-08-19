## This code does the following:
     ## - takes the cluster types annotated object from database, for example Human Primary Cell Atlas (HPCA) used here.
     ## - makes two-two pairs from the unique clusters for differential gene expression(dge)
     ## - performs dge between those combinations (one cluster vs the other)
     ## - performs the filtration for significant degs (differentially expressed genes)


## Input: Seurat object with cell type annotations from hpca
## Output: Dge between all the clusters in a sample

library(Seurat)
library(dplyr)
library(openxlsx)


hpca.data <- readRDS("../21-H-6_D5_GEX.rds") 
x <- as.data.frame(table(Idents(hpca.data)))

x_filt <- x[x[2] > 50,]
x_filt <- subset(x_filt, !grepl("Unknown", x_filt$Var1))

cell_types <- as.vector(x_filt[['Var1']])
combi <- as.data.frame(combn(cell_types, 2))
dge_df = data.frame(S.No. = numeric(),Groups = numeric(), Upregulated = numeric(), Downregulated= numeric())
S.No. <- data.frame("S.No." = seq.int(ncol(combi)))

print(combi)
for (k in seq(length(combi))){
  cell_type1 <- combi[1,k]
  cell_type2 <- combi[2,k] 
  #print("Dge started between",cell_type1 , " and ", cell_type2)
  print(cell_type1)
  print(cell_type2)
  DefaultAssay(hpca.data) <- "RNA"
  res <- FindMarkers(hpca.data, ident.1 = cell_type1, ident.2 = cell_type2)
  
  res2 <- cbind(gene = rownames(res), res)
  row.names(res2) <- NULL
  arranged <- subset(res2,select = c(gene,pct.1,pct.2,p_val,p_val_adj,avg_log2FC))
  arranged[is.na(arranged)]="N/A"
  arranged <- arranged %>% dplyr::arrange(avg_log2FC)
  file_name <- gsub(" ","-",paste0(cell_type1,"_vs_",cell_type2))
  file_name1 <- paste0("Differential_Expression_",file_name,".xlsx")
  print(file_name1)
  
  signi <- arranged[which(arranged$p_val_adj < 0.01),]
  signi1 <- signi[which(signi$avg_log2FC >=1 | signi$avg_log2FC <=-1),]
  signi1 <- signi1 %>% dplyr::arrange(avg_log2FC)
  
  up <- sum(signi1$avg_log2FC > 0)
  down <- sum(signi1$avg_log2FC < 0)
  dge_df <- rbind(dge_df, data.frame(file_name,up,down))
  
  dge_header <- read.csv("dge_header.txt",sep="\t")
  wb <- createWorkbook()
  addWorksheet(wb, sheetName = "Header Information")
  writeData(wb, "Header Information",dge_header)
  header <- read.csv("dge_header.txt",sep="\t")
  addWorksheet(wb, sheetName = "All.Genes")
  writeData(wb, "All.Genes",arranged)
  addWorksheet(wb, sheetName = "p_val_adj < 0.01 & FC +-2")
  writeData(wb,"p_val_adj < 0.01 & FC +-2",signi1)
  saveWorkbook(wb, file = file_name1, overwrite = T)
  
  upreg <- signi[which(signi$avg_log2FC > 0),]
  downreg <- signi[which(signi$avg_log2FC <0),]
  
  upreg <- upreg %>% dplyr::arrange(desc(avg_log2FC))
  downreg <- downreg %>% dplyr::arrange(avg_log2FC)
  
  top25up <- upreg[1:25,]
  top25down <- downreg[1:25,]
  genes_heatmap <-  c(top25up$gene,top25down$gene)
  

   print(file_name)
 
 
}
  
final_table <- cbind(S.No., dge_df)
colnames(final_table)[2] <- "Cluster comparisons"
colnames(final_table)[3] <- "Upregulated genes"
colnames(final_table)[4] <- "Downregulated genes"
final_table<-final_table[-1]
write.csv(final_table,"dge_counts.csv")

