## This code does the following:
     ## - takes the cluster types annotated from database, for example Human Primary Cell Atlas (HPCA) used here.
     ## - makes two-two pairs from the unique clusters for differential gene expression(dge)
     ## - performs dge between those combinations (one cluster vs the other)
     ## - performs the filtration for significant degs (differentially expressed genes)

x <- as.data.frame(table(Idents(hpca.data)))
x_filt <- x[x[2] > 20,]
cell_types <- as.vector(x_filt[['Var1']])
combi <- as.data.frame(combn(cell_types, 2))
for (k in seq(length(combi))){
  cell_type1 <- combi[1,k]
  cell_type2 <- combi[2,k]  
  print(cell_type1)
  print(cell_type2)
  res <- FindMarkers(hpca.data, ident.1 = cell_type1, ident.2 = cell_type2,test.use = "MAST")
  res2 <- cbind(gene = rownames(res), res)
  row.names(res2) <- NULL
  arranged <- subset(res2,select = c(gene,pct.1,pct.2,p_val,p_val_adj,avg_log2FC))
  arranged[is.na(arranged)]="N/A"
  arranged <- arranged %>% dplyr::arrange(avg_log2FC)
  file_name <- gsub(" ","-",paste0(cell_type1,"_vs_",cell_type2))
  file_name1 <- paste0("Differential_Expression_",file_name,".xlsx")
  print(file_name1)
  write.xlsx(arranged, file=file_name1,sheetName="All.Genes",append=FALSE,row.names=FALSE)
  signi <- arranged[which(arranged$p_val_adj < 0.01),]
  signi1 <- signi[which(signi$avg_log2FC >=1 | signi$avg_log2FC <=-1),]
  signi1 <- signi1 %>% dplyr::arrange(avg_log2FC)
  
  
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
}
