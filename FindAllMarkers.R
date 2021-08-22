## This function finds all the markers in every cluster type and will make three excel sheets in the output file
## Input : Seurat annotated object
## Output : Heatmap of top 5 differentially expressed genes per cluster type in the sample


library(Seurat)
library(dplyr)

data<-readRDS("data.rds")
pbmc.markers <- FindAllMarkers(data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

arranged_all <- subset(pbmc.markers,select = c(gene,p_val,p_val_adj,avg_log2FC,pct.1,pct.2,cluster))
arranged_all <- arranged_all %>% group_by(cluster)

dge_header<-read.csv("dge_header.txt",sep="\t")
wb <- createWorkbook()
addWorksheet(wb, sheetName = "Header Information")
writeData(wb, "Header Information",dge_header)
addWorksheet(wb, sheetName = "All.Genes")
writeData(wb, "All.Genes",arranged_all)

top10.markers <- pbmc.markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC)
	
top10.markers <- subset(top10.markers,select = c(gene,p_val,p_val_adj,avg_log2FC,pct.1,pct.2,cluster))
top10.markers <- top10.markers %>% dplyr::group_by(cluster)
	
addWorksheet(wb, sheetName = "Top 10 markers")
writeData(wb,"Top 10 markers",top10.markers)
saveWorkbook(wb, file = "21-H-4_D0_dge_results.xlsx", overwrite = T)
