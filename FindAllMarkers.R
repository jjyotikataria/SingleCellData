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


## Heatmap for 21 clusters for example:

png("21-H-4_D0_GEX.png",width=10600,height=9300,res=400) 	
#DoHeatmap(data, features = top10.markers$gene,draw.lines=F, size=1, angle=90 )+ scale_fill_gradientn(guide="colourbar",colors = c("blue", "white", "red"))+ theme(axis.title.x.top =element_text(angle = 90, vjust = 0.5, hjust=1),plot.title = element_text(hjust = 0.5))
colors30<- c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a","#ffff99","#b15928","#01665e","#c51b7d","#67001f","#053061","#4d4d4d","#f768a1","#276419","#8e0152","#7a0177","#80cdc1","#f1b6da","#e0e0e0","#f46d43","#878787","#045a8d","#fcc5c0","#1a1a1a","#4d004b")

DoHeatmap(data, features = top10.markers$gene,draw.lines=F, size=1, angle=90, group.colors=colors30)+ scale_fill_gradientn(guide="colourbar",colors = c("blue", "white", "red"))+ theme(axis.title.x.top =element_text(angle = 90, vjust = 0.5, hjust=1),plot.title = element_text(hjust = 0.5))
dev.off()
