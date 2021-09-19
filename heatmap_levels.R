---
title: "scsorter annotations"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```


```{r cars, echo=FALSE}

setwd("F:/NEW/Projects/FINAL_VALIDATION")
full_path<-getwd()
pacman::p_load(dplyr, Seurat, patchwork, scSorter, SingleR, celldex, ReactomeGSA,writexl, ggplot2,stringr,DoubletFinder,cowplot,RColorBrewer,gridExtra,ggpubr,purrr)

files<-list.files(path = full_path, pattern = "*_scSorter.rds$",full.names = TRUE)

for (i in files){
data<-readRDS(i)
pbmc.markers <- FindAllMarkers(data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>%
    group_by(cluster) %>%
    top_n(n = 5, wt = avg_log2FC) -> top5.markers
sam_name<-gsub("*_scSorter.rds","",i)
sam_name<-gsub(".*/","",sam_name)
title_name<-paste0("Top 5 differentially expressed markers for ",sam_name)

check<-as.data.frame(sort(table(Idents(data))))

colors30<- c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a","#ffff99","#b15928","#01665e","#c51b7d","#67001f","#053061","#4d4d4d","#f768a1","#276419","#8e0152","#7a0177","#80cdc1","#f1b6da","#e0e0e0","#f46d43","#878787","#045a8d","#fcc5c0","#1a1a1a","#4d004b")

check_sorted<-data.frame(Var1=character(),Freq=numeric())

for(i in 1:(nrow(check)/2)){

if(i!=(nrow(check)-i+1)){check_sorted<-rbind(check_sorted,c(as.character(check$Var1[nrow(check)-i+1]),check$Freq[nrow(check)-i+1]))}
check_sorted<-rbind(check_sorted,c(as.character(check$Var1[i]),check$Freq[i]))

}

colnames(check_sorted)<-c("Var1","Freq")

levels(data)<-check_sorted$Var1
p13<-DoHeatmap(data, features = top5.markers$gene,draw.lines=F, size=1, angle=90, group.colors=colors30)+ scale_fill_gradientn(guide="colourbar",colors = c("blue", "white", "red"))+ theme(axis.title.x.top =element_text(angle = 90, vjust = 0.5, hjust=1),plot.title = element_text(hjust = 0.5))+ ggtitle(title_name) + theme(plot.title = element_text(hjust = 0.5, vjust=0.5))
print(p13)

}


```
