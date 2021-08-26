## input  - user subsetted Seurat objects
## Output - combined PCA plot for all the samples


a<-readRDS("day1_subsetted_obj.rds")
b<-readRDS("day2_subsetted_obj.rds")
c<-readRDS("day3_GEX_subsetted_obj.rds")
d<-readRDS("day4_subsetted_obj.rds")
e<-readRDS("day5_subsetted_obj.rds")
f<-readRDS("day6_subsetted_obj.rds")
g<-readRDS("day7_subsetted_obj.rds")
h<-readRDS("day8_subsetted_obj.rds")


pbmc<-merge(a,c(b,c,d,e,f,g,h), add.cell.ids = c("day1","day2","day3","day4","day5","day6","day7","day8"), project = "10X")

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 3000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
png("merged_pca.png",width=3000,height=2000,res=400)
DimPlot(pbmc, reduction = "pca")
dev.off()
