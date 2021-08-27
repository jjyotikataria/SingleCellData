## Input: User subsetted filtered object
## Output: scSorter annotated Seurat object


time_zero <- readRDS("time_zero_filt_var.rds")  
time_three <- readRDS("time_three_filt_var.rds")  
time_seven_alk <- readRDS("time_seven_alk_filt_var.rds")  
time_seven_dm <- readRDS("time_seven_dm_filt_var.rds")

#----------------
all.genes.time_zero <- rownames(time_zero)
time_zero <- ScaleData(time_zero, features = all.genes.time_zero)


all.genes.time_three <- rownames(time_three)
time_three <- ScaleData(time_three, features = all.genes.time_three)

all.genes.time_seven_dm <- rownames(time_seven_dm)
time_seven_dm <- ScaleData(time_seven_dm, features = all.genes.time_seven_dm)


all.genes.time_seven_alk <- rownames(time_seven_alk)
time_seven_alk <- ScaleData(time_seven_alk, features = all.genes.time_seven_alk)

#---

time_zero <- RunPCA(time_zero, features = VariableFeatures(object = time_zero))
time_three <- RunPCA(time_three, features = VariableFeatures(object = time_three))
time_seven_dm <- RunPCA(time_seven_dm, features = VariableFeatures(object = time_seven_dm))
time_seven_alk <- RunPCA(time_seven_alk, features = VariableFeatures(object = time_seven_alk))

#---

time_zero <- RunUMAP(time_zero, reduction = "pca", dims = 1:20)
time_three <- RunUMAP(time_three, reduction = "pca", dims = 1:20)
time_seven_dm <- RunUMAP(time_seven_dm, reduction = "pca", dims = 1:20)
time_seven_alk <- RunUMAP(time_seven_alk, reduction = "pca", dims = 1:20)

--

anno=read.csv(file="markers_lungs.csv",header=TRUE)

data_topgenes_time_zero <- head(VariableFeatures(time_zero), 3000)
data_topgenes_time_three <- head(VariableFeatures(time_three), 3000)
data_topgenes_time_seven_dm <- head(VariableFeatures(time_seven_dm), 3000)
data_topgenes_time_seven_alk <- head(VariableFeatures(time_seven_alk), 3000)

count_time_zero = GetAssayData(time_zero)
count_time_three = GetAssayData(time_three)
count_time_seven_dm = GetAssayData(time_seven_dm)
count_time_seven_alk = GetAssayData(time_seven_alk)

picked_genes1 = unique(c(anno$Marker, data_topgenes_time_zero))
picked_genes2 = unique(c(anno$Marker, data_topgenes_time_three))
picked_genes3 = unique(c(anno$Marker, data_topgenes_time_seven_dm))
picked_genes4 = unique(c(anno$Marker, data_topgenes_time_seven_alk))


count_time_zero = count_time_zero[rownames(count_time_zero) %in% picked_genes1, ]
data_rts1 <- scSorter(count_time_zero, anno)

count_time_three = count_time_three[rownames(count_time_three) %in% picked_genes2, ]
data_rts2 <- scSorter(count_time_three, anno)

count_time_seven_dm = count_time_seven_dm[rownames(count_time_seven_dm) %in% picked_genes3, ]
data_rts3 <- scSorter(count_time_seven_dm, anno)

count_time_seven_alk = count_time_seven_alk[rownames(count_time_seven_alk) %in% picked_genes4, ]
data_rts4 <- scSorter(count_time_seven_alk, anno)

data_cells_name1 <- data_rts1$Pred_Type
data_cells_name2 <- data_rts2$Pred_Type
data_cells_name3 <- data_rts3$Pred_Type
data_cells_name4 <- data_rts4$Pred_Type


time_zero <- AddMetaData(object = time_zero,metadata = data_cells_name1,col.name = 'type')
time_three <- AddMetaData(object = time_three,metadata = data_cells_name2,col.name = 'type')
time_seven_dm <- AddMetaData(object = time_seven_dm,metadata = data_cells_name3,col.name = 'type')
time_seven_alk <- AddMetaData(object = time_seven_alk,metadata = data_cells_name4,col.name = 'type')


Idents(time_zero) <- time_zero@meta.data$type
Idents(time_three) <- time_three@meta.data$type
Idents(time_seven_dm) <- time_seven_dm@meta.data$type
Idents(time_seven_alk) <- time_seven_alk@meta.data$type
