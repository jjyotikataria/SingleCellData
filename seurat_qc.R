## Input
## List of filtered_feature_bc_matrix directories
## Violin plot for data qc via Seurat


## Output


files_list <- list.files(path = getwd(), pattern = "_results_GEX", include.dirs = TRUE)
filtered_list <- paste0(files_list,"/","filtered_feature_bc_matrix")

merged_obj <- NULL

for (file in filtered_list) {
    seurat_data <- Read10X(file, strip.suffix=T)
    sample_name <- gsub("_results_GEX/filtered_feature_bc_matrix", "",file)
    seurat_obj <- CreateSeuratObject(counts = seurat_data,
                                     min.cells = 3,
                                     min.features = 200,
                                     project = sample_name)
    assign(sample_name, seurat_obj)

    if(is.null(merged_obj)){
        merged_obj <- seurat_obj
    }
    else{merged_obj <- merge(x=merged_obj, y=seurat_obj)}


}

merged_obj <- PercentageFeatureSet(merged_obj, "^MT-", col.name = "percent_mito")
merged_obj <- PercentageFeatureSet(merged_obj, "^RP[SL]", col.name = "percent_ribo")


suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))

merged_obj@meta.data <- merged_obj@meta.data %>%
     dplyr::rename(nUMI = nCount_RNA,
                   nGene = nFeature_RNA,  
				   percent.mt = percent_mito,
                   percent.rb = percent_ribo)

feats <- c("nGene", "nUMI", "percent.mt", "percent.rb")

print(params$folder_path)

png("seurat_qc.png",width=3000,height=3500,res=400)
VlnPlot(merged_obj, group.by = "orig.ident", features = feats, pt.size = 0.1, ncol = 2) & theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5))
dev.off()
