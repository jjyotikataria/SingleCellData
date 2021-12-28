## Cell cycle scoring from merged subsetted object

library(RCurl)
library(AnnotationHub)

data <-readRDS(paste0(input_files,"/","user_defined_merged_obj.rds"))


if (species == "Human") {
cc_file <- getURL("https://raw.githubusercontent.com/hbc/tinyatlas/master/cell_cycle/Homo_sapiens.csv")

} else if (species=="Mouse") {
cc_file <- getURL("https://raw.githubusercontent.com/hbc/tinyatlas/master/cell_cycle/Mus_musculus.csv")

}


## Cell cycle scoring
cell_cycle_genes <- read.csv(text = cc_file)
# Connect to AnnotationHub
ah <- AnnotationHub()

# Access the Ensembl database for organism
ahDb <- query(ah,
              pattern = c("Homo sapiens", "EnsDb"),
              ignore.case = TRUE)

# Acquire the latest annotation files
id <- ahDb %>%
        mcols() %>%
        rownames() %>%
        tail(n = 1)

# Download the appropriate Ensembldb database
edb <- ah[[id]]

# Extract gene-level information from database
annotations <- genes(edb,
                     return.type = "data.frame")

# Select annotations of interest
annotations <- annotations %>%
        dplyr::select(gene_id, gene_name, seq_name, gene_biotype, description)

# Get gene names for Ensembl IDs for each gene
cell_cycle_markers <- dplyr::left_join(cell_cycle_genes, annotations, by = c("geneID" = "gene_id"))

# Acquire the S phase genes
s_genes <- cell_cycle_markers %>%
        dplyr::filter(phase == "S") %>%
        pull("gene_name")

# Acquire the G2M phase genes
g2m_genes <- cell_cycle_markers %>%
        dplyr::filter(phase == "G2/M") %>%
        pull("gene_name")

## PErform cell scoring
data <- CellCycleScoring(object = data, g2m.features = cc.genes$g2m.genes,
    s.features = cc.genes$s.genes)

## Violin plot distribution
jpeg("cell_cycle_scoring.jpeg",width=2500,height=2000,res=400,units = "px", pointsize = 12)
VlnPlot(data, features = c("S.Score", "G2M.Score"), group.by = "orig.ident",
    ncol = 2, pt.size = 0.1) & theme(axis.title.x=element_blank()) & theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5))

dev.off()

## Doublet Finder

#data.filt <- data
#data.filt = ScaleData(data.filt, vars.to.regress = c("nGene", "percent.mt"), verbose = F)
#data.filt = RunPCA(data.filt, verbose = F, npcs = 20)
#data.filt = RunUMAP(data.filt, dims = 1:10, verbose = F)


#nExp <- round(ncol(data.filt) * 0.04)  # expect 4% doublets
#data.filt <- doubletFinder_v3(data.filt, pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:10)
#DF.name = colnames(data.filt@meta.data)[grepl("DF.classification", colnames(data.filt@meta.data))]



#cowplot::plot_grid(ncol = 2, DimPlot(data.filt, group.by = "orig.ident") + NoAxes(), DimPlot(data.filt, group.by = DF.name) + NoAxes())
#VlnPlot(data.filt, features = "nFeature_RNA", group.by = DF.name, pt.size = 0.1)

saveRDS(data,"merged_final_cell_scoring.rds")
