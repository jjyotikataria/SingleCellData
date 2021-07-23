## Reading the Seurat object




## Saving Seurat object into a dataframe

library(data.table)
data_to_write_out <- as.data.frame(as.matrix(seuratObject@scale.data))
fwrite(x = data_to_write_out, file = "outfile.csv")
 fwrite(x = data_to_write_out, row.names = TRUE, file = "outfile.csv")
