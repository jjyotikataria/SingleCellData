## Input: Annotated Seurat object
## Output: Table of cell types and their corresponding numbers 

library(Seurat)
library(dplyr)

files <- list.files(path = full_path, pattern = "*_mid.annotation.rds$")

data<-c()
tables<-c()
for (i in seq(length(files))){
sam<-gsub("_mid.annotation.rds","",files[i])
data[[i]]<-readRDS(files[i])
tables[[i]]<-as.data.frame(table(Idents(data[[i]])))
colnames(tables[[i]])<-c("Var1",sam)
if(i>1){
final<-full_join(final,tables[[i]], by="Var1")
}
else{
final<-as.data.frame(tables[[i]])
}
}
