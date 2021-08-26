## Input : Annotated Seurat objects
## Output : Merged table of cell counts across samples

full_path<-getwd()
files <- list.files(path = full_path, pattern = "*_hpca.annotation.rds$")

data<-c()
tables<-c()

for (i in seq(length(files))){
sam<-gsub("_hpca.annotation.rds","",files[i])
data[[i]]<-readRDS(files[i])
tables[[i]]<-as.data.frame(table(Idents(data[[i]])))
colnames(tables[[i]])<-c("Var1",sam)

if(i>1){
final<-full_join(final,tables[[i]], by="Var1")
} else {
final<-as.data.frame(tables[[i]])
}
}


print(final)
                   Var1   Sam1          Sam2            Sam3
1               T cells   8898         5126            8898
2              NK cells    338          300             338
3             Monocytes     10            8              10
4               B cells     15            5              15
5  Pre-B cells CD34 Neg     25           64              25
6                   BMs      2            7               2
7  Pro-B cells CD34 Pos      2           NA               2
8     Endothelial cells      1           NA               1
9                  GMPs      1            3               1
10      Dendritic cells      1            1               1
11            HSC G-CSF     NA            4              NA
12                 MEPs     NA            1              NA
