## Input
         Cell_Type         sample1        sample2
1       NK T cells          127              95
2      CD8 T cells         3222            5342
3      CD4 T cells         2062            3779
4        Monocytes           20               5
5         NK cells           34              40

## Output
     Cell_Type DMSO_LIB_GEX PEPTIDE_LIB_GEX
1       NK T cells    127(2.3%)       95(1.02%)
2      CD8 T cells 3222(58.38%)    5342(57.48%)
3      CD4 T cells 2062(37.36%)    3779(40.67%)
4        Monocytes    20(0.36%)        5(0.05%)
5         NK cells    34(0.62%)       40(0.43%)
6           Total         5465            9261


final<-read.csv("input.csv")
final<-final[-1]
change<-names(final[-1])
final[change]<-sapply(final[change],as.numeric)
final2<-final[-1] %>% mutate_all(~paste0(.x, '(', round(.x * 100 / sum(.x),2), '%)'))
final3<-cbind(final[1],final2)

join<-tail(rbind(final, data.frame(Cell_Type = "Total", t(colSums(final[, -1])))),1)
final4<-rbind(final3,join)
