library(Seurat)
library(decoupleR)
sc=readRDS("cao_normalized.rds")
library(pheatmap)
library(tibble)
library(tidyr)
library(dplyr)

library("biomaRt")
mat <- as.data.frame(sc@assays$RNA@data)
mat[["short"]]=sapply(rownames(mat),function(i){
  return(strsplit(as.character(i),"[.]")[[1]][1])
})
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <-  mat$short
gene_IDs <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),
                  values = genes, mart= mart)

mat=merge(mat,gene_IDs,by.y="ensembl_gene_id",by.x="short")
excl=table(mat$hgnc_symbol)[table(mat$hgnc_symbol)!=1]
mat=mat[!mat$hgnc_symbol %in% names(excl),]
rownames(mat)=mat$hgnc_symbol
remove=c(1,dim(mat)[2])
mat=mat[,-c(remove)]
mat=as.matrix(mat)


net <- readRDS("/lab-share/Gene-Lee-e2/Public/home/fatimagr/obesity/gen3g/child_dev_net_5.5.rds")
acts <- run_ulm(mat=mat, net=net,
        .source='pathway', .target='hgnc_symbol',
        .mor='stat', minsize = 10)
saveRDS(acts,"cells_progeny_child_dev_5.5.rds")


net <- get_progeny(organism = 'human', top = 500)
acts <- run_mlm(mat=mat, net=net, .source='source', .target='target',
                .mor='weight', minsize = 5)
saveRDS(acts,"cells_progeny.rds")

acts=readRDS("cells_progeny.rds")
hypoxia=subset(acts,source=="Hypoxia")

acts=readRDS("cells_progeny_child_dev_5.5.rds")
SDQ=subset(acts,source=="SDQ_3")
ASEBA=subset(acts,source=="ASEBA")
SDQ_5=subset(acts,source=="SDQ_5")

toget=merge(hypoxia,SDQ,by="condition")
toget=merge(toget,ASEBA,by="condition")
toget=merge(toget,SDQ_5,by="condition")

toget=toget[,c(1,4,8,12,16)]
colnames(toget)=c("barcode","hypoxia_score","SDQ_3_score","ASEBA_5_score","SDQ_5_score")

library(Seurat)
sc=readRDS("cao_normalized.rds")

test=FetchData(sc,vars=c("cluster_name","fetus_id"))
test$barcode=rownames(test)

toget=merge(test,toget,by="barcode")

