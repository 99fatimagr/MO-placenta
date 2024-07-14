require(Matrix)
require(Seurat) 
require(data.table)
require(ggplot2)

##batch 1
Dirs <- c("1130761","1130764","1217426","1219693")
prefixes=c("61","64","26","93")
##for each apply the function below

##batch 2
Dirs <- c("1135670","1217423","1219699","1305906",
          "1305912","1312819","1313275","1313281")
prefixes=c("70","23","99","06","12","19","75","81")

###batch 3
Dirs <- c("1135667","1219476","1228994","1229000","1232081",
          "1232087","1308756","1308762")
prefixes=c("67","76","94","9000","2081","2087","56","62")


seurat_list=sapply(1:length(Dirs),function(i){
  #make Seurat object
  name=Dirs[i]
  #counts=as.data.frame(fread(paste0('filtered/',name,"_gene_expression.csv")))
  counts=as.data.frame(fread(paste0('filtered_expr_matrices/',name,"_gene_expression.csv")))
  rownames(counts)=counts$V1
  counts=counts[,-1]
  counts=t(counts)
  print(name)
  print(dim(counts))
  meta.data=as.data.frame(fread(paste0("filtered_expr_matrices/",name,"_metadata.csv")))
  rownames(meta.data)=meta.data$V1
  meta.data=meta.data[,-1]
  sc=CreateSeuratObject(counts=counts,meta.data=meta.data)
  #rename barcodes
  sc=RenameCells(sc, add.cell.id=prefixes[i])
  return(sc)
})

##then save
##batch 1
sc_combined=merge(seurat_list[[1]],y=c(seurat_list[[2]],seurat_list[[3]],seurat_list[[4]]))

saveRDS(sc_combined,"batch_1_seurat_object.rds")
##batch 2
sc_combined=merge(seurat_list[[1]],y=c(seurat_list[[2]],seurat_list[[3]],seurat_list[[4]],
                                       seurat_list[[5]],seurat_list[[6]],seurat_list[[7]],
                                       seurat_list[[8]]))
saveRDS(sc_combined,"batch_2_seurat_object.rds")

##batch 3
sc_combined=merge(seurat_list[[1]],y=c(seurat_list[[2]],seurat_list[[3]],seurat_list[[4]],
                                       seurat_list[[5]],seurat_list[[6]],seurat_list[[7]],
                                       seurat_list[[8]]))

saveRDS(sc_combined,"batch_3_seurat_object.rds")

###seurat integration script
pilot_expression_matrix=readRDS("batch_1_seurat_object.rds")
batch_2_expression_matrix=readRDS("batch_2_seurat_object.rds")
batch_3_expression_matrix=readRDS("batch_3_seurat_object.rds")
seurat.list=as.list(c(pilot_expression_matrix,batch_2_expression_matrix,batch_3_expression_matrix))

seurat.list <- lapply(X = seurat.list, FUN = function(x) {
  print(x)
  x=subset(x,nFeature_RNA > 200 & nCount_RNA>500)
  print(x)
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = seurat.list,nfeatures = 2000)
anchors <- FindIntegrationAnchors(object.list = 
                                           seurat.list, anchor.features = features)
sample_tree=matrix(c(-2, 1, -3, -1), ncol = 2)
seurat.combined <- IntegrateData(anchorset = anchors,sample.tree=sample_tree,
                                 normalization.method="LogNormalize")
saveRDS(seurat.combined,"seurat_combined.rds")


seurat.combined=readRDS("seurat_combined.rds")
DefaultAssay(seurat.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
seurat.combined <- ScaleData(seurat.combined, verbose = FALSE)
seurat.combined <- RunPCA(seurat.combined, npcs = 30, verbose = FALSE)
Elbow_plot=ElbowPlot(seurat.combined)
ggsave("elbow.png",Elbow_plot)
seurat.combined <- RunUMAP(seurat.combined, reduction = "pca", dims = 1:30)
seurat.combined <- FindNeighbors(seurat.combined, reduction = "pca", dims = 1:30)

seurat.combined[["orig.ident"]]<-factor(unlist(lapply(rownames(seurat.combined@meta.data), function(j){
  strsplit(as.character(j),"_")[[1]][1]
})))

###add metadata info on batch, obesity status, maternal age, fetal sex, etc.
metadata=read.csv("../metadata.txt",sep="\t")
meta_from_file=function(num){
  name=colnames(metadata)[num]
  return(sapply(colnames(seurat.combined),function(i){
    prefix=strsplit(i,"_")[[1]][1]
    index=which(metadata$Name==as.numeric(prefix))
    return(metadata[index,num])
  }))
}

for(i in 2:13){
  data=meta_from_file(i)
  col_name=colnames(metadata)[i]
  seurat.combined=AddMetaData(object=seurat.combined,data,
                              col.name=col_name)
  print(head(seurat.combined@meta.data))
}

saveRDS(seurat.combined,"post_UMAP.rds")


# For performing differential expression after integration, we switch back to the original
# data, and use an SingleCellExperiment object to interface with scDblFinder package 
DefaultAssay(seurat.combined) <- "RNA"
sce <- as.SingleCellExperiment(seurat.combined)
saveRDS(sce,"RNA_sce.rds")

library(scDblFinder)
library(BiocParallel)
sce=readRDS("RNA_sce.rds")
sce <- scDblFinder(sce,samples="orig.ident", BPPARAM=MulticoreParam(5))
saveRDS(sce,"scDblFinder_scores.rds")

###back to seurat
library(Seurat)
seurat.combined=readRDS("post_UMAP.rds")
dbl=readRDS("scDblFinder_scores.rds")
seurat.combined$scDblFinder.score <- dbl$scDblFinder.score
seurat.combined$scDblFinder.class <- dbl$scDblFinder.class
seurat.combined$cxds_score <- dbl$scDblFinder.cxds_score

seurat.combined <- FindClusters(seurat.combined, algorithm=4, res=1.5, method='igraph')
library(SingleR)
ref=readRDS("path/to/first_trim/reference.rds")
ref.sce=as.SingleCellExperiment(ref)
DefaultAssay(seurat.combined)="RNA"
pred.ref <- SingleR(test=as.SingleCellExperiment(seurat.combined),
                      ref=ref, labels=ref.sce$annotated_cluster,
                      de.method="wilcox")
clust=pred.ref$labels
names(clust)=rownames(pred.ref)
seurat.combined=AddMetaData(seurat.combined,clust,col.name="teich_cluster_singleR")

#remove all doublets
seurat.combined=subset(seurat.combined,scDblFinder.class=="singlet")
#remove doublet dominated clusters
seurat.combined=subset(seurat.combined,integrated_snn_res.1.5!=35)
seurat.combined=subset(seurat.combined,integrated_snn_res.1.5!=37)
saveRDS(seurat.combined,"post_dblf.rds")

