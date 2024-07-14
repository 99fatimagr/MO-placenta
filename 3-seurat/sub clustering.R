library(Seurat)
library(ggplot2)
seurat.combined=readRDS("post_dblf.rds")
troph=subset(seurat.combined,integrated_snn_res.1.5 %in% c(2,3,4,6,8,9,10,11,12,13,14,15,16,
                                                           17,21,26,27,36))
immune=subset(seurat.combined,integrated_snn_res.1.5 %in% c(1,18,19,20,22,28,30,31))
endo_stromal=subset(seurat.combined,integrated_snn_res.1.5 %in% c(5,7,23,24,25,29,32,33,34,38))

min.cells=0

broad_groups=list(troph,immune,endo_stromal)
names=c("troph","immune","endo_stromal")
for(i in 1:3){
  pilot_expression_matrix=readRDS("batch_1_seurat_object.rds")
  batch_2_expression_matrix=readRDS("batch_2_seurat_object.rds")
  batch_3_expression_matrix=readRDS("batch_3_seurat_object.rds")
  seurat.list=as.list(c(pilot_expression_matrix,batch_2_expression_matrix,batch_3_expression_matrix))
  
  seurat.list <- lapply(X = seurat.list, FUN = function(x) {
    print(x)
    x=x[,colnames(x) %in% rownames(broad_groups[[i]]@meta.data)]
    print(x)
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  })
  
  # select features that are repeatedly variable across datasets for integration
  features <- SelectIntegrationFeatures(object.list = seurat.list,nfeatures=2000)
  anchors <- FindIntegrationAnchors(object.list = 
                                      seurat.list, anchor.features = features)
  sample_tree=matrix(c(-2, 1, -3, -1), ncol = 2)
  seurat.combined <- IntegrateData(anchorset = anchors,sample.tree=sample_tree,
                                   normalization.method="LogNormalize")
  
  DefaultAssay(seurat.combined) <- "integrated"
  seurat.combined=AddMetaData(seurat.combined,broad_groups[[i]]@meta.data)
  seurat.combined <- ScaleData(seurat.combined, verbose = FALSE)
  seurat.combined <- RunPCA(seurat.combined, npcs = 30, verbose = FALSE)
  seurat.combined <- RunUMAP(seurat.combined, reduction = "pca", dims = 1:15)
  seurat.combined <- FindNeighbors(seurat.combined, reduction = "pca", dims = 1:15)
  a=DimPlot(seurat.combined,group.by="teich_cluster_singleR")
  ggsave(paste0(names[i],"_pcs_15_UMAP.pdf"),a)
  saveRDS(seurat.combined,paste0(names[i],"_pcs_15_integrated.rds"))}
  
  
troph=readRDS("troph_pcs_15_integrated.rds")
immune=readRDS("immune_pcs_15_integrated.rds")
endo_stromal=readRDS("endo_stromal_pcs_15_integrated.rds")

broad_groups=list(troph,immune,endo_stromal)
names=c("troph","immune","endo_stromal")
for(i in 1:3){
  use=broad_groups[[i]]
  
    use <- FindClusters(use, algorithm=4, res=0.8, method='igraph')
    markers=FindAllMarkers(use,logfc.threshold=0.5)
    write.csv(markers,paste0("markers_res_0.8_",names[i],".csv"))
    a=DimPlot(use,label=T)
    ggsave(paste0(names[i],"_0.8.pdf"),a)
  saveRDS(use,paste0(names[i],"_clustered.rds"))}

  
library(Seurat)
troph=readRDS("sub_group_processing/troph_clustered.rds")
Idents(troph)="integrated_snn_res.0.8"
new_idents=c("SCT_FLT1","SCT","SCT","SCT","SCT","SCT_ECM","SCT","SCT_FLT1","npVCT",
             "SCT","SCT_BACE2","SCT",
             "SCT_PTCHD4","SCT","pVCT","EVT","fusing_VCT")
names(new_idents)=levels(Idents(troph))
troph=RenameIdents(troph,new_idents)
troph@meta.data[['cell_type']]=Idents(troph)
troph$large_group="trophoblast"


es=readRDS("sub_group_processing/endo_stromal_clustered.rds")
Idents(es)="integrated_snn_res.0.8"
new_idents=c("Endo","FB1","FB2","Endo","Doublet","Endo","FB1","Endo","FB1","FB2","Endo","FB2",
             "EGC","Endo_M","dS","Doublet","LED")
names(new_idents)=levels(Idents(es))
es=RenameIdents(es,new_idents)
es@meta.data[['cell_type']]=Idents(es)
es$large_group="endothelial_stromal"

immune=readRDS("sub_group_processing/immune_clustered.rds")
Idents(immune)="integrated_snn_res.0.8"
new_idents=c("FM","T","FM","FM",
             "NK_CD16+","dM2","dM1","Doublet_8",
             "BM","dM1","FM","dNK","T",
             "B_cells","Doublet_15")
names(new_idents)=levels(Idents(immune))
immune=RenameIdents(immune,new_idents)
immune@meta.data[['cell_type']]=Idents(immune)


immune$large_group="immune"

data=lapply(list(troph,es,immune),function(i){
  a=FetchData(i,c("cell_type","large_group"))
  return(a)
})

data=do.call(rbind,data)
sc=readRDS("whole_object_processing/post_dblf.rds")
sc=AddMetaData(sc,data)
sc=subset(sc,cell_type != "Doublet")
sc=subset(sc,cell_type!="Doublet_8")
sc=subset(sc,cell_type!="Doublet_15")

saveRDS(sc,"sc_with_labels_mar_24.rds")