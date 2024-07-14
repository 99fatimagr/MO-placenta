#!/usr/bin/env Rscript
###note: decidua/ dec refers to maternal side; villous/ vil refers to fetal side
library(Seurat)
library(DESeq2)
library(scuttle)
library(BiocParallel)
register(MulticoreParam(10))

sc=readRDS("sc_with_labels_mar_24.rds")
Idents(sc)="cell_type"
metadata=read.csv("../metadata_pb.csv")
rownames(metadata)=metadata$"Number"

##decidua only
sc_dec=subset(sc,Organ=="decidua")
sc_list=SplitObject(sc_dec, split.by = "ident")

diff.analysis <- lapply(c(1:22), FUN = function(i) {
  sc_use=sc_list[[i]]
  print(i)
  print(names(sc_list)[i])
  Idents(sc_use)="Number"
  DefaultAssay(sc_use)="RNA"
  sc_sce=as.SingleCellExperiment(sc_use)
  summed <- aggregateAcrossCells(sc_sce, ids=sc_sce$Number,use.assay.type=1)
  summed=assays(summed)$counts
  meta_use=subset(metadata,Organ=="decidua")
  meta_use=meta_use[rownames(meta_use) %in% colnames(summed),]
  reorder_idx <- match(colnames(summed),rownames(meta_use))
  meta_use=meta_use[reorder_idx,]
  meta_use$Maternal_age=scale(meta_use$Maternal_age,center=T,scale=T)
  meta_use$Batch=as.character(meta_use$Batch)

  dds <- DESeqDataSetFromMatrix(countData = summed,
                                colData = meta_use,
                                design= ~BMI_cat+Batch+Fetal_sex+Delivery.Mode+Maternal_age)


  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep,]

  dds <- DESeq(dds,parallel=T)
  saveRDS(dds,paste0("DESeq_outs/",names(sc_list)[i],"_dec_DESeq.rds"))
  res <- results(dds, name="BMI_cat_obese_vs_control")
  res_df<-as.data.frame(res)
  return(res)
})

names(diff.analysis) <- names(sc_list)[1:22]
saveRDS(diff.analysis,"dec_pbDEGs.rds")

###villous
sc_vil=subset(sc,Organ=="villous")
sc_list=SplitObject(sc_vil, split.by = "ident")

diff.analysis <- lapply(c(1:16,18), FUN = function(i) {
  sc_use=sc_list[[i]]
  print(i)
  print(names(sc_list)[i])
  Idents(sc_use)="Number"
  DefaultAssay(sc_use)="RNA"
  sc_sce=as.SingleCellExperiment(sc_use)
  summed <- aggregateAcrossCells(sc_sce, ids=sc_sce$Number,use.assay.type=1)
  summed=assays(summed)$counts
  meta_use=subset(metadata,Organ=="villous")
  meta_use=meta_use[rownames(meta_use) %in% colnames(summed),]
  reorder_idx <- match(colnames(summed),rownames(meta_use)) 
  meta_use=meta_use[reorder_idx,]
  meta_use$Maternal_age=scale(meta_use$Maternal_age,center=T,scale=T)
  meta_use$Batch=as.character(meta_use$Batch)
  
  dds <- DESeqDataSetFromMatrix(countData = summed,
                                colData = meta_use,
                                design= ~BMI_cat+Batch+Delivery.Mode+Maternal_age)
  
  
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep,]
  
  dds <- DESeq(dds,parallel=T)
  saveRDS(dds,paste0("DESeq_outs/",names(sc_list)[i],"_vil_DESeq.rds"))
  res <- results(dds, name="BMI_cat_obese_vs_control")
  res_df<-as.data.frame(res)
  return(res)
})

names(diff.analysis) <- names(sc_list)[c(1:16,18)]
saveRDS(diff.analysis,"vil_pbDEGs.rds")


library(decoupleR)
dec=readRDS("dec_pbDEGs.rds")
vil=readRDS("vil_pbDEGs.rds")


net <- get_progeny(organism = 'human', top = 500)
progeny_results=list()
for(i in 1:length(dec)){
  print(i)
  DEGs=dec[[i]]
  contrast_acts <- run_mlm(mat=DEGs[, 'stat', drop=FALSE], net=net,
                           .source='source', .target='target',
                           .mor='weight', minsize = 5)
  progeny_results[[paste0("dec_",names(dec))[i]]]=contrast_acts
}

for(i in 1:length(vil)){
  print(i)
  DEGs=vil[[i]]
  contrast_acts <- run_mlm(mat=DEGs[, 'stat', drop=FALSE], net=net, 
                           .source='source', .target='target',
                           .mor='weight', minsize = 5)
  
  progeny_results[[paste0("vil_",names(vil))[i]]]=contrast_acts
}


progeny_list=lapply(1:length(progeny_results),function(j){
  i=progeny_results[[j]]
  print(i)
  i$FDR=p.adjust(i$p_value,method="fdr")
  i$cell_type=names(progeny_results)[j]
  return(i)
})

progeny_list=do.call(rbind,progeny_list)

progeny_list[["side"]]=sapply(progeny_list$cell_type,function(i){
  strsplit(i,"_")[[1]][1]
})

for(z in c("dec","vil")){
  progeny_sig=subset(progeny_list,FDR<0.05&side==z)
  
  table(progeny_sig$cell_type)
  library(data.table)
  progeny_hm=dcast(progeny_sig,source~cell_type,value.var ="score",fill=0)
  library(pheatmap)
  
  rownames(progeny_hm)=progeny_hm[,1]
  progeny_hm=progeny_hm[,-1]
  
  palette_length = 100
  my_color = colorRampPalette(c("Darkblue", "white","red"))(palette_length)
  
  my_breaks <- c(seq(min(progeny_hm), 0, length.out=ceiling(palette_length/2) + 1),
                 seq(0.05, max(progeny_hm), length.out=floor(palette_length/2)))
  
  # Plot
  pdf(paste0(z,"_progeny_DEG_enrichment.pdf"),height=7)
  pheatmap(progeny_hm, border_color = NA, color=my_color, breaks = my_breaks) 
  dev.off()}

