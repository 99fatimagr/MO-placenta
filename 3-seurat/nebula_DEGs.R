#!/usr/bin/env Rscript
###note: decidua/ dec refers to maternal side; villous/ vil refers to fetal side
library(Seurat)
library(DESeq2)
library(scuttle)
library(BiocParallel)
register(MulticoreParam(10))

sc=readRDS("sc_with_labels_mar_24.rds")
Idents(sc)="cell_type"

###calculate DEGs
sc_dec=subset(sc,Organ=="decidua")
sc_list=SplitObject(sc_dec, split.by = "ident")

diff.analysis <- lapply(c(1:22), FUN = function(i) {
  sc_use=sc_list[[i]]
  print(i)
  print(names(sc_list)[i])
  counts=sc_use@assays$RNA@counts
  meta <- sc_use@meta.data
  subjects <- sc_use@meta.data$Patient
  group_levels <- c("control","obese")
  meta$BMI_cat=factor(meta$BMI_cat, levels = group_levels)
  meta$offset=colSums(counts)
  meta$Batch=as.character(meta$Batch)
  delivery_levels <- c("C","VS","VI")
  meta$Delivery.Mode=factor(meta$Delivery.Mode, levels = delivery_levels)
  df <- model.matrix(~BMI_cat+Batch+Delivery.Mode+Fetal_sex+Maternal_age, data=meta)
  re = nebula(counts,subjects,pred=df,model='NBLMM',reml=1,method='HL',ncore=6)
  re$summary$padj <- p.adjust(re$summary$p_BMI_catobese, method = "fdr")
  return(re)
})

names(diff.analysis) <- names(sc_list)
saveRDS(diff.analysis,"dec_nebula.rds")

###villous
sc_vil=subset(sc,Organ=="villous")
sc_list=SplitObject(sc_vil, split.by = "ident")
diff.analysis <- lapply(c(1:16,18), FUN = function(i) {
  sc_use=sc_list[[i]]
  print(i)
  print(names(sc_list)[i])
  counts=sc_use@assays$RNA@counts
  meta <- sc_use@meta.data
  subjects <- sc_use@meta.data$Patient
  group_levels <- c("control","obese")
  meta$BMI_cat=factor(meta$BMI_cat, levels = group_levels)
  meta$offset=colSums(counts)
  meta$Batch=as.character(meta$Batch)
  delivery_levels <- c("C","VS","VI")
  meta$Delivery.Mode=factor(meta$Delivery.Mode, levels = delivery_levels)
  #fetal sex is collinear with delivery mode in villous only
  df <- model.matrix(~BMI_cat+Batch+Delivery.Mode+Maternal_age, data=meta)
  re = nebula(counts,subjects,pred=df,model='NBLMM',reml=1,method='HL',ncore=6)
  re$summary$padj <- p.adjust(re$summary$p_BMI_catobese, method = "fdr")
  return(re)
})
names(diff.analysis) <- names(sc_list)[c(1:16,18)]
saveRDS(diff.analysis,"vil_nebula.rds")

###PROGENy enrichment
##decidua only
library(decoupleR)
dec=readRDS("dec_nebula.rds")
vil=readRDS("vil_nebula.rds")

net <- get_progeny(organism = 'human', top = 500)
progeny_results=list()

for(i in 1:length(dec)){
  print(i)
  use=dec[[i]]$summary
  use$rank <- sign(use$logFC_BMI_catobese)*-log10(use$p_BMI_catobese+1e-300)

  gsea.rank <- as.data.frame(use[,c("gene","rank")])
  rownames(gsea.rank)=gsea.rank[,1]
  gsea.rank=gsea.rank[dec[[i]]$convergence>(-1)*(20),]
  gsea.rank[,1]=0
  contrast_acts <- run_mlm(mat=gsea.rank[, 'rank', drop=FALSE], net=net, .source='source', .target='target',
                           .mor='weight', minsize = 5)
  progeny_results[[paste0("dec_",names(dec))[i]]]=contrast_acts
}

for(i in 1:length(vil)){
  print(i)
  use=vil[[i]]$summary
  use$rank <- sign(use$logFC_BMI_catobese)*-log10(use$p_BMI_catobese+1e-300)
  
  gsea.rank <- as.data.frame(use[,c("gene","rank")])
  rownames(gsea.rank)=gsea.rank[,1]
  gsea.rank=gsea.rank[vil[[i]]$convergence>(-1)*(20),]
  gsea.rank[,1]=0
  contrast_acts <- run_mlm(mat=gsea.rank[, 'rank', drop=FALSE], net=net, .source='source', .target='target',
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
  pdf(paste0(z,"_progeny_nebula_DEG_enrichment.pdf"),height=7)
  pheatmap(progeny_hm, color=my_color, breaks = my_breaks) 
  dev.off()}
