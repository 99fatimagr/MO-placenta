library(data.table)
metadata=read.table("../metadata_combined.tsv",sep="\t")
metadata=subset(metadata,status=="included")

metadata=subset(metadata,Side=="F")
rownames(metadata)=metadata$ID
RNA_seq_data=as.data.frame(fread("../placenta_rnaseq.v1.1.clean.fetal_side.gene_reads.gct",header=T))

rownames(RNA_seq_data)=RNA_seq_data$Name

metadata[["BMI_cat"]]=sapply(metadata$BMI_V1,function(i){
  result=0
  if(i<18.5){
    result="under"
  }
  else if(i<25){
    result="normal"
  }
  else if(i<30){
    result="over"
  }
  else{result="obese"}
  return(result)
})


metadata=subset(metadata,GDM==0 & preeclampsia_hypertension==0 & 
                  Accouchement_Sem_gest_accouch>=37 & smoking==0 &
                  Accouchement_Poids_zscore>(-1)*1.28&Accouchement_Evacuation==2)

RNA_seq_data=RNA_seq_data[,colnames(RNA_seq_data) %in% metadata$ID]
reorder_idx <- match(colnames(RNA_seq_data),rownames(metadata)) 
metadata=metadata[reorder_idx,]

metadata$Patient_age_calc_V1=scale(metadata$Patient_age_calc_V1,center=T,scale=T)
metadata$offspring_sex=as.character(metadata$offspring_sex)
metadata$Accouchement_Mode=as.character(metadata$Accouchement_Mode)

###use DESEq2 to get variance stabilized counts
library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = RNA_seq_data,
                              colData = metadata,
                              design= ~offspring_sex+Patient_age_calc_V1+Sample_Kit+BMI_cat)


keep <- rowSums(counts(dds)) >= 10
sum(keep)
dds <- dds[keep,]
vsd <- vst(dds, blind=TRUE)
counts=as.data.frame(assay(vsd))

library("biomaRt")
counts[["short"]]=sapply(rownames(counts),function(i){
  strsplit(as.character(i),"[.]")[[1]][1]
})
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <-  counts$short
gene_IDs <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),
                  values = genes, mart= mart)

colnames(counts)[dim(counts)[2]]="ensembl_gene_id"
res_gene=merge(counts,gene_IDs,by="ensembl_gene_id")
res_gene_sm=subset(res_gene,hgnc_symbol!='')
head(res_gene)

#remove duplicate genes
library(dplyr)
counts=res_gene %>%
  distinct(hgnc_symbol, .keep_all = TRUE)
rownames(counts)=counts$hgnc_symbol
counts=counts[,-c(1,ncol(counts))]

library(decoupleR)
net <- get_progeny(organism = 'human', top = 500)
net


sample_acts <- run_mlm(mat=counts, net=net, .source='source', .target='target',
                       .mor='weight', minsize = 5)
with_meta=merge(metadata,sample_acts,by.x="ID",by.y="condition")
ndd=read.csv("ndd_data.csv")
with_meta=merge(with_meta,ndd,by.x="ID",by.y="Patient_No_etude",all.x=T)

