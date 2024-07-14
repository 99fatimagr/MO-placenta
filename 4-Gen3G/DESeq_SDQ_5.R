library(data.table)
metadata=read.table("metadata_combined.tsv",sep="\t")
metadata=subset(metadata,status=="included")
ndd=read.csv("progeny/more_ndd.csv")
metadata=merge(metadata,ndd,by.x="ID",by.y="Patient_No_etude",all.x=T)

metadata=subset(metadata,Side=="M")
rownames(metadata)=metadata$ID
RNA_seq_data=as.data.frame(fread("placenta_rnaseq.v1.1.clean.maternal_side.gene_reads.gct",header=T))

rownames(RNA_seq_data)=RNA_seq_data$Name


metadata=subset(metadata,GDM==0 & preeclampsia_hypertension==0 & 
                  Accouchement_Sem_gest_accouch>=37 & smoking==0 & 
                  Accouchement_Poids_zscore>(-1)*1.28&Accouchement_Evacuation==2)
metadata=metadata[!is.na(metadata$score_SDQ_total_5y),]

RNA_seq_data=RNA_seq_data[,colnames(RNA_seq_data) %in% metadata$ID]
reorder_idx <- match(colnames(RNA_seq_data),rownames(metadata)) 
metadata=metadata[reorder_idx,]


metadata$Patient_age_calc_V1=scale(metadata$Patient_age_calc_V1,center=T,scale=T)
metadata$offspring_sex=as.character(metadata$offspring_sex)
metadata$Accouchement_Mode=as.character(metadata$Accouchement_Mode)

metadata$score_SDQ_total_5y=scale(metadata$score_SDQ_total_5y,center=T,scale=T)
metadata$niveau_scol_5y=as.character(metadata$niveau_scol_5y)

library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = RNA_seq_data,
                              colData = metadata,
                              design= ~offspring_sex+Patient_age_calc_V1+Accouchement_Mode+
                              Sample_Kit+niveau_scol_5y+score_SDQ_total_5y)

dds$niveau_scol_5y <- relevel(dds$niveau_scol_5y, ref = "4")

keep <- rowSums(counts(dds)) >= 10
sum(keep)
dds <- dds[keep,]

dds <- DESeq(dds)


saveRDS(dds,"maternal_side_total_SDQ_5.rds")

resultsNames(dds)
res <- results(dds, name="score_SDQ_total_5y")
res_df<-as.data.frame(res)

library("biomaRt")
res_df[["short"]]=sapply(rownames(res_df),function(i){
  strsplit(as.character(i),"[.]")[[1]][1]
})
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <-  res_df$short
gene_IDs <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),
                  values = genes, mart= mart)
head(res_df)
colnames(res_df)=c("baseMean","log2FoldChange","lfcSE","stat","pvalue",
                   "padj","ensembl_gene_id")
res_gene=merge(res_df,gene_IDs,by="ensembl_gene_id")
res_gene_sm=subset(res_gene,hgnc_symbol!='')
head(res_gene)

write.csv(res_gene,"maternal_side_SDQ_5_total.csv")
