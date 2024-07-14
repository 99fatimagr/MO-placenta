library(Seurat)
library(decoupleR)
library(pheatmap)
library(tibble)
library(tidyr)
library(dplyr)
sc=readRDS("sc_with_labels_mar_24.rds")

##calcualte scores from progeny 
net <- get_progeny(organism = 'human', top = 500)
mat <- as.matrix(sc@assays$RNA@data)
acts <- run_mlm(mat=mat, net=net, .source='source', .target='target',
                .mor='weight', minsize = 5)
saveRDS(acts,"cells_progeny.rds")

##calculate NDIGSs
net <- readRDS("/lab-share/Gene-Lee-e2/Public/home/fatimagr/obesity/gen3g/child_dev_net_5.5.rds")
mat <- as.matrix(sc@assays$RNA@data)
acts <- run_ulm(mat=mat, net=net, 
        .source='pathway', .target='hgnc_symbol',
        .mor='stat', minsize = 10)
saveRDS(acts,"cells_progeny_child_dev_5.5.rds")

###after processing
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
sc=readRDS("sc_with_labels_mar_24.rds")
DefaultAssay(sc)="RNA"
test=FetchData(sc,vars=c("BMI","cell_type","Organ","Number","Fetal_sex","Batch"))
test$barcode=rownames(test)

toget=merge(test,toget,by="barcode")
head(toget)
saveRDS(toget,"child_dev_gene_scores_with_hypoxia_5.5.rds")
saveRDS(subset(toget,Organ=="decidua"),"child_dev_gene_scores_with_hypoxia_dec.rds")

gene_scores=toget
for_graph=as.data.frame(matrix(nrow=0,ncol=7))

for(i in unique(gene_scores$cell_type)){
  use=subset(gene_scores,cell_type==i)
  a=cor.test(use$hypoxia_score,use$SDQ_3,method="spearman")
  b=cor.test(use$hypoxia_score,use$ASEBA_5,method="spearman")
  c=cor.test(use$hypoxia_score,use$SDQ_5,method="spearman")
  result=c(i,a$estimate,a$p.value,b$estimate,b$p.value,c$estimate,c$p.value)
  for_graph=rbind(for_graph,result)
}

colnames(for_graph)=c("Cell_type" ,   "SDQ_3_cor"  ,"SDQ_3_p", "ASEBA_5_cor","ASEBA_5_p",
                      "SDQ_5_cor","SDQ_5_p")


##mediation analyses
library(mediation)
df=subset(toget,Organ=="decidua"&cell_type=="EVT")

##SDQ 3
fit.totaleffect=lm(SDQ_3_score~BMI,df)
summary(fit.totaleffect)
fit.mediator=lm(hypoxia_score~BMI,df)
summary(fit.mediator)
fit.dv=lm(SDQ_3_score~hypoxia_score+BMI,df)
summary(fit.dv)
results = mediate(fit.mediator, fit.dv, 
                  treat='BMI', mediator='hypoxia_score', boot=T)

###SDQ 5
fit.totaleffect=lm(SDQ_5_score~BMI,df)
summary(fit.totaleffect)
fit.mediator=lm(hypoxia_score~BMI,df)
summary(fit.mediator)
fit.dv=lm(SDQ_5_score~hypoxia_score+BMI,df)
summary(fit.dv)
results = mediate(fit.mediator, fit.dv, 
                  treat='BMI', mediator='hypoxia_score', boot=T)


###ASEBA
fit.totaleffect=lm(ASEBA_5_score~BMI,df)
summary(fit.totaleffect)
fit.mediator=lm(hypoxia_score~BMI,df)
summary(fit.mediator)
fit.dv=lm(ASEBA_5_score~hypoxia_score+BMI,df)
summary(fit.dv)
results = mediate(fit.mediator, fit.dv, 
                  treat='BMI', mediator='hypoxia_score', boot=T)

