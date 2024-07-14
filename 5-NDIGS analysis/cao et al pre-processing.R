library(loomR)
lfile <- connect(filename ="GSE156793_S3_gene_count.loom", mode = "r+")

### extract cell attributes (organs to extract placenta and unique sample id)
cells_df = data.frame(Organ=lfile$col.attrs$Organ[], sample=lfile$col.attrs$sample[])

### extract placenta only
placenta_cell_attr = subset(cells_df, Organ=="Placenta")

### extract # of rows for placenta cells
placenta_cell_idx = as.numeric(rownames(placenta_cell_attr))
placenta.matrix <- lfile$matrix[placenta_cell_idx, ]

### add colnames and rownames to the dataframe
placenta_df = data.frame(placenta.matrix)
rownames(placenta_df) = placenta_cell_attr$sample
colnames(placenta_df) = lfile$row.attrs$gene_id[]

saveRDS(placenta_df, file="GSE156793_S3_gene_count.placenta.RDS")


###extract cell types
cell_types = data.frame(Organ=lfile$col.attrs$Organ[], sample=lfile$col.attrs$sample[],
                        cell_type=lfile$col.attrs$Organ_cell_lineage[],
                        fetus_id=lfile$col.attrs$Fetus_id[],
                        cluster_name=lfile$col.attrs$Main_cluster_name[])

placenta_meta=subset(cell_types,Organ=="Placenta")
saveRDS(placenta_meta,"placenta_meta_data.rds")


require(Matrix)
require(Seurat) ###Seurat 4
require(data.table)
require(ggplot2)


metadata=readRDS("placenta_meta_data.rds")
expr_matrix=t(readRDS("GSE156793_S3_gene_count.placenta.RDS"))

min.cells <- 1 
x=expr_matrix
x<- CreateSeuratObject(x,
                       min.cells = min.cells)
x <- NormalizeData(x)

rownames(metadata)=metadata$sample

x=AddMetaData(x,metadata)

Idents(x)="cluster_name"
saveRDS(x,"cao_normalized.rds")