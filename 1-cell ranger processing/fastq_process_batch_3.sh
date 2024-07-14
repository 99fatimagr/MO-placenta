#!/bin/bash
#SBATCH -n 16                  
#SBATCH -t 120:00:00             
#SBATCH -p bch-compute         
#SBATCH --mem-per-cpu=8G       
#SBATCH --open-mode=append      
#SBATCH -o %j.out               
#SBATCH -e %j.err               

source /programs/biogrids.shrc



cellranger count --id 1229000 --sample LIB056564_CRN00242294 --include-introns --transcriptome=/lab-share/Gene-Lee-e2/Public/home/fatimagr/refdata-gex-GRCh38-2020-A --fastqs=/path/to/fastqs

cellranger count --id 1308762 --sample LIB056564_CRN00242295 --include-introns --transcriptome=/lab-share/Gene-Lee-e2/Public/home/fatimagr/refdata-gex-GRCh38-2020-A --fastqs=/path/to/fastqs

cellranger count --id 1219476 --sample LIB056564_CRN00242296 --include-introns --transcriptome=/lab-share/Gene-Lee-e2/Public/home/fatimagr/refdata-gex-GRCh38-2020-A --fastqs=/path/to/fastqs

cellranger count --id 1232087 --sample LIB056564_CRN00242297 --include-introns --transcriptome=/lab-share/Gene-Lee-e2/Public/home/fatimagr/refdata-gex-GRCh38-2020-A --fastqs=/path/to/fastqs

cellranger count --id 1308756 --sample LIB056564_CRN00242298 --include-introns --transcriptome=/lab-share/Gene-Lee-e2/Public/home/fatimagr/refdata-gex-GRCh38-2020-A --fastqs=/path/to/fastqs

 cellranger count --id 1135667 --sample LIB056564_CRN00242299 --include-introns --transcriptome=/lab-share/Gene-Lee-e2/Public/home/fatimagr/refdata-gex-GRCh38-2020-A --fastqs=/path/to/fastqs

cellranger count --id 1232081 --sample  LIB056564_CRN00242300 --include-introns --transcriptome=/lab-share/Gene-Lee-e2/Public/home/fatimagr/refdata-gex-GRCh38-2020-A --fastqs=/path/to/fastqs

cellranger count --id 1228994 --sample LIB056564_CRN00242301 --include-introns --transcriptome=/lab-share/Gene-Lee-e2/Public/home/fatimagr/refdata-gex-GRCh38-2020-A --fastqs=/path/to/fastqs
