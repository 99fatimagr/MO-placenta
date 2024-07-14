#!/bin/bash
#SBATCH -n 16                  
#SBATCH -t 120:00:00             
#SBATCH -p bch-compute         
#SBATCH --mem-per-cpu=8G       
#SBATCH --open-mode=append      
#SBATCH -o %j.out               
#SBATCH -e %j.err               

source /programs/biogrids.shrc



cellranger count --id 1305912 --sample LIB056126_CRN00238067 --include-introns --transcriptome=/lab-share/Gene-Lee-e2/Public/home/fatimagr/refdata-gex-GRCh38-2020-A --fastqs=/path/to/fastqs

cellranger count --id 1135670 --sample LIB056126_CRN00238068 --include-introns --transcriptome=/lab-share/Gene-Lee-e2/Public/home/fatimagr/refdata-gex-GRCh38-2020-A --fastqs=/path/to/fastqs

cellranger count --id 1313281 --sample LIB056126_CRN00238069 --include-introns --transcriptome=/lab-share/Gene-Lee-e2/Public/home/fatimagr/refdata-gex-GRCh38-2020-A --fastqs=/path/to/fastqs

cellranger count --id 1219699 --sample LIB056126_CRN00238070 --include-introns --transcriptome=/lab-share/Gene-Lee-e2/Public/home/fatimagr/refdata-gex-GRCh38-2020-A --fastqs=/path/to/fastqs

cellranger count --id 1305906 --sample LIB056126_CRN00238071 --include-introns --transcriptome=/lab-share/Gene-Lee-e2/Public/home/fatimagr/refdata-gex-GRCh38-2020-A --fastqs=/path/to/fastqs

cellranger count --id 1312819 --sample LIB056126_CRN00238072 --include-introns --transcriptome=/lab-share/Gene-Lee-e2/Public/home/fatimagr/refdata-gex-GRCh38-2020-A --fastqs=/path/to/fastqs

cellranger count --id 1217423 --sample  LIB056126_CRN00238073 --include-introns --transcriptome=/lab-share/Gene-Lee-e2/Public/home/fatimagr/refdata-gex-GRCh38-2020-A --fastqs=/path/to/fastqs

cellranger count --id 1313275 --sample LIB056126_CRN00238074 --include-introns --transcriptome=/lab-share/Gene-Lee-e2/Public/home/fatimagr/refdata-gex-GRCh38-2020-A --fastqs=/path/to/fastqs
