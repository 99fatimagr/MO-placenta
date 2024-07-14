#!/bin/bash
#SBATCH -n 16                    
#SBATCH -t 120:00:00             
#SBATCH -p bch-compute         
#SBATCH --mem-per-cpu=8G       
#SBATCH --open-mode=append      
#SBATCH -o %j.out               
#SBATCH -e %j.err               

source /programs/biogrids.shrc



cellranger count --id 1217426 --sample LIB055279_CRN00231318 --include-introns --transcriptome=/lab-share/Gene-Lee-e2/Public/home/fatimagr/refdata-gex-GRCh38-2020-A --fastqs=/path/to/fastqs

cellranger count --id 1130764 --sample LIB055279_CRN00231319 --include-introns --transcriptome=/lab-share/Gene-Lee-e2/Public/home/fatimagr/refdata-gex-GRCh38-2020-A --fastqs=/path/to/fastqs

cellranger count --id 1130761 --sample LIB055279_CRN00231320 --include-introns --transcriptome=/lab-share/Gene-Lee-e2/Public/home/fatimagr/refdata-gex-GRCh38-2020-A --fastqs=/path/to/fastqs

cellranger count --id 1219693 --sample LIB055279_CRN00231321 --include-introns --transcriptome=/lab-share/Gene-Lee-e2/Public/home/fatimagr/refdata-gex-GRCh38-2020-A --fastqs=./path/to/fastqs
