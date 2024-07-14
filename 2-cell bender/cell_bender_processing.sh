#!/bin/bash
#SBATCH -n 1                    
#SBATCH -t 10:00:00             
#SBATCH --partition=mghpcc-gpu            
#SBATCH --open-mode=append      
#SBATCH -o %j.out               
#SBATCH -e %j.err               
#SBATCH --gres=gpu:NVIDIA_A40:1 
#SBATCH --mem-per-cpu=120G

module load cuda
source /programs/biogrids.shrc
names_list="1305912 1313281 1219699 1305906 1312819 1217423 1313275 1135667 1228994 1135670 1232081 1308756 1219476 1229000 1232087 1308762 1130761 1130764 1217426 1219693"
for i in $names_list;
  do
  mkdir ${i}
  cd ${i}
  cellbender remove-background --cuda --input ../${i}_raw_feature_bc_matrix.h5 --fpr 0 --output cb_${i}_0_fpr.h5
  cd ..
done