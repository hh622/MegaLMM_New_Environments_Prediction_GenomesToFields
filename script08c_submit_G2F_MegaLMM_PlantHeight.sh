#!/bin/bash -l
#SBATCH -D /home/hxhu/P5G2Fv3/scripts
#SBATCH -o /home/hxhu/P5G2Fv3/logs/07_G2F_MegaLMM_%A_%a.out
#SBATCH -e /home/hxhu/P5G2Fv3/logs/07_G2F_MegaLMM_%A_%a.err
#SBATCH -J G2F_MegaLMM
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=20G
#SBATCH --time=24:00:00
#SBATCH -A runciegrp 
#SBATCH -p med2

module load R

echo "###### 07_G2F_MegaLMM_EnvNew_$1_$2_$3_$4_$Pred_Type_$SLURM_ARRAY_TASK_ID.Rout"
echo "###### G2F_$1_Xenv_$2_CV2_Fold0$4_$3_Fold0$SLURM_ARRAY_TASK_ID"

Rscript script07_G2F_MegaLMM_2023-10-02.R $1 $2 $3 $4 $SLURM_ARRAY_TASK_ID  > ../logs/07_G2F_MegaLMM_EnvNew_$1_$2_$3_$4_$Pred_Type_$SLURM_ARRAY_TASK_ID.Rout

# run array jobs
# sbatch --array=1-20 script08c_submit_G2F_MegaLMM.sh Plant_Height_cm Wea CV3NewState 1
# sbatch --array=1-20 script08c_submit_G2F_MegaLMM.sh Plant_Height_cm Wea+TesterGRM CV3NewState 1

