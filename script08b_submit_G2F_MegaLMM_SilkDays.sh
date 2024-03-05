#!/bin/bash -l
#SBATCH -D /home/hxhu/P5G2Fv3/scripts
#SBATCH -o /home/hxhu/P5G2Fv3/logs_07b/07_G2F_MegaLMM_%A_%a.out
#SBATCH -e /home/hxhu/P5G2Fv3/logs_07b/07_G2F_MegaLMM_%A_%a.err
#SBATCH -J G2F_MegaLMM
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=20G
#SBATCH --time=2-10:00:00
#SBATCH -A runciegrp 
#SBATCH -p med2

module load R

echo "###### G2F_$1_Xenv_$2_CV2_Fold0$4_$3_Fold0$SLURM_ARRAY_TASK_ID"
Rscript script07b_G2F_MegaLMM_SilkDays_2023-10-02.R $1 $2 $3 $4 $SLURM_ARRAY_TASK_ID  > ../logs_07b/07b_G2F_MegaLMM_$1_$2_$3_$4_$Pred_Type_$SLURM_ARRAY_TASK_ID.Rout

# run array jobs
sbatch --array=1-10 script08b_submit_G2F_MegaLMM.sh Silk_DAP_days Wea CV3NewState 1
sbatch --array=1-10 script08b_submit_G2F_MegaLMM.sh Silk_DAP_days Wea+TesterGRM CV3NewState 1


