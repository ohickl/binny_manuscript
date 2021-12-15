#!/bin/bash -l
#SBATCH -J c2_t_urogen_pool
#SBATCH --mail-type=fail
#SBATCH --mail-user=oskar.hickl@uni.lu
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --time=0-12:00:00
#SBATCH -p batch
#SBATCH --qos=normal
#SBATCH --array 0-9%10 #Num input files for non c2_tg: 39
#SBATCH -o "slurm_logs/mapping/%A-%a-%x-%j.out"

echo "== Starting run at $(date)"
echo "== Job ID: ${SLURM_JOBID}"
echo "== Node list: ${SLURM_NODELIST}"
echo "== Submit dir. : ${SLURM_SUBMIT_DIR}"

shopt -s nullglob
#IN_ARR=(cami_data/cami1_toy_high/H_S00*/)
#IN_ARR=(cami_data/short_read/2017.12.04_18.45.54_sample_*/)
IN_ARR=(cami_data/c2_t_urogen/short_read/201*_sample_*/)
#IN_ARR=(img_data/*/)
#IN_ARR=(${IN_ARR[*]:100:10})
shopt -u nullglob

conda activate mapping

srun cami_mapping_v2.sh ${IN_ARR[$SLURM_ARRAY_TASK_ID]}
# srun cami_mapping_mm2_mb_depth.sh ${IN_ARR[$SLURM_ARRAY_TASK_ID]}
