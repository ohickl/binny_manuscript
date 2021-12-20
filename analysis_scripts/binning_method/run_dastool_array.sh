#!/bin/bash -l
#SBATCH -J binny_eval_dastool_c2_t_urogen_wo_bi #
#SBATCH --mail-type=fail
#SBATCH --mail-user=oskar.hickl@uni.lu
#SBATCH -n 1
#SBATCH -c 7
#SBATCH --time=0-01:30:00
#SBATCH -p batch
#SBATCH --qos=normal
#SBATCH --array 0-9%10
#SBATCH -o "slurm_logs/das_tool/%A-%a-%x-%j.out"

echo "== Starting run at $(date)"
echo "== Job ID: ${SLURM_JOBID}"
echo "== Node list: ${SLURM_NODELIST}"
echo "== Submit dir. : ${SLURM_SUBMIT_DIR}"

shopt -s nullglob
#IN_ARR=(cami_data/c1_t_high/short_read/H_S00*/)
IN_ARR=(cami_data/c2_t_urogen/short_read/201*_sample_*/)
#IN_ARR=(100s/)
shopt -u nullglob

srun das_tool.sh ${IN_ARR[$SLURM_ARRAY_TASK_ID]}
