#!/bin/bash -l
#SBATCH -J binny_eval_metawrap_c2_t_urogen #
#SBATCH --mail-type=fail
#SBATCH --mail-user=oskar.hickl@uni.lu
#SBATCH -n 1
#SBATCH -c 4
#SBUTCH --mem-per-cpu=40GB
#SBATCH --time=1-00:00:00
#SBATCH -p bigmem
#SBATCH --array 0-9%10
#SBATCH -o "slurm_logs/metawrap_refiner/%A-%a-%x-%j.out"

echo "== Starting run at $(date)"
echo "== Job ID: ${SLURM_JOBID}"
echo "== Node list: ${SLURM_NODELIST}"
echo "== Submit dir. : ${SLURM_SUBMIT_DIR}"

shopt -s nullglob
#IN_ARR=(cami_data/c1_t_high/short_read/H_S00*/)
IN_ARR=(cami_data/c2_t_urogen/short_read/201*_sample_*/)
#IN_ARR=(100s/)
shopt -u nullglob
## Because sample 0 is already processed.
#IN_ARR=("${IN_ARR[@]:1}")
## Run samples that failed becaue of time limit.
#IN_ARR2=("${IN_ARR[1]}" "${IN_ARR[2]}" "${IN_ARR[6]}" "${IN_ARR[7]}")

srun metawrap_refine.sh ${IN_ARR[$SLURM_ARRAY_TASK_ID]}
#srun metawrap_refine_cc_only.sh ${IN_ARR[$SLURM_ARRAY_TASK_ID]}
