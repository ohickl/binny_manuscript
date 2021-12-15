#!/bin/bash -l
#SBATCH -J binning_eval_c2_t_urogen_binny
#SBATCH --mail-type=fail
#SBATCH --mail-user=oskar.hickl@uni.lu
#SBATCH -n 1
#SBATCH -c 128
#SBUTCH --mem-per-cpu=40GB
#SBATCH --time=0-12:00:00
#SBATCH -p batch
# #SBATCH --qos=normal
#SBATCH --array 0-9%10
#SBATCH -o "slurm_logs/checkm/%A-%a-%x-%j.out"

echo "== Starting run at $(date)"
echo "== Job ID: ${SLURM_JOBID}"
echo "== Node lis: ${SLURM_NODELIST}"
echo "== Submit dir. : ${SLURM_SUBMIT_DIR}"

shopt -s nullglob
#DATA_SET='mb_genomes'  # CAMI1_toy_high, CAMI2_toy_gut mb_genomes

# IN_ARR=(${1}_out/c2_t_urogen/201*_sample_*/)

#IN_ARR=(binny_out/c2_t_airways/201*_sample_*/)
#IN_ARR=(binny_out/c1_t_high/H_S00?/)

IN_ARR=(binny_out/img_genomes/*_exp_binny/)
IN_ARR=(${IN_ARR[*]:0:10})

# For vamb leftover
#IN_ARR=( "$@" )

#IN_ARR=(binny_out/img_genomes/KBSS*/)

# IN_ARR=(binny_out/img_genomes/June2015DPH_4_15_FD_exp_binny/)

#DATA_SET='100s'
#IN_ARR=(concoct_out/100s concoct_out/100s concoct_out/100s)
#DATA_SET='CAMI1_toy_high'
#DATA_SET='CAMI1_toy_high_pool'
#IN_ARR=(concoct_out/c1_th_pool)
#DATA_SET='terabase_soil'

shopt -u nullglob

srun checkm.sh ${IN_ARR[$SLURM_ARRAY_TASK_ID]} $2 $3
