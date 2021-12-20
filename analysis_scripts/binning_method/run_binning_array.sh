#!/bin/bash -l
#SBATCH -J binny_eval_binning_c2_t_urogen  #c1_th_pool #terabase_soil_v001 # cami2_toy_gut_binny_perp40
#SBATCH --mail-type=fail
#SBATCH --mail-user=oskar.hickl@uni.lu
#SBATCH -n 1
#SBATCH -c 128
# #SBATCH --mem-per-cpu=78GB
#SBATCH --time=1-12:00:00
#SBATCH -p batch
# #SBATCH --qos=normal
#SBATCH --array 0-9%10
#SBATCH -o "slurm_logs/binny/%A-%a-%x-%j.out"

echo "== Starting run at $(date)"
echo "== Job ID: ${SLURM_JOBID}"
echo "== Node list: ${SLURM_NODELIST}"
echo "== Submit dir. : ${SLURM_SUBMIT_DIR}"

shopt -s nullglob
#IN_ARR=(cami_data/c1_t_high/short_read/H_S00?/)
#IN_ARR=(cami_data/cami1_toy_high/) # For pool

# IN_ARR=(cami_data/c2_t_urogen/short_read/201*_sample_*/)

#IN_ARR=(cami_data/c2_t_airways/short_read/2017.12.04_18.56.22_sample_10/)
#IN_ARR=(cami_data/c2_t_oral/short_read/2017.12.04_18.45.54_sample_18/)
#IN_ARR=(cami_data/c2_t_urogen/short_read/) # For VAMB
#IN_ARR=(cami_data/short_read/2017.12.04_18.45.54_sample_*/)

# bulk img runs
IN_ARR=(img_data/*/)
IN_ARR=(${IN_ARR[*]:0:10})

#IN_ARR=(img_data/) # VAMB
#IN_ARR=(cami_data/XXXX/short_read/) # VAMB pool

# IN_ARR=(img_data/KBSSwiS52_FD/)

# IN_ARR=(img_data/KBSS*/)

# IN_ARR=(img_data/KBSSwiS52_FD/)

#IN_ARR=(100s/)
#IN_ARR=(real_world_data/MetaHipMer_Assembly_WA-TmG.1.0/WA-?/)
shopt -u nullglob
#IN_ARR=("${IN_ARR[@]:0:2}" "${IN_ARR[@]:4:2}" "${IN_ARR[9]}")

if [ "$1" == "binny" ] 
then
  echo "Running Binny."
  srun runBinny.sh ${IN_ARR[$SLURM_ARRAY_TASK_ID]} $2 $3
elif [ "$1" == "concoct" ]
then
  echo "Running CONCOCT." 
  srun concoct.sh ${IN_ARR[$SLURM_ARRAY_TASK_ID]}
elif [ "$1" == "maxbin" ]
then
  echo "Running MaxBin."
  srun maxbin.sh ${IN_ARR[$SLURM_ARRAY_TASK_ID]}
elif [ "$1" == "metabat" ]
then
  echo "Running MetaBAT."
  srun metabat.sh ${IN_ARR[$SLURM_ARRAY_TASK_ID]} $2
elif [ "$1" == "vamb" ]
then
  echo "Running VAMB."
  if [[ ! -z $2 ]]; then
    shift # Get rid of binner var
    srun vamb.sh ${IN_ARR[$SLURM_ARRAY_TASK_ID]} $@
  else
    srun vamb.sh ${IN_ARR[$SLURM_ARRAY_TASK_ID]}
  fi
else
  echo "Unsupported argument. Exiting."
  exit 1
fi
