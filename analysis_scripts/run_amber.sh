#!/bin/bash -l
#SBATCH -J binny_eval_amber_cami_all_final_binny_tests_v01 #
#SBATCH --mail-type=fail
#SBATCH --mail-user=oskar.hickl@uni.lu
#SBATCH -n 1
#SBATCH -c 7
#SBATCH --time=0-04:45:00
#SBATCH -p batch
#SBATCH --qos=normal
#SBATCH -o "slurm_logs/amber/%j-%x.out"

echo "== Starting run at $(date)"
echo "== Job ID: ${SLURM_JOBID}"
echo "== Node list: ${SLURM_NODELIST}"
echo "== Submit dir. : ${SLURM_SUBMIT_DIR}"

conda activate amber_v2

AMBER_DIR="/scratch/users/ohickl/binning/tools/amber_v2.0.3"
EVAL_DIR="/scratch/users/ohickl/binning/binny_eval"
OUT_DIR="$EVAL_DIR/amber_out/cami_all_final_binny_tests_v2"
paper_mode=false
GS_MAP="$OUT_DIR/gsa_mappings.tsv"

gs_amber="${OUT_DIR}/Gold-Standard_bins_amber.tsv"

if [[ -f ${gs_amber} ]]; then
  AMBER_INPUT=(${gs_amber})
else
    AMBER_INPUT=()
fi

if [[ $paper_mode == true ]]; then
  echo "Paper mode on."
  for i in binny binny\[CheckM\] binny_1df binny_1df\[CheckM\] binny_no_mask \
           binny_no_mask\[CheckM\] CONCOCT CONCOCT\[CheckM\] MaxBin2 MaxBin2\[CheckM\] \
           MetaBAT2 MetaBAT2\[CheckM\] VAMB VAMB\[CheckM\]; do
             binner_amber="${OUT_DIR}/${i}_bins_amber.tsv"
             if [[ -f ${binner_amber} ]]; then
               AMBER_INPUT+=(${binner_amber})
             fi
  done
else
  AMBER_INPUT=$(ls ${OUT_DIR}/*amber.tsv)
fi

echo ${AMBER_INPUT[*]}

declare -a dt_names_int

for i in ${AMBER_INPUT[*]}; do
  ref_out=${i##*/}
  am_names_int+=("${ref_out%%_bins*}")
done

IFS=\, eval 'am_names="${am_names_int[*]}"'

echo $am_names

${AMBER_DIR}/amber.py -g $GS_MAP \
         -l "${am_names[*]}" \
         -k "circular element" \
         -x "50,70,90" \
         --ncbi_dir "${EVAL_DIR}/20210528_ncbi_taxdump" \
         ${AMBER_INPUT[*]} \
         -o $OUT_DIR

