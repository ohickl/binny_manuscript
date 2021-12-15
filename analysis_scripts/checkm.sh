#!/bin/bash -l

BINNER=${1%%_out*}
# SAMPLE="${1%/*}"
SAMPLE="${1%_exp_binny*}"
SAMPLE="${SAMPLE##*/}"
DATA_SET="$2"

if [[ $DATA_SET == 'CAMI2_toy_gut' ]]; then
  BENCHMARK_DATA_SET="${1##*_out/}"
  BENCHMARK_DATA_SET=${BENCHMARK_DATA_SET%%/*}
  PREFIX="cami2_toy_${BENCHMARK_DATA_SET##*_}"
elif [[ $DATA_SET == 'CAMI1_toy_high' ]]; then
  BENCHMARK_DATA_SET="${1##*_out/}"
  BENCHMARK_DATA_SET=${BENCHMARK_DATA_SET%%/*}
  PREFIX="cami1_toy_${BENCHMARK_DATA_SET##*_}"
elif [[ $DATA_SET == 'CAMI1_toy_high_pool' ]]; then
  SAMPLE="c1_th_pool"
  PREFIX="cami1_toy_high"
elif [[ $DATA_SET == 'pool' ]]; then
  SAMPLE="${1%%/pool/*}"
  SAMPLE="${SAMPLE##*_out/}_pool"
  BENCHMARK_DATA_SET="${1##*_out/}"
  echo BENCHMARK_DATA_SET: $BENCHMARK_DATA_SET
  BENCHMARK_DATA_SET=${BENCHMARK_DATA_SET%%/*}
  echo BENCHMARK_DATA_SET: $BENCHMARK_DATA_SET
  CAMI_DS=${BENCHMARK_DATA_SET%%_*}
  echo CAMI_DS: $CAMI_DS
  if [[ $CAMI_DS == 'c1' ]]; then
    PREFIX_DS="cami1"
  elif [[ $CAMI_DS == 'c2' ]]; then
    PREFIX_DS="cami2"
  fi
  PREFIX="${PREFIX_DS}_toy_${BENCHMARK_DATA_SET##*_}"
elif [[ $DATA_SET == 'mb_genomes' ]]; then
  ASSEMBLY="$(ls img_data/${SAMPLE}/Assembly/IMG_Data/*.fna)"
  ASSEMBLY_MAP="img_data/${SAMPLE}/Assembly/IMG_Data/mapping_pe_sorted.bam"
  CONTIG_DEPTH="img_data/${SAMPLE}/Assembly/IMG_Data/contig_depth.txt"
  BENCHMARK_DATA_SET="img_genomes"
  PREFIX="img_genomes"
elif [[ $DATA_SET == '100s' ]]; then
  PREFIX="100s"
elif [[ $DATA_SET == 'terabase_soil' ]]; then
  PREFIX="terabase_soil"
fi

conda activate metawrap

output_folder=${1%/*} # _exp_binny

echo Sample: ${SAMPLE}
echo CheckM output table: ${output_folder}/checkm/output.txt
echo bins folder: ${output_folder}/bins
echo CheckM folder: ${output_folder}/checkm
echo DATA_SET: $DATA_SET
echo BENCHMARK_DATA_SET: $BENCHMARK_DATA_SET
echo CAMI_DS: $CAMI_DS
echo PREFIX_DS: $PREFIX_DS
echo Cores: $3

if [[ -d ${output_folder}/checkm ]]; then
  echo "Checkm folder already present. Exiting."
  exit 0
fi

if [[ $DATA_SET == 'CAMI1_toy_high' ]] || [[ $CAMI_DS == 'c1' ]]; then
  echo "CAMI 1 data set detected, circumventing | naming scheme for contigs."
  for bin_fasta in ${output_folder}/bins/*.fasta; do
    bin_name=${bin_fasta##*/}
    mv ${bin_fasta} ${output_folder}/bins/${bin_name//|/___}
  done
fi

INPUT_BINS="${output_folder}/bins"

core_mem=1.75
cores=$3
pplacer_threads=$(python -c "print (int(${core_mem}*${cores}/40))")

checkm lineage_wf -t ${cores} \
                  --pplacer_threads ${pplacer_threads} \
                  -x fasta \
                  -f ${output_folder}/checkm/output.txt \
                  ${INPUT_BINS} ${output_folder}/checkm

./get_checkm_bins.py ${output_folder}/

if [[ $DATA_SET == 'CAMI1_toy_high' ]] || [[ $CAMI_DS == 'c1' ]]; then
  for bin_fasta in ${output_folder}/bins_checkm/*.fasta ${output_folder}/bins/*.fasta; do
    bin_name=${bin_fasta##*/}
    mv ${bin_fasta} ${bin_fasta%/*}/${bin_name//___/|}
  done
  sed -i -e 's/___/|/g' ${output_folder}/checkm/output.txt
fi

if [[ $DATA_SET == 'CAMI2_toy_gut' || $DATA_SET == 'CAMI1_toy_high' || $DATA_SET == 'CAMI1_toy_high_pool' || $DATA_SET == 'pool' ]]; then
  # Create amber input file
  echo "Writing Amber format output table."
  conda activate amber_v2
  python3 convert_fasta_bins_to_biobox_format.py ${output_folder}/bins_checkm/*.fasta -o ${output_folder}/${BINNER}_checkm_bins_amber.tsv
  # Name cleanup for cami2 toy gut samples
  if [[ $PREFIX == cami2_toy_${BENCHMARK_DATA_SET##*_} ]]
  then
    SAMPLE_NR=$(echo $SAMPLE | cut -d "_" -f 3-4)
  elif [[ $PREFIX == cami1_toy_${BENCHMARK_DATA_SET##*_} ]]; then
    SAMPLE_NR="sample_$(echo $SAMPLE | cut -d "0" -f 3)"
  else
    SAMPLE_NR=$SAMPLE
  fi

  if [[ $DATA_SET == 'pool' ]]; then
    SAMPLE_NR="sample_pooled"
    GSA_TAX_MAP="cami_data/${BENCHMARK_DATA_SET}/short_read/gsa_anonymous_pooled_amber_hr_tax.tsv"
    FIRST_TAX=(cami_data/${BENCHMARK_DATA_SET}/short_read/taxonomic_profile_*.txt)
    TAX_FILE=${FIRST_TAX[0]}
  else
    GSA_TAX_MAP="cami_data/${BENCHMARK_DATA_SET}/short_read/${SAMPLE}/contigs/gsa_mapping_amber_hr_tax.tsv"
    TAX_FILE="cami_data/${BENCHMARK_DATA_SET}/short_read/taxonomic_profile_${SAMPLE_NR##*_}.txt"
  fi

  sed -r -i "s|_SAMPLEID_|gsa_${PREFIX}_${SAMPLE_NR}|g" ${1}${BINNER}_checkm_bins_amber.tsv

  ./get_amber_taxid_v2.py "False" \
                          ${GSA_TAX_MAP} \
                          ${TAX_FILE} \
                          "${output_folder}/${BINNER}_checkm_bins_amber.tsv"
fi

