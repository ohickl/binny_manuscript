#!/bin/bash -l

DATA_SET='CAMI2_toy_gut'

if [[ $DATA_SET == 'CAMI2_toy_gut' ]]
then
  SAMPLE="${1##*short_read/}"
  SAMPLE=${SAMPLE%/*}
  ASSEMBLY="${1}contigs/anonymous_gsa.fasta"
  ASSEMBLY_MAP="${1}contigs/mapping_pe_sorted.bam"
  CONTIG_DEPTH="${1}contigs/anonymous_gsa_contig_depth.txt"
  BENCHMARK_DATA_SET="${1##*cami_data/}"
  BENCHMARK_DATA_SET=${BENCHMARK_DATA_SET%%/*}
  PREFIX="cami2_toy_${BENCHMARK_DATA_SET##*_}"
elif [[ $DATA_SET == 'CAMI1_toy_high' ]]
then
  SAMPLE="${1##*short_read/}"
  SAMPLE=${SAMPLE%/*}
  ASSEMBLY="${1}contigs/anonymous_gsa.fasta"
  ASSEMBLY_MAP="${1}contigs/mapping_pe_sorted.bam"
  CONTIG_DEPTH="${1}contigs/anonymous_gsa_contig_depth_v2.txt"
  BENCHMARK_DATA_SET="${1##*cami_data/}"
  BENCHMARK_DATA_SET=${BENCHMARK_DATA_SET%%/*}
  PREFIX="cami1_toy_${BENCHMARK_DATA_SET##*_}"
elif [[ $DATA_SET == '100s' ]]
then
  SAMPLE=${1%/*}
  ASSEMBLY="${1}anonymous_gsa.fasta"
  ASSEMBLY_MAP="${1}mapping_pe_sorted.bam"
  CONTIG_DEPTH="${1}anonymous_gsa_contig_depth_v2.txt"
  PREFIX="100s"
fi

OUT_NAME="das_tool_co_ma_me_va"
OUT_DIR="das_tool_out/${BENCHMARK_DATA_SET}/${SAMPLE}/${OUT_NAME}"
THREADS=14

if [[ ! -d $OUT_DIR ]]
then
  echo "Creating out dir."
  mkdir -p ${OUT_DIR}/bins
fi

shopt -s nullglob
# binny_out maxbin_out metabat_out concoct_out vamb_out
BIN_ARR=(concoct_out maxbin_out metabat_out vamb_out)
shopt -u nullglob

declare -a dt_input_int
declare -a dt_names_int

for i in ${BIN_ARR[*]}; do
  if [[ $i == 'binny_out' ]]; then
    dt_input_int+=("${i}/${BENCHMARK_DATA_SET}/${SAMPLE}/${i%_out}.scaffolds2bin.tsv")
    dt_names_int+=("${i%_out}")
  else
    dt_input_int+=("${i}/${BENCHMARK_DATA_SET}/${SAMPLE}/${i%_out}.scaffolds2bin.tsv")
    dt_names_int+=("${i%_out}")
  fi
done

IFS=\, eval 'dt_input="${dt_input_int[*]}"'
IFS=\, eval 'dt_names="${dt_names_int[*]}"'

conda activate dastool

# Prepare input table for DAS_tool.
for i in ${BIN_ARR[@]}; do
  if [[ $i == 'binny_out' ]]; then
    Fasta_to_Scaffolds2Bin.sh -i ${i}/${BENCHMARK_DATA_SET}/${SAMPLE}/bins \
                              -e fasta >  ${i}/${BENCHMARK_DATA_SET}/${SAMPLE}/${i%%_*}.scaffolds2bin.tsv
  else
    Fasta_to_Scaffolds2Bin.sh -i ${i}/${BENCHMARK_DATA_SET}/${SAMPLE}/bins \
                              -e fasta >  ${i}/${BENCHMARK_DATA_SET}/${SAMPLE}/${i%%_*}.scaffolds2bin.tsv
  fi
done

# concoct_out/${SAMPLE}/concoct.scaffolds2bin.tsv
# binny_out/${SAMPLE}/binny.scaffolds2bin.tsv
# binny_out/${SAMPLE}/binny.scaffolds2bin.tsv,maxbin_out/${SAMPLE}/maxbin.scaffolds2bin.tsv,metabat_out/${SAMPLE}/metabat.scaffolds2bin.tsv

DAS_Tool -i ${dt_input[*]} \
         -l ${dt_names[*]} \
         -c $ASSEMBLY \
         -o $OUT_DIR/${OUT_NAME} \
         --threads $THREADS \
         --search_engine diamond

# Get das_tool bin fastas.
conda activate base
#mkdir $OUT_DIR/bins
python get_das_tool_bins.py $ASSEMBLY $OUT_DIR/${OUT_NAME}_DASTool_scaffolds2bin.txt $OUT_DIR/bins/

if [[ $DATA_SET == 'CAMI2_toy_gut' || $DATA_SET == 'CAMI1_toy_high' || $DATA_SET == 'CAMI1_toy_high_pool' ]]; then
  # Create amber input file
  echo "Writing Amber format output table."
  conda activate amber_v2
  python3 convert_fasta_bins_to_biobox_format.py $OUT_DIR/bins/*.fasta -o $OUT_DIR/${OUT_NAME}_bins_amber.tsv
  # Name cleanup for cami2 toy gut samples
  if [[ $PREFIX == cami2_toy_${BENCHMARK_DATA_SET##*_} ]]
  then
    SAMPLE_NR=$(echo $SAMPLE | cut -d "_" -f 3-4)
  elif [[ $PREFIX == cami1_toy_${BENCHMARK_DATA_SET##*_} ]]; then
    SAMPLE_NR="sample_$(echo $SAMPLE | cut -d "0" -f 3)"
  else
    SAMPLE_NR=$SAMPLE
  fi
  sed -r "s|_SAMPLEID_|gsa_${PREFIX}_${SAMPLE_NR}|g" $OUT_DIR/${OUT_NAME}_bins_amber.tsv > $OUT_DIR/tmp_${SAMPLE} && mv $OUT_DIR/tmp_${SAMPLE} $OUT_DIR/${OUT_NAME}_bins_amber.tsv

  ./get_amber_taxid_v2.py "False" \
                         "cami_data/${BENCHMARK_DATA_SET}/short_read/${SAMPLE}/contigs/gsa_mapping_amber_hr_tax.tsv"\
                         "cami_data/${BENCHMARK_DATA_SET}/short_read/taxonomic_profile_${SAMPLE_NR##*_}.txt" \
                         "$OUT_DIR/${OUT_NAME}_bins_amber.tsv"
fi

