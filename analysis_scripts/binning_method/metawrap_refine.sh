#!/bin/bash -l

PROJ_DIR="/scratch/users/ohickl/binning/binny_eval"

DATA_SET='CAMI1_toy_high'
#BENCHMARK_DATA_SET='urogen'

if [[ $DATA_SET == 'CAMI2_toy_gut' ]]
then
  SAMPLE="${1##*short_read/}"
  SAMPLE=${SAMPLE%/*}
  ASSEMBLY="${1}contigs/anonymous_gsa.fasta" #_50
  ASSEMBLY_MAP="${1}contigs/mapping_pe_sorted.bam" #_mm2
  CONTIG_DEPTH="${1}contigs/anonymous_gsa_contig_depth.txt" #_50 #_mm2
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
elif [[ $DATA_SET == 'CAMI1_toy_high_pool' ]]
then
  SAMPLE="c1_th_pool"
  ASSEMBLY="${1}H_pooled_gsa_anonymous.fasta"
  ASSEMBLY_MAP="${1}mapping_pe_sorted.bam"
  CONTIG_DEPTH="${1}H_pooled_contig_depth.txt"
  PREFIX="cami1_toy_high"
elif [[ $DATA_SET == '100s' ]]
then
  SAMPLE=${1%/*}
  ASSEMBLY="${1}anonymous_gsa.fasta"
  ASSEMBLY_MAP="${1}mapping_pe_sorted.bam"
  CONTIG_DEPTH="${1}anonymous_gsa_contig_depth.txt"
  PREFIX="100s"
elif [[ $DATA_SET == 'terabase_soil' ]]
then
  SAMPLE="${1##*_WA-TmG.1.0/}"
  SAMPLE=${SAMPLE%/*}
  ASSEMBLY="${1}assembly.fasta"
  ASSEMBLY_MAP="${1}mapping_pe_sorted.bam"
  CONTIG_DEPTH="${1}contig_depth.txt"
  PREFIX="terabase_soil"
fi


OUT_DIR="${PROJ_DIR}/metawrap_out/${BENCHMARK_DATA_SET}/${SAMPLE}"
PREFIX_1="metawrap_co_ma_me"
PREFIX_2="metawrap_co_ma_me_va"
COMPLETENESS=70
CONTAMINATION=10
BIN_FOLDER="metawrap_${COMPLETENESS}_${CONTAMINATION}_bins"
THREADS=4
MEMORY=280


# binny_out concoct_out maxbin_out metabat_out vamb_out
shopt -s nullglob
BIN_ARR=(concoct_out maxbin_out metabat_out vamb_out)
shopt -u nullglob

if [[ ${#BIN_ARR[*]} == 1 ]] 
then
  mw_input="-A ${BIN_ARR[0]}/${BENCHMARK_DATA_SET}/${SAMPLE}/bins/"
elif [[ ${#BIN_ARR[*]} == 2 ]]
then
  mw_input="-A ${BIN_ARR[0]}/${BENCHMARK_DATA_SET}/${SAMPLE}/bins/ "
  mw_input+="-B ${BIN_ARR[1]}/${BENCHMARK_DATA_SET}/${SAMPLE}/bins/"
elif [[ ${#BIN_ARR[*]} > 2 ]]
then
  mw_input="-A ${BIN_ARR[0]}/${BENCHMARK_DATA_SET}/${SAMPLE}/bins/ "
  mw_input+="-B ${BIN_ARR[1]}/${BENCHMARK_DATA_SET}/${SAMPLE}/bins/ "
  mw_input+="-C ${BIN_ARR[2]}/${BENCHMARK_DATA_SET}/${SAMPLE}/bins/"
fi

# For testing the combination multiple binny setups
# mw_input+=" -B ${BIN_ARR[0]}/${SAMPLE}/bins_no_filter/"


if [[ ! -d $OUT_DIR ]]
then
  echo "Creating out dir."
  mkdir -p $OUT_DIR
fi

conda activate metawrap
echo "Running metawrap refinement with ${COMPLETENESS}% completion threshold."
metawrap bin_refinement -o $OUT_DIR/$PREFIX_1/ \
                        -t $THREADS \
                        -m $MEMORY \
                        ${mw_input} \
                        -c ${COMPLETENESS} \
                        -x ${CONTAMINATION}

# Rename bin fastas back to .fasta ....
for i in $OUT_DIR/$PREFIX_1/${BIN_FOLDER}/*.fa; do
  BIN=$(echo ${i##*/} | cut -d "." -f 2)
  mv $i ${i%/*}/metawrap.${BIN}.fasta
done

# Compress work files
cd $OUT_DIR/$PREFIX_1
rm -r bins/ bins.stats
zip -rm work_files.zip work_files/
cd -

# Create amber input file
echo "Writing Amber format output table."
conda activate amber_v2
python3 convert_fasta_bins_to_biobox_format.py $OUT_DIR/$PREFIX_1/${BIN_FOLDER}/*.fasta -o $OUT_DIR/${PREFIX_1}/${PREFIX_1}_bins_amber.tsv 
# Name cleanup for cami2 toy gut samples
if [[ $PREFIX == cami2_toy_${BENCHMARK_DATA_SET##*_} ]]
then
  SAMPLE_NR=$(echo $SAMPLE | cut -d "_" -f 3-4)
elif [[ $PREFIX == cami1_toy_${BENCHMARK_DATA_SET##*_} ]]; then
  SAMPLE_NR="sample_$(echo $SAMPLE | cut -d "0" -f 3)"
else
  SAMPLE_NR=$SAMPLE
fi
sed -r "s|_SAMPLEID_|gsa_${PREFIX}_${SAMPLE_NR}|g" $OUT_DIR/${PREFIX_1}/${PREFIX_1}_bins_amber.tsv > $OUT_DIR/${PREFIX_1}/tmp_${SAMPLE} && mv $OUT_DIR/${PREFIX_1}/tmp_${SAMPLE} $OUT_DIR/${PREFIX_1}/${PREFIX_1}_bins_amber.tsv 

./get_amber_taxid_v2.py "False" \
                       "cami_data/${BENCHMARK_DATA_SET}/short_read/${SAMPLE}/contigs/gsa_mapping_amber_hr_tax.tsv"\
                       "cami_data/${BENCHMARK_DATA_SET}/short_read/taxonomic_profile_${SAMPLE_NR##*_}.txt" \
                       "$OUT_DIR/${PREFIX_1}/${PREFIX_1}_bins_amber.tsv"

if [[ ${#BIN_ARR[*]} > 3 ]] 
then
  if [[ ${#BIN_ARR[*]} == 4 ]]; then
    mw_input2="-B ${BIN_ARR[3]}/${BENCHMARK_DATA_SET}/${SAMPLE}/bins/"
  elif [[ ${#BIN_ARR[*]} == 5 ]]; then
    mw_input2+="-B ${BIN_ARR[3]}/${BENCHMARK_DATA_SET}/${SAMPLE}/bins/ "
    mw_input2+="-C ${BIN_ARR[4]}/${BENCHMARK_DATA_SET}/${SAMPLE}/bins/"
  fi

  conda activate metawrap
  
  echo "####################################################################################################################################################################################################################################################"
  echo "Running metawrap refinement with ${COMPLETENESS}% completeness threshold, adding leftover binners"
  metawrap bin_refinement -o $OUT_DIR/${PREFIX_2}/ \
                          -t $THREADS \
                          -m $MEMORY \
                          -A $OUT_DIR/$PREFIX_1/${BIN_FOLDER}/ \
                          ${mw_input2} \
                          -c ${COMPLETENESS} \
                          -x ${CONTAMINATION}
  
  # Rename bin fastas back to .fasta ....
  for i in $OUT_DIR/${PREFIX_2}/${BIN_FOLDER}/*.fa; do
    BIN=$(echo ${i##*/} | cut -d "." -f 2)
    mv $i ${i%/*}/metawrap.${BIN}.fasta
  done
  
  # Clean up and compress work files
  cd $OUT_DIR/$PREFIX_2
  rm -r bins/ bins.stats ${BIN_FOLDER}/binsA/
  zip -rm work_files.zip work_files/
  cd -
  
  # Create amber input file
  echo "Writing Amber format output table."
  conda activate amber_v2
  python3 convert_fasta_bins_to_biobox_format.py $OUT_DIR/$PREFIX_2/${BIN_FOLDER}/*.fasta -o $OUT_DIR/${PREFIX_2}/${PREFIX_2}_bins_amber.tsv
  # Name cleanup for cami2 toy gut samples
  if [[ $PREFIX == cami2_toy_${BENCHMARK_DATA_SET##*_} ]]
  then
    SAMPLE_NR=$(echo $SAMPLE | cut -d "_" -f 3-4)
  elif [[ $PREFIX == cami1_toy_${BENCHMARK_DATA_SET##*_} ]]; then
    SAMPLE_NR="sample_$(echo $SAMPLE | cut -d "0" -f 3)"
  else
    SAMPLE_NR=$SAMPLE
  fi
  sed -r "s|_SAMPLEID_|gsa_${PREFIX}_${SAMPLE_NR}|g" $OUT_DIR/${PREFIX_2}/${PREFIX_2}_bins_amber.tsv > $OUT_DIR/${PREFIX_2}/tmp_${SAMPLE} && mv $OUT_DIR/${PREFIX_2}/tmp_${SAMPLE} $OUT_DIR/${PREFIX_2}/${PREFIX_2}_bins_amber.tsv

  ./get_amber_taxid_v2.py "False" \
                         "cami_data/${BENCHMARK_DATA_SET}/short_read/${SAMPLE}/contigs/gsa_mapping_amber_hr_tax.tsv"\
                         "cami_data/${BENCHMARK_DATA_SET}/short_read/taxonomic_profile_${SAMPLE_NR##*_}.txt" \
                         "$OUT_DIR/${PREFIX_2}/${PREFIX_2}_bins_amber.tsv"
fi

echo "Metawrap refinement finished."
