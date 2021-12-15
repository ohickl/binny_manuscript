#!/bin/bash -l

PROJ_DIR="/scratch/users/ohickl/binning/binny_eval"

DATA_SET="${2}"

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
elif [[ $DATA_SET == 'CAMI1_toy_high_pool' ]]
then
  SAMPLE="c1_th_pool"
  ASSEMBLY="${1}H_pooled_gsa_anonymous.fasta"
  ASSEMBLY_MAP="${1}mapping_pe_sorted.bam"
  CONTIG_DEPTH="${1}H_pooled_contig_depth.txt"
  PREFIX="cami1_toy_high"
elif [[ $DATA_SET == 'pool' ]]; then
  echo "Data set is $DATA_SET"
  SAMPLE="${1%%/short_read/*}"
  SAMPLE="${SAMPLE##*cami_data/}_pool"
  ASSEMBLY="${1}gsa_anonymous_pooled.fasta"
  ASSEMBLY_MAP="${1}mapping_pe_sorted.bam"
  CONTIG_DEPTH="${1}anonymous_gsa_pooled_contig_depth.txt"
  BENCHMARK_DATA_SET="${1##*cami_data/}"
  BENCHMARK_DATA_SET=${BENCHMARK_DATA_SET%%/*}
  CAMI_DS=${BENCHMARK_DATA_SET%%_*}
  if [[ $CAMI_DS == 'c1' ]]; then
    PREFIX_DS="cami1"
  elif [[ $CAMI_DS == 'c2' ]]; then
    PREFIX_DS="cami2"
  fi
  PREFIX="${PREFIX_DS}_toy_${BENCHMARK_DATA_SET##*_}"
elif [[ $DATA_SET == 'mb_genomes' ]]
then
  SAMPLE="${1##*img_data/}"
  SAMPLE=${SAMPLE%/*}
  ASSEMBLY="$(ls ${1}Assembly/IMG_Data/*.fna)"
  ASSEMBLY_MAP="${1}Assembly/IMG_Data/mapping_pe_sorted.bam"
  CONTIG_DEPTH="${1}Assembly/IMG_Data/contig_depth.txt"
  BENCHMARK_DATA_SET="img_genomes"
  PREFIX="img_genomes"
elif [[ $DATA_SET == '100s' ]]
then
  SAMPLE=${1%/*}
  ASSEMBLY="${1}anonymous_gsa_50.fasta"
  ASSEMBLY_MAP="${1}mapping_pe_sorted.bam"
  CONTIG_DEPTH="${1}anonymous_gsa_contig_depth_50.txt"
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

THREADS=128

#OUT_DIR="${PROJ_DIR}/metabat_out/${BENCHMARK_DATA_SET}/${SAMPLE}"
OUT_DIR="${PROJ_DIR}/metabat_out/${BENCHMARK_DATA_SET}/pool"
#OUT_DIR="${PROJ_DIR}/metabat_out/${SAMPLE}_def_min_cont"

echo $DATA_SET
echo $SAMPLE
echo $ASSEMBLY
echo $ASSEMBLY_MAP
echo $CONTIG_DEPTH
echo $BENCHMARK_DATA_SET
echo $PREFIX

if [[ ! -d $OUT_DIR ]]
then
  echo "Creating out dir."
  mkdir -p $OUT_DIR/bins
fi

# Create tmp dir
mkdir $OUT_DIR/tmp
TMP_FILE=$(mktemp --tmpdir=$OUT_DIR/tmp -t "metabat_depth_XXXXXX.tsv")
echo "\# fakeheader" > $TMP_FILE
cat "$CONTIG_DEPTH"  >> $TMP_FILE

# Run MetaBat
conda activate metabat

metabat2 -i $ASSEMBLY \
         -a $TMP_FILE \
         --cvExt \
         --saveCls \
         -o $OUT_DIR/metabat \
         -t ${THREADS} \
         -m 2500 \
         -v

ln -s $(readlink -f $OUT_DIR/metabat) $OUT_DIR/metabat.scaffolds2bin.tsv

# Move bin fastas and rename.
for i in $OUT_DIR/*.fa; do
  mv $i ${i%.*}.fasta && mv ${i%.*}.fasta $OUT_DIR/bins
done 

# Exit if no gsa data is processed.
if [[ $PREFIX == 'terabase_soil' ]]
then
  exit 0
fi

if [[ $DATA_SET == 'CAMI2_toy_gut' || $DATA_SET == 'CAMI1_toy_high' || $DATA_SET == 'CAMI1_toy_high_pool' || $DATA_SET == 'pool' ]]; then
  # Create amber input file
  echo "Writing Amber format output table."
  conda activate amber_v2
  python3 convert_fasta_bins_to_biobox_format.py $OUT_DIR/bins/*.fasta -o $OUT_DIR/metabat_bins_amber.tsv
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

  sed -r "s|_SAMPLEID_|gsa_${PREFIX}_${SAMPLE_NR}|g" $OUT_DIR/metabat_bins_amber.tsv > $OUT_DIR/tmp_${SAMPLE} && mv $OUT_DIR/tmp_${SAMPLE} $OUT_DIR/metabat_bins_amber.tsv

  ./get_amber_taxid_v2.py "False" \
                          ${GSA_TAX_MAP} \
                          ${TAX_FILE} \
                          "$OUT_DIR/metabat_bins_amber.tsv"
fi
