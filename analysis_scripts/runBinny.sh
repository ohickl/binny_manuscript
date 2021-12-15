#!/bin/bash -l

#ROOT_PATH="/work/projects/ecosystem_biology/local_tools/IMP3" # IMP binny
binny_ver="20210922" # 20210624 20210922
ROOT_PATH="/scratch/users/ohickl/binning/tools/binny_alpha_006_${binny_ver}" # Binny standalone 20210624

PROJ_DIR="/scratch/users/ohickl/binning/binny_eval"

if [[ $ROOT_PATH == "/scratch/users/ohickl/binning/tools/binny_alpha_006_${binny_ver}" ]]
then
  BINNY_CONFIG_SOURCE="${PROJ_DIR}/config_binny_alpha_007.yaml" # Python Binny
else
  BINNY_CONFIG_SOURCE="${PROJ_DIR}/config.binny.yaml" # IMP binny
fi

DATA_SET="${2}"
#BENCHMARK_DATA_SET='urogen'

if [[ $DATA_SET == 'CAMI2_toy_gut' ]]; then
  SAMPLE="${1##*short_read/}"
  SAMPLE=${SAMPLE%/*}
  ASSEMBLY="${1}contigs/anonymous_gsa.fasta" #_50
  ASSEMBLY_MAP="${1}contigs/mapping_pe_sorted.bam" #_mm2
  CONTIG_DEPTH="${1}contigs/anonymous_gsa_contig_depth.txt" #_50 #_mm2
  BENCHMARK_DATA_SET="${1##*cami_data/}"
  BENCHMARK_DATA_SET=${BENCHMARK_DATA_SET%%/*}
  PREFIX="cami2_toy_${BENCHMARK_DATA_SET##*_}"
elif [[ $DATA_SET == 'CAMI1_toy_high' ]]; then
  SAMPLE="${1##*short_read/}"
  SAMPLE=${SAMPLE%/*}
  ASSEMBLY="${1}contigs/anonymous_gsa.fasta"
  ASSEMBLY_MAP="${1}contigs/mapping_pe_sorted.bam"
  CONTIG_DEPTH="${1}contigs/anonymous_gsa_contig_depth_v2.txt"
  BENCHMARK_DATA_SET="${1##*cami_data/}"
  BENCHMARK_DATA_SET=${BENCHMARK_DATA_SET%%/*}
  PREFIX="cami1_toy_${BENCHMARK_DATA_SET##*_}"
elif [[ $DATA_SET == 'CAMI1_toy_high_pool' ]]; then
  SAMPLE="c1_th_pool"
  ASSEMBLY="${1}H_pooled_gsa_anonymous.fasta"
  ASSEMBLY_MAP="${1}mapping_pe_sorted.bam"
  CONTIG_DEPTH="${1}H_pooled_contig_depth.txt"
  PREFIX="cami1_toy_high"
elif [[ $DATA_SET == 'pool' ]]; then
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
elif [[ $DATA_SET == 'mb_genomes' ]]; then
  SAMPLE="${1##*img_data/}"
  SAMPLE=${SAMPLE%/*}
  ASSEMBLY="$(ls ${1}Assembly/IMG_Data/*.fna)"
  ASSEMBLY_MAP="${1}Assembly/IMG_Data/mapping_pe_sorted.bam"
  CONTIG_DEPTH="${1}Assembly/IMG_Data/contig_depth.txt"
  BENCHMARK_DATA_SET="img_genomes"
  PREFIX="img_genomes"
elif [[ $DATA_SET == '100s' ]]; then
  SAMPLE=${1%/*}
  ASSEMBLY="${1}anonymous_gsa.fasta"
  ASSEMBLY_MAP="${1}mapping_pe_sorted.bam"
  CONTIG_DEPTH="${1}anonymous_gsa_contig_depth.txt"
  PREFIX="100s"
elif [[ $DATA_SET == 'terabase_soil' ]]; then
  SAMPLE="${1##*_WA-TmG.1.0/}"
  SAMPLE=${SAMPLE%/*}
  ASSEMBLY="${1}assembly.fasta"
  ASSEMBLY_MAP="${1}mapping_pe_sorted.bam"
  CONTIG_DEPTH="${1}contig_depth.txt"
  PREFIX="terabase_soil"
fi

THREADS="${3}"


if [[ $DATA_SET == 'pool' ]]; then
  OUT_DIR="${PROJ_DIR}/binny_out/${BENCHMARK_DATA_SET}/pool_exp_binny"
else
  OUT_DIR="${PROJ_DIR}/binny_out/${BENCHMARK_DATA_SET}/${SAMPLE}_exp_binny"
fi
#OUT_DIR="${PROJ_DIR}/binny_out/${BENCHMARK_DATA_SET}/${SAMPLE}_exp_binny"
#OUT_DIR="${PROJ_DIR}/binny_out/${BENCHMARK_DATA_SET}/pool_exp_binny"
#OUT_DIR="${PROJ_DIR}/binny_out/c2_t_${BENCHMARK_DATA_SET}/${SAMPLE}_alpha_006"
#OUT_DIR="${PROJ_DIR}/binny_out/${SAMPLE}_alpha_006"

if [[ ! -d $OUT_DIR ]]
then
  echo "Creating out dir."
  mkdir -p $OUT_DIR
#  mkdir -p $OUT_DIR/bins
fi


#    -e "s|pk: 10|pk: 10|g" \
#    -e "s|perp: 30|perp: 30|g" \
#    -e "s|nn: 4|nn: 4|g" \
#    -e "s|Contig_depth: \"\"|Contig_depth: \"$PROJ_DIR/$CONTIG_DEPTH\"|g" \
# Setup config file for sample.
sed -e "s|__SAMPLE__|$SAMPLE|g" \
    -e "s|__PROJ_DIR__|$PROJ_DIR|g" \
    -e "s|__ASSEMBLY__|$ASSEMBLY|g" \
    -e "s|__ASSEMBLY_MAP__|$ASSEMBLY_MAP|g" \
    -e "s|Contig_depth: \"\"|Contig_depth: \"$PROJ_DIR/$CONTIG_DEPTH\"|g" \
    -e "s|__OUT_DIR__|$OUT_DIR|g" \
    -e "s|kmers: '2,3,4'|kmers: '2,3,4'|g" \
    -e "s|cutoff: 1000|cutoff: 1000|g" \
    -e "s|completeness: 80|completeness: 80|g" \
    -e "s|purity: 85|purity: 85|g" \
    -e "s|big_mem_avail: FALSE|big_mem_avail: FALSE|g" \
    -e "s|big_mem_per_core_gb: 26|big_mem_per_core_gb: 26|g" \
    -e "s|__PREFIX__|$PREFIX|g" ${BINNY_CONFIG_SOURCE} > $OUT_DIR/${SAMPLE}_config.yaml

SNAKEFILE="binny_Snakefile"

echo $1

# Run binny.
if [[ $ROOT_PATH == "/scratch/users/ohickl/binning/tools/binny" ]]
then
  conda activate $ROOT_PATH/conda/snakemake_env
  SNAKEFILE="Snakefile"
elif [[ $ROOT_PATH == "/scratch/users/ohickl/binning/tools/binny_alpha_006_${binny_ver}" ]]
then
  conda activate $ROOT_PATH/conda/snakemake_env
  SNAKEFILE="Snakefile"
fi

snakemake --verbose --unlock -s $ROOT_PATH/${SNAKEFILE} --configfile $OUT_DIR/${SAMPLE}_config.yaml -j $THREADS --use-conda --conda-prefix $ROOT_PATH/conda

# For IMP binny: binny_Snakefile
snakemake --verbose --omit-from "zip_output" -s $ROOT_PATH/${SNAKEFILE} --configfile $OUT_DIR/${SAMPLE}_config.yaml -j $THREADS --use-conda --conda-prefix $ROOT_PATH/conda
# --unlock --rerun-incomplete --cleanup-metadata 'intermediary/prokka.faa.markers.hmmscan'
# --omit-from "zip_output"

# Get bin fastas. <-- Obsolete, now integrated into Binny.
#./get_binny_bins_v2.py $OUT_DIR/assembly.fa $OUT_DIR/contigs2bin.tsv $OUT_DIR/bins/

# Exit if no gsa data is processed.
if [[ $PREFIX == 'terabase_soil' ]]
then
  exit 0
fi

if [[ $DATA_SET == 'CAMI2_toy_gut' || $DATA_SET == 'CAMI1_toy_high' || $DATA_SET == 'CAMI1_toy_high_pool' || $DATA_SET == 'pool' ]]; then
  # Create amber input file
  echo "Writing Amber format output table."
  conda activate amber_v2
  python3 convert_fasta_bins_to_biobox_format.py $OUT_DIR/bins/*.fasta -o $OUT_DIR/binny_bins_amber.tsv
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

  sed -r "s|_SAMPLEID_|gsa_${PREFIX}_${SAMPLE_NR}|g" $OUT_DIR/binny_bins_amber.tsv > $OUT_DIR/tmp_${SAMPLE} && mv $OUT_DIR/tmp_${SAMPLE} $OUT_DIR/binny_bins_amber.tsv
  
  ./get_amber_taxid_v2.py "False" \
                          ${GSA_TAX_MAP} \
                          ${TAX_FILE} \
                          "$OUT_DIR/binny_bins_amber.tsv"
fi

