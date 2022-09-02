#!/bin/bash -l

PROJ_DIR="/scratch/users/ohickl/binning/binny_eval"

echo $ROOT_PATH
echo $BINNY_CONFIG_SOURCE

DATA_SET="${2}"

if [[ $DATA_SET == 'CAMI2_toy_gut' ]]; then
  SAMPLE="${1##*short_read/}"
  SAMPLE=${SAMPLE%/*}
  ASSEMBLY="${1}contigs/anonymous_gsa.fasta"
  ASSEMBLY_MAP="${1}contigs/mapping_pe_sorted.bam"
  CONTIG_DEPTH="${1}contigs/anonymous_gsa_contig_depth.txt"
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
  ASSEMBLY_MAP=(${1}*mapping_pe_sorted.bam)
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
  OUT_DIR="${PROJ_DIR}/semibin_out/${BENCHMARK_DATA_SET}/pool"
else
  OUT_DIR="${PROJ_DIR}/semibin_out/${BENCHMARK_DATA_SET}/${SAMPLE}"
fi

if [[ ! -d $OUT_DIR ]]
then
  echo "Creating out dir."
  mkdir -p $OUT_DIR
fi


if [[ -d "${OUT_DIR}/bins" ]]; then
  echo "Bin dir ${OUT_DIR}/bins already exists. Delete dir to rerun."
  exit 1
fi

# IMG env arrays
ww_array=('KloWWTslC11_HNv2_FD' 'AD_JPNAS1_MetaG_FD' 'AD_JPNAS3_MetaG_FD' 'AD_JPNHG2_MetaG_FD' 'AD_UKC034_MetaG_FD' \
          'AD_UKC035_MetaG_FD' 'EBPBio102011_DNA_FD' 'EBPBio122009_DNA_FD' 'EBPBio132013_DNA_FD' 'EBPBio232012_DNA_FD' \
          'EBPBio242008_DNA_FD' 'TNRBio252014_DNA_FD')
oc_array=('DEB19_08_DNA_4_FD' 'Colrivmet1449B02_FD' 'Colrivmet1555A02_FD' 'Colrivmeta1546C3_FD' 'SaaInlV_135m_DNA_FD' \
          'MalvireepMed_906_FD' 'ETNvirome_2_1000_FD' 'Browns_ThreeSqA__2_FD' 'China_Galinas_PW_3_FD' 'China_Galinas_PW_6_FD')
so_array=('EasRivMeERMZT100_FD' 'MedBloCAOM2_O1_FD' 'HubBroMetageWRM3_FD' 'RifOxyC2_FD' 'NGEPerone_FD' 'Sb_taG_9_FD' 'Coa_R1_10_FD' \
          'Coa_R1_7_FD' 'Coa_R2_2_FD' 'Coa_R2_FD' 'KBSS13metaG_FD' 'KBSSwiS32_FD' 'KBSSwiS52_FD')

conda activate semibin

if [[ $DATA_SET == 'pool' ]]; then
  model=""
elif [[ $BENCHMARK_DATA_SET == 'c2_t_gi' ]]; then
  model="--environment human_gut"
elif [[ $BENCHMARK_DATA_SET == 'c2_t_oral' ]]; then
  model="--environment human_oral"
elif [[ " ${ww_array[*]} " =~ " ${SAMPLE} " ]]; then
  model="--environment wastewater"
elif [[ " ${oc_array[*]} " =~ " ${SAMPLE} " ]]; then
  model="--environment ocean"
elif [[ " ${so_array[*]} " =~ " ${SAMPLE} " ]]; then
  model="--environment soil"
else
  model="--environment global"
fi

db_dir="/mnt/isilon/projects/ecosystem_biology/user/ohickl/imp3_test_db/gtdb_semibin/"

# Overwrite 64 core limit
NUMEXPR_MAX_THREADS=$THREADS
export NUMEXPR_MAX_THREADS

# Run SemiBin
set -x
SemiBin -v

# --min-len 2500 # for oom samples

SemiBin \
  single_easy_bin \
    --input-fasta $ASSEMBLY \
    --input-bam ${ASSEMBLY_MAP[@]} \
    $model \
    --threads $THREADS \
    --random-seed 0 \
    --reference-db-data-dir $db_dir \
    --output $OUT_DIR
set +x

# Move bin fastas and rename.
mv ${OUT_DIR}/output_recluster_bins ${OUT_DIR}/bins
shopt -s nullglob
for i in ${OUT_DIR}/bins/*.fa; do
  bin_name=$(basename -s .fa ${i})
  mv ${i} ${OUT_DIR}/bins/${bin_name}.fasta
done
shopt -u nullglob

if [[ $PREFIX == 'terabase_soil' ]]
then
  exit 0
fi

if [[ $DATA_SET == 'CAMI2_toy_gut' || $DATA_SET == 'CAMI1_toy_high' || $DATA_SET == 'CAMI1_toy_high_pool' || $DATA_SET == 'pool' ]]; then
  # Create amber input file
  echo "Writing Amber format output table."
  conda activate amber_v2
  python3 convert_fasta_bins_to_biobox_format.py $OUT_DIR/bins/*.fasta -o $OUT_DIR/semibin_bins_amber.tsv
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

  sed -r "s|_SAMPLEID_|gsa_${PREFIX}_${SAMPLE_NR}|g" $OUT_DIR/semibin_bins_amber.tsv > $OUT_DIR/tmp_${SAMPLE} && mv $OUT_DIR/tmp_${SAMPLE} $OUT_DIR/semibin_bins_amber.tsv

  ./get_amber_taxid_v2.py "False" \
                          ${GSA_TAX_MAP} \
                          ${TAX_FILE} \
                          "$OUT_DIR/semibin_bins_amber.tsv"
fi
