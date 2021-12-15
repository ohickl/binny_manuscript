#!/bin/bash

DATA_SET='CAMI2_toy_gut'
COASSEMBLY=FALSE

if [[ $DATA_SET == 'CAMI2_toy_gut' ]]
then
  SAMPLE="${1##*short_read/}"
  SAMPLE=${SAMPLE%/*}
  ASSEMBLY="${1}contigs/anonymous_gsa.fasta"
  ASSEMBLY_MAP="${1}contigs/mapping_pe_sorted.bam"
  R_1="${1}reads/anonymous_reads_1.fq"
  R_2="${1}reads/anonymous_reads_2.fq"
  PREFIX="cami2_toy_gut"
  OUT_DIR="${1}contigs"
elif [[ $DATA_SET == 'CAMI1_toy_high' ]]
then
  SAMPLE="${1##*toy_high/}"
  SAMPLE=${SAMPLE%/*}
  ASSEMBLY="${1}anonymous_gsa.fasta"
  ASSEMBLY_MAP="${1}mapping_pe_sorted.bam"
  R_1="${1}anonymous_reads_1.fq"
  R_2="${1}anonymous_reads_2.fq"
  PREFIX="cami1_toy_high"
  OUT_DIR="${1}"
elif [[ $DATA_SET == 'mb_genomes' ]]
then
  SAMPLE="${1##*img_data/}"
  SAMPLE=${SAMPLE%/*}
  ASSEMBLY="$(ls ${1}Assembly/IMG_Data/*.fna)"
  ASSEMBLY_MAP="${1}Assembly/IMG_Data/mapping_pe_sorted.bam"
  if [[ (-d "${1}Sequence/Filtered_Raw_Data") ]]; then
    raw_data="${1}Sequence/Filtered_Raw_Data"
  else
    raw_data="${1}Sequence/QA_Filtered_Raw_Data"
  fi
  R_1="${raw_data}/reads_1.fq"
  R_2="${raw_data}/reads_2.fq"
  PREFIX="${1}Assembly/IMG_Data/mapping"
  OUT_DIR="${1}Assembly/IMG_Data"
fi

if [[ $DATA_SET == 'mb_genomes' ]]; then
  OUT_DEPTH_FILE="contig_depth.txt"
else
  OUT_DEPTH_FILE="anonymous_gsa_contig_depth.txt"
fi

OUTPUT="$OUT_DIR/$OUT_DEPTH_FILE"


if [[ $COASSEMBLY == TRUE && $DATA_SET == 'CAMI2_toy_gut' ]]; then
  BENCHMARK_DATA_SET="${1##*cami_data/}"
  BENCHMARK_DATA_SET=${BENCHMARK_DATA_SET%%/*}
  ASSEMBLY="${1%/20*}/gsa_anonymous_pooled.fasta"
  SAMPLE="${SAMPLE}_pooled"
  PREFIX="${1%/20*}/${SAMPLE}_mapping"
elif [[ $COASSEMBLY == TRUE ]]; then
  ASSEMBLY="${1%/H_*}/H_pooled_gsa_anonymous.fasta"
  SAMPLE="H_pooled_${SAMPLE}"
  PREFIX="${1%/H_*}/H_pooled_${SAMPLE}"
  ASSEMBLY_MAP="${PREFIX}_mapping_pe_sorted.bam"
  OUTPUT="${PREFIX}_contig_depth.txt"
fi

echo $SAMPLE
# Check if dirs are present
if [[ ! -d $OUT_DIR ]]
then
  echo "Creating $OUT_DIR."
  mkdir -p $OUT_DIR
fi

#conda activate bedtools

echo "Running BEDTools for average depth in each position"
TMP_DEPTH=$(mktemp --tmpdir=$OUT_DIR/tmp -t "depth_file_XXXXXX.txt")
genomeCoverageBed -ibam "$ASSEMBLY_MAP" | grep -v "genome" > $TMP_DEPTH
echo "Depth calculation done"

## This method of depth calculation was adapted and modified from the CONCOCT code
./calcAvgCoverage.pl $TMP_DEPTH $ASSEMBLY > $OUTPUT 

echo "Removing the temporary file"
rm $TMP_DEPTH
