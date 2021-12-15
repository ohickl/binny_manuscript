#!/bin/bash -l

DATA_SET='CAMI1_toy_high'

if [[ $DATA_SET == 'CAMI2_toy_gut' ]]
then
  SAMPLE="${1##*short_read/}"
  SAMPLE=${SAMPLE%/*}
  ASSEMBLY="${1}contigs/anonymous_gsa.fasta"
  CONTIG_TO_GENOME="${1}contigs/gsa_mapping.tsv"
  ASSEMBLY_MAP="${1}contigs/mapping_pe_sorted.bam"
  CONTIG_DEPTH="${1}contigs/anonymous_gsa_contig_depth_v2.txt"
  PREFIX="cami2_toy_gut"
elif [[ $DATA_SET == 'CAMI1_toy_high' ]]
then
  SAMPLE="${1##*toy_high/}"
  SAMPLE=${SAMPLE%/*}
  ASSEMBLY="${1}anonymous_gsa.fasta"
  CONTIG_TO_GENOME="${1}gsa_mapping.tsv"
  ASSEMBLY_MAP="${1}mapping_pe_sorted.bam"
  CONTIG_DEPTH="${1}anonymous_gsa_contig_depth_v2.txt"
  PREFIX="cami1_toy_high"
elif [[ $DATA_SET == '100s' ]]
then
  SAMPLE=${1%/*}
  ASSEMBLY="${1}anonymous_gsa.fasta"
  CONTIG_TO_GENOME="${1}gsa_mapping.tsv"
  ASSEMBLY_MAP="${1}mapping_pe_sorted.bam"
  CONTIG_DEPTH="${1}anonymous_gsa_contig_depth_v2.txt"
  PREFIX="100s"
fi

OUT_DIR="cami_gsa_fastas/${SAMPLE}"

if [[ ! -d $OUT_DIR ]]
then
  echo "Creating out dir."
  mkdir -p ${OUT_DIR}/bins
fi

./get_cami_gsa_fastas.py ${ASSEMBLY} ${CONTIG_TO_GENOME} ${OUT_DIR}/bins/ 

