#!/bin/bash -l

conda activate mapping

DATA_SET='CAMI2_toy_gut'
COASSEMBLY=TRUE

if [[ $DATA_SET == 'CAMI2_toy_gut' ]]
then
  SAMPLE="${1##*short_read/}"
  SAMPLE=${SAMPLE%/*}
  ASSEMBLY="${1}contigs/anonymous_gsa.fasta"
  ASSEMBLY_MAP="${1}contigs/mapping_pe_sorted.bam"
  R_1="${1}reads/anonymous_reads_1.fq"
  R_2="${1}reads/anonymous_reads_2.fq"
  PREFIX="${1}contigs/mapping"
  OUT_DIR="${1}contigs"
elif [[ $DATA_SET == 'CAMI1_toy_high' ]]
then
  SAMPLE="${1##*toy_high/}"
  SAMPLE=${SAMPLE%/*}
  ASSEMBLY="${1}anonymous_gsa.fasta"
  ASSEMBLY_MAP="${1}mapping_pe_sorted.bam"
  R_1="${1}anonymous_reads_1.fq"
  R_2="${1}anonymous_reads_2.fq"
  PREFIX="${1}mapping"
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

if [[ $COASSEMBLY == TRUE && $DATA_SET == 'CAMI2_toy_gut' ]]; then
  BENCHMARK_DATA_SET="${1##*cami_data/}"
  BENCHMARK_DATA_SET=${BENCHMARK_DATA_SET%%/*}
  ASSEMBLY="${1%/20*}/gsa_anonymous_pooled.fasta"
  SAMPLE="${SAMPLE}_pooled"
  PREFIX="${1%/20*}/${SAMPLE}_mapping"
elif [[ $COASSEMBLY == TRUE ]]; then
  ASSEMBLY="${1%/H_*}/H_pooled_gsa_anonymous.fasta"
  SAMPLE="H_pooled_${SAMPLE}"
  PREFIX="${1%/H_*}/H_pooled_${SAMPLE}_mapping"
fi

echo $ASSEMBLY
echo $PREFIX
echo $SAMPLE


LOG="${1}mapping_log"
SAMHEADER="@RG\tID:${SAMPLE}\tSM:MG"

THREADS=4

if [[ ! -f ${ASSEMBLY}.bwt ]]
then
  echo "Indexing assembly."
  bwa index $ASSEMBLY
  echo "Indexing done."
else
  echo "Index already present."
fi

# merge paired and se if se exists
if [[ -f $R_SE ]]; then
  echo Beginning mapping of pe+se reads and merging.
  samtools merge --threads $THREADS -f ${PREFIX}_merged.bam \
                 <(bwa mem -v 1 -t $THREADS -M -R "$SAMHEADER" $ASSEMBLY $R_1 $R_2 2>> $LOG| \
                 samtools view --threads $THREADS -bS -) \
                 <(bwa mem -v 1 -t $THREADS -M -R "$SAMHEADER" $ASSEMBLY $R_SE 2> $LOG| \
                 samtools view --threads $THREADS -bS -) 2>> $LOG
  echo Mapping of pe+se reads and merging done.
else
  echo Beginning mapping of pe-only reads.
  bwa mem -v 1 -t $THREADS -M -R "$SAMHEADER" $ASSEMBLY $R_1 $R_2 2>> $LOG| samtools view --threads $THREADS -bS - > ${PREFIX}_pe.bam
  echo Mapping of pe-only reads done.
fi

# sort
if [[ -f ${PREFIX}_merged.bam ]]; then
  echo Sorting merged pe+se bam file.
  samtools sort --threads $THREADS ${PREFIX}_merged.bam > ${PREFIX}_merged_sorted.bam 2>> $LOG && rm ${PREFIX}_merged.bam
  echo Sorting of merged pe+se bam file done.
else
  echo Sorting pe-only bam file.
  samtools sort --threads $THREADS ${PREFIX}_pe.bam > ${PREFIX}_pe_sorted.bam 2>> $LOG && rm ${PREFIX}_pe.bam
  echo Sorting of pe-only bam file done.
fi

if [[ -f ${PREFIX}_pe_sorted.bam ]]; then
  if [[ -f ${PREFIX}_pe.bam ]]; then
    rm ${PREFIX}_pe.bam
  fi
  echo "Run finished succesfully."
else
  rm ${PREFIX}_pe.bam
  echo "Something went wrong. Final output is missing."
  exit 1
fi
