#!/bin/bash -l

PROJ_DIR="/scratch/users/ohickl/binning/binny_eval"

DATA_SET='mb_genomes'

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
  ASSEMBLY="${1}H_pooled_gsa_anonymous.fasta" # _no_pipe
  ASSEMBLY_MAP="${1}mapping_pe_sorted.bam"
  CONTIG_DEPTH="${1}H_pooled_contig_depth.txt" # _no_pipe
  PREFIX="cami1_toy_high"
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
  PREFIX="terabase_soi"l
fi

THREADS=10
CALC_DEPTH="FALSE"

OUT_DIR="${PROJ_DIR}/concoct_out/${BENCHMARK_DATA_SET}/${SAMPLE}"
#OUT_DIR="${PROJ_DIR}/concoct_out/${SAMPLE}"

if [[ ! -d $OUT_DIR ]]
then
  echo "Creating out dir."
  mkdir -p $OUT_DIR/bins
fi

if [[ CALC_DEPTH  == "TRUE"  ]]
then
  if [[ ! -f ${ASSEMBLY_MAP}.bai ]]
  then
    echo "Indexing bam file."
    conda activate mapping
    samtools index -@ 4 ${ASSEMBLY_MAP}
  fi
  
  conda activate concoct
  
  # Cut contigs into smaller parts.
  cut_up_fasta.py $ASSEMBLY \
                  -c 10000 \
                  -o 0 \
                  --merge_last \
                  -b $OUT_DIR/contigs_10K.bed > $OUT_DIR/contigs_10K.fa
  
  # Generate table with coverage depth information per sample and subcontig.
  concoct_coverage_table.py $OUT_DIR/contigs_10K.bed ${1}contigs/*_pe_sorted.bam > $OUT_DIR/coverage_table.tsv
  
  # Run CONCOCT.
  concoct --threads $THREADS \
          --composition_file $OUT_DIR/contigs_10K.fa \
          --coverage_file $OUT_DIR/coverage_table.tsv \
          -b $OUT_DIR/
  
  # Merge subcontig clustering into original contig clustering.
  merge_cutup_clustering.py $OUT_DIR/clustering_gt1000.csv > $OUT_DIR/clustering_merged.csv
  
  # Extract bins as individual FASTA
  extract_fasta_bins.py $ASSEMBLY $OUT_DIR/clustering_merged.csv --output_path $OUT_DIR/bins
else
  conda activate concoct2
  # Format depth file
  if [[ $DATA_SET == 'CAMI1_toy_high_pool' ]]; then
    cat <(printf "contig\tcov_mean_sample_${SAMPLE}_H_S001\tcov_mean_sample_${SAMPLE}_H_S002\tcov_mean_sample_${SAMPLE}_H_S003\tcov_mean_sample_${SAMPLE}_H_S004\tcov_mean_sample_${SAMPLE}_H_S005\n") ${CONTIG_DEPTH} > ${OUT_DIR}/coverage_table.tsv
  else
    cat <(printf "contig\tcov_mean_sample_${SAMPLE}\n") ${CONTIG_DEPTH} > ${OUT_DIR}/coverage_table.tsv
  fi
  # Run CONCOCT.
  concoct --threads ${THREADS} \
          --composition_file ${ASSEMBLY} \
          --coverage_file ${OUT_DIR}/coverage_table.tsv \
          -b $OUT_DIR/
  # Extract bins as individual FASTA
  extract_fasta_bins.py $ASSEMBLY $OUT_DIR/clustering_gt1000.csv --output_path $OUT_DIR/bins
fi

# Exit if no gsa data is processed.
if [[ $PREFIX == 'terabase_soil' ]]
then
  exit 0
fi

for i in $OUT_DIR/bins/*.fa; do
  mv $i ${i%.*}.fasta
done

if [[ $DATA_SET == 'CAMI2_toy_gut' || $DATA_SET == 'CAMI1_toy_high' || $DATA_SET == 'CAMI1_toy_high_pool' ]]; then
  # Create amber input file
  conda activate amber_v2
  python3 convert_fasta_bins_to_biobox_format.py $OUT_DIR/bins/*.fasta -o $OUT_DIR/concoct_bins_amber.tsv
  # Name cleanup for cami2 toy gut samples
  if [[ $PREFIX == cami2_toy_${BENCHMARK_DATA_SET##*_} ]]
  then
    SAMPLE_NR=$(echo $SAMPLE | cut -d "_" -f 3-4)
  elif [[ $PREFIX == cami1_toy_${BENCHMARK_DATA_SET##*_} ]]; then
    SAMPLE_NR="sample_$(echo $SAMPLE | cut -d "0" -f 3)"
  else
    SAMPLE_NR=$SAMPLE
  fi
  sed -r "s|_SAMPLEID_|gsa_${PREFIX}_${SAMPLE_NR}|g" $OUT_DIR/concoct_bins_amber.tsv > $OUT_DIR/tmp_${SAMPLE} && mv $OUT_DIR/tmp_${SAMPLE} $OUT_DIR/concoct_bins_amber.tsv
  
  ./get_amber_taxid_v2.py "False" \
                         "cami_data/${BENCHMARK_DATA_SET}/short_read/${SAMPLE}/contigs/gsa_mapping_amber_hr_tax.tsv"\
                         "cami_data/${BENCHMARK_DATA_SET}/short_read/taxonomic_profile_${SAMPLE_NR##*_}.txt" \
                         "$OUT_DIR/concoct_bins_amber.tsv"
fi

