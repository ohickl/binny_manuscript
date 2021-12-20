#!/bin/bash -l

PROJ_DIR="/scratch/users/ohickl/binning/binny_eval"

DATA_SET='mb_genomes'

if [[ $DATA_SET == 'CAMI2_toy_gut' ]]
then
  SAMPLE="${1##*short_read/}"
  SAMPLE=${SAMPLE%/*}
  ASSEMBLY="${1}contigs/anonymous_gsa.fasta"
  ASSEMBLY_MAP="${1}contigs/mapping_pe_sorted.bam"
  CONTIG_DEPTH="${1}contigs/anonymous_gsa_contig_depth_v2.txt"
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
  CONTIG_DEPTH="${1}H_pooled_H_S00?_contig_depth.txt"
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
  PREFIX="terabase_soil"
fi

THREADS=14

OUT_DIR="${PROJ_DIR}/maxbin_out/${BENCHMARK_DATA_SET}/${SAMPLE}"
#OUT_DIR="${PROJ_DIR}/maxbin_out/${SAMPLE}"


if [[ ! -d $OUT_DIR ]]
then
  echo "Creating out dir."
  mkdir -p $OUT_DIR/bins
fi

conda activate maxbin

if [[ $DATA_SET == 'CAMI1_toy_high_pool' ]]; then
  ls $CONTIG_DEPTH > $OUT_DIR/anonymous_gsa_contig_depth.0.txt
  ABUND_PARAM="-abund_list $OUT_DIR/anonymous_gsa_contig_depth.0.txt"
else
  cat $CONTIG_DEPTH \
      <(awk 'BEGIN {OFS="\t"}; {print $0, 0}' \
      <(diff --new-line-format="" --unchanged-line-format="" \
      <(grep "^>" $ASSEMBLY |sed -e 's/>//g' | sort) \
      <(cut -f1 $CONTIG_DEPTH | sort))) > $OUT_DIR/anonymous_gsa_contig_depth.0.txt
  ABUND_PARAM="-abund $OUT_DIR/anonymous_gsa_contig_depth.0.txt"
fi

# -abund $OUT_DIR/anonymous_gsa_contig_depth.0.txt \
run_MaxBin.pl -contig $ASSEMBLY \
              $ABUND_PARAM \
              -out $OUT_DIR/maxbin \
              -thread $THREADS \
              -min_contig_length 1000 #\
#              -verbose
mv $OUT_DIR/*.fasta $OUT_DIR/bins/

#Create maxbin contig2bin file
for f in $OUT_DIR/bins/*.fasta
do
  sed -e "s|$|\t$f|g" <(grep "^>" $f | sed -e 's|>||g') >> $OUT_DIR/maxbin_contig2bin.txt
done
ln -s $(readlink -f $OUT_DIR/maxbin_contig2bin.txt) $OUT_DIR/maxbin_scaffold2bin.tsv
echo "Contig to bin mapping done"

# Exit if no gsa data is processed.
if [[ $PREFIX == 'terabase_soil' ]]
then
  exit 0
fi

if [[ $DATA_SET == 'CAMI2_toy_gut' || $DATA_SET == 'CAMI1_toy_high' || $DATA_SET == 'CAMI1_toy_high_pool' ]]; then
  # Create amber input file
  conda activate amber
  python3 convert_fasta_bins_to_biobox_format.py $OUT_DIR/bins/*.fasta -o $OUT_DIR/maxbin_bins_amber.tsv
  # Name cleanup for cami2 toy gut samples
  if [[ $PREFIX == cami2_toy_${BENCHMARK_DATA_SET##*_} ]]
  then
    SAMPLE_NR=$(echo $SAMPLE | cut -d "_" -f 3-4)
  elif [[ $PREFIX == cami1_toy_${BENCHMARK_DATA_SET##*_} ]]; then
    SAMPLE_NR="sample_$(echo $SAMPLE | cut -d "0" -f 3)"
  else
    SAMPLE_NR=$SAMPLE
  fi
  sed -r "s|_SAMPLEID_|gsa_${PREFIX}_${SAMPLE_NR}|g" $OUT_DIR/maxbin_bins_amber.tsv > $OUT_DIR/tmp_${SAMPLE} && mv $OUT_DIR/tmp_${SAMPLE} $OUT_DIR/maxbin_bins_amber.tsv
  
  ./get_amber_taxid_v2.py "False" \
                         "cami_data/${BENCHMARK_DATA_SET}/short_read/${SAMPLE}/contigs/gsa_mapping_amber_hr_tax.tsv"\
                         "cami_data/${BENCHMARK_DATA_SET}/short_read/taxonomic_profile_${SAMPLE_NR##*_}.txt" \
                         "$OUT_DIR/maxbin_bins_amber.tsv"
fi

