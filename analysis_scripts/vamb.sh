#!/bin/bash -l

PROJ_DIR="/scratch/users/ohickl/binning/binny_eval"

DATA_SET='pool'
SNAKEMAKE=true

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
elif [[ $DATA_SET == 'CAMI2_toy_gut_all' ]]
then
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
elif [[ $DATA_SET == 'CAMI1_toy_high_all' ]]
then
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

echo $DATA_SET
echo $SAMPLE
echo $ASSEMBLY
echo $ASSEMBLY_MAP
echo $CONTIG_DEPTH
echo $BENCHMARK_DATA_SET
echo $PREFIX

THREADS=128

arg1=$1

if [ ${DATA_SET} == "mb_genomes"  ]; then
  shift # Get rid of input dir var
  proj_samples=( "$@" )
  combined_out=""
  for i in ${proj_samples[*]}; do
    combined_out+="${i}_"
  done
  combined_out+="combined"
  OUT_DIR="${PROJ_DIR}/vamb_out/${BENCHMARK_DATA_SET}/${combined_out}"
else
  OUT_DIR="${PROJ_DIR}/vamb_out/${BENCHMARK_DATA_SET}/combined"
fi
#OUT_DIR="${PROJ_DIR}/vamb_out/${BENCHMARK_DATA_SET}/${SAMPLE}_sm"
#OUT_DIR="${PROJ_DIR}/vamb_out/${SAMPLE}"


if [[ $SNAKEMAKE && $DATA_SET != 'pool' ]]; then
  if [[ ! -d $OUT_DIR ]]
  then
    echo "Creating out dir."
    mkdir -p $OUT_DIR
  fi
  
  if [ ${DATA_SET%%_*} == "CAMI1" ]; then
    SAMPLE_NAMING_SCHEME="H_S00?"
  elif [ ${DATA_SET%%_*} == "CAMI2" ]; then
    SAMPLE_NAMING_SCHEME="201*_sample_*"
  elif [ ${DATA_SET} == "mb_genomes" ]; then
    SAMPLE_NAMING_SCHEME="{"
    for i in ${proj_samples[@]}; do
      SAMPLE_NAMING_SCHEME+="${i},"
    done
    SAMPLE_NAMING_SCHEME=${SAMPLE_NAMING_SCHEME%,*}
    SAMPLE_NAMING_SCHEME+="}"
  fi
  
  cd $OUT_DIR
  VAMB_PATH="/mnt/lscratch/users/ohickl/binning/tools/vamb/workflow"
  SNAKEFILE="${VAMB_PATH}/vamb.snake.conda.py"
  VAMB_ENVS="${VAMB_PATH}/envs"
  MEM=$(expr $THREADS \* 4)

  CONTIG_FILE="${OUT_DIR}/${BENCHMARK_DATA_SET}_contigs.txt"
  S2D_FILE="${OUT_DIR}/${BENCHMARK_DATA_SET}_samples2data.txt"
  CONFIG_FILE="${OUT_DIR}/${BENCHMARK_DATA_SET}_config.json"

  echo "\${proj_dir[*]}: ${proj_samples[*]}"
  echo "\$SAMPLE_NAMING_SCHEME: ${SAMPLE_NAMING_SCHEME}"
  echo "\$OUT_DIR: $OUT_DIR"
  
  if [ ${DATA_SET} == "mb_genomes" ]; then
    sample_wc_string=""
    for i in ${proj_samples[@]}; do
      sample_wc_string+="${PROJ_DIR}/${arg1}${i} "
    done
    sample_wc_string=${sample_wc_string% *}
  else
    sample_wc_string=${PROJ_DIR}/${arg1}${SAMPLE_NAMING_SCHEME}
  fi

  if [[ ! -f $CONTIG_FILE && ! -f $S2D_FILE ]]; then
   for sample in $sample_wc_string; do
      echo "${sample}"
      if [[ ${DATA_SET%%_*} == "CAMI1"  || ${DATA_SET%%_*} == "CAMI2" ]]; then
        SAMPLE="${sample##*short_read/}"
        SAMPLE=${SAMPLE%/*}
        R_1="${sample}/reads/anonymous_reads_1.fq"
        R_2="${sample}/reads/anonymous_reads_2.fq"
        ASSEMBLY="${sample}/contigs/anonymous_gsa.fasta"
      elif [ ${DATA_SET} == "mb_genomes"  ]; then
        #SAMPLE="${sample##*img_data/}"
        # SAMPLE=${SAMPLE%/*}
        SAMPLE=${sample##*/}
        ASSEMBLY="$(ls ${sample}/Assembly/IMG_Data/*.fna)"
        if [[ (-d ${sample}/Sequence/Filtered_Raw_Data) ]]; then
          raw_data="${sample}/Sequence/Filtered_Raw_Data"
        elif [[ (-d ${sample}/Sequence/QC_and_Genome_Assembly) ]]; then
          raw_data="${sample}/Sequence/QC_and_Genome_Assembly"
        else
          raw_data="${sample}/Sequence/QA_Filtered_Raw_Data"
        fi
        R_1="${raw_data}/reads_1.fq"
        R_2="${raw_data}/reads_2.fq"
      fi
      printf "${ASSEMBLY}\n" >> $CONTIG_FILE
      printf "${SAMPLE}\t${R_1}\t${R_2}\n" >> $S2D_FILE
    done
  fi

  # -m 2000 --minfasta 500000
  printf "%s\n" \
         "{" \
         "  \"contigs\": \"${CONTIG_FILE}\"," \
         "  \"sample_data\": \"${S2D_FILE}\"," \
         "  \"index_size\": \"$(echo ${MEM} / 2.5 | bc)G\"," \
         "  \"minimap_mem\": \"${MEM}gb\"," \
         "  \"minimap_ppn\": \"${THREADS}\"," \
         "  \"vamb_mem\": \"${MEM}gb\"," \
         "  \"vamb_ppn\": \"${THREADS}\"," \
         "  \"checkm_mem\": \"${MEM}gb\"," \
         "  \"checkm_ppn\": \"${THREADS}\"," \
         "  \"vamb_params\": \"-o C -m 2000 --minfasta 500000 --outdir ${OUT_DIR}/vamb\"" \
         "}" > $CONFIG_FILE

  # conda activate vamb_sm
  snakemake --cores ${THREADS} --use-conda --conda-prefix $VAMB_ENVS --configfile $CONFIG_FILE --snakefile $SNAKEFILE

  echo 'Restoring original bin as well as contig names.' 
  for bin in ${OUT_DIR}/vamb/bins/*.fna; do
    bin_sample=$(head -n 1 ${bin})
    bin_cluster="${bin##*/}"
    bin_cluster="${bin_cluster%.*}"
    bin_cluster="C${bin_cluster##*C}"
    if [ ${DATA_SET%%_*} == "CAMI1"  ]; then
      bin_sample="S${bin_sample##*CH|S}"
      bin_sample="${bin_sample%|*}"
      sed -i -e "s/S[^C]\{1,2\}CH|S[^C]\{1,2\}/H|${bin_sample}|/g" $bin && mv $bin ${OUT_DIR}/vamb/bins/${bin_sample##*|}${bin_cluster}.fasta
    elif [ ${DATA_SET%%_*} == "CAMI2"  ]; then
      bin_sample="S${bin_sample##*CS}"
      bin_sample="${bin_sample%%C*}" 
      sed -i -e "s/S[^C]\{1,2\}CS[^C]\{1,2\}/${bin_sample}/g" $bin && mv $bin ${OUT_DIR}/vamb/bins/${bin_sample}${bin_cluster}.fasta
    else
      cp $bin ${bin%.*}.fasta 
    fi
  done

  echo 'Moving bins to sample directories.'
  counter=0
  for sample in $sample_wc_string; do
    counter=$(( counter + 1 ))
    if [[ ${DATA_SET%%_*} == "CAMI1"  || ${DATA_SET%%_*} == "CAMI2" ]]; then
      SAMPLE="${sample##*short_read/}"
      SAMPLE=${SAMPLE%/*}
      if [ ${DATA_SET%%_*} == "CAMI1"  ]; then
        SAMPLE_NR=${SAMPLE##*00}
        BIN_PREFIX="S${SAMPLE_NR}C"
      elif [ ${DATA_SET%%_*} == "CAMI2"  ]; then
        SAMPLE_NR=${SAMPLE##*_}
        BIN_PREFIX="S${SAMPLE_NR}C"
      fi
    elif [ ${DATA_SET} == "mb_genomes"  ]; then
      SAMPLE="${sample##*img_data/}"
      SAMPLE=${SAMPLE%/*}
      BIN_PREFIX="S${counter}"
    fi
    SAMPLE_OUT_DIR=${OUT_DIR%/*}/${SAMPLE}
    echo test ${SAMPLE_OUT_DIR} ${OUT_DIR%/*} ${SAMPLE} ${OUT_DIR}
    if [[ ! -d $SAMPLE_OUT_DIR ]]; then
      echo "Creating sample out dir for ${SAMPLE}."
      mkdir -p ${SAMPLE_OUT_DIR}/bins
    fi
    for bin in ${OUT_DIR}/vamb/bins/${BIN_PREFIX}*.fasta; do
      bin_name=${bin##*/}
      cp $bin ${SAMPLE_OUT_DIR}/bins/${bin_name}
    done
       
    if [[ $DATA_SET == 'CAMI2_toy_gut' || $DATA_SET == 'CAMI1_toy_high' || $DATA_SET == 'CAMI1_toy_high_pool' ]]; then
      # Create amber input file
      conda activate amber_v2
      python3 convert_fasta_bins_to_biobox_format.py $SAMPLE_OUT_DIR/bins/*.fasta -o $SAMPLE_OUT_DIR/vamb_bins_amber.tsv
      conda deactivate
      # Name cleanup for cami2 toy gut samples
      if [[ $PREFIX == cami2_toy_${BENCHMARK_DATA_SET##*_} ]]; then
        SAMPLE_NR=$(echo $SAMPLE | cut -d "_" -f 3-4)
      elif [[ $PREFIX == cami1_toy_${BENCHMARK_DATA_SET##*_} ]]; then
        SAMPLE_NR="sample_$(echo $SAMPLE | cut -d "0" -f 3)"
      else
        SAMPLE_NR=$SAMPLE
      fi
      sed -r "s|_SAMPLEID_|gsa_${PREFIX}_${SAMPLE_NR}|g" $SAMPLE_OUT_DIR/vamb_bins_amber.tsv > $SAMPLE_OUT_DIR/tmp_${SAMPLE} && mv $SAMPLE_OUT_DIR/tmp_${SAMPLE} $SAMPLE_OUT_DIR/vamb_bins_amber.tsv
      
      ./get_amber_taxid_v2.py "False" \
                              "cami_data/${BENCHMARK_DATA_SET}/short_read/${SAMPLE}/contigs/gsa_mapping_amber_hr_tax.tsv"\
                              "cami_data/${BENCHMARK_DATA_SET}/short_read/taxonomic_profile_${SAMPLE_NR##*_}.txt" \
                              "$SAMPLE_OUT_DIR/vamb_bins_amber.tsv"
    fi
  done
  
  # zip -r ${OUT_DIR}/vamb/bins.zip ${OUT_DIR}/vamb/bins && rm -r ${OUT_DIR}/vamb/bins
  
  exit 0
elif [[ $SNAKEMAKE && $DATA_SET == 'pool' ]]; then
  OUT_DIR="${PROJ_DIR}/vamb_out/${BENCHMARK_DATA_SET}/pool"
  if [[ ! -d $OUT_DIR ]]
  then
    echo "Creating out dir."
    mkdir -p $OUT_DIR
  fi

  if [ $CAMI_DS == 'c1' ]; then
    SAMPLE_NAMING_SCHEME="H_S00?/"
  elif [ $CAMI_DS == 'c2' ]; then
    SAMPLE_NAMING_SCHEME="201*_sample_*/"
  fi

  cd $OUT_DIR
  VAMB_PATH="/mnt/lscratch/users/ohickl/binning/tools/vamb/workflow"
  SNAKEFILE="${VAMB_PATH}/vamb.snake.conda.py"
  VAMB_ENVS="${VAMB_PATH}/envs"
  MEM=$(expr $THREADS \* 2)

  CONTIG_FILE="${OUT_DIR}/${BENCHMARK_DATA_SET}_contigs.txt"
  S2D_FILE="${OUT_DIR}/${BENCHMARK_DATA_SET}_samples2data.txt"
  CONFIG_FILE="${OUT_DIR}/${BENCHMARK_DATA_SET}_config.json"

  echo "\${proj_dir[*]}: ${proj_samples[*]}"
  echo "\$SAMPLE_NAMING_SCHEME: ${SAMPLE_NAMING_SCHEME}"
  echo "\$OUT_DIR: $OUT_DIR"

  sample_wc_string=${PROJ_DIR}/${arg1}${SAMPLE_NAMING_SCHEME}

  echo "\$sample_wc_string: $sample_wc_string"

  if [[ ! -f $CONTIG_FILE && ! -f $S2D_FILE ]]; then
    for sample in $sample_wc_string; do
      echo "${sample}"
      if [[ $CAMI_DS == 'c1'  || $CAMI_DS == 'c2' ]]; then
        SAMPLE="${sample##*short_read/}"
        SAMPLE=${SAMPLE%/*}
        R_1="${sample}reads/anonymous_reads_1.fq"
        R_2="${sample}reads/anonymous_reads_2.fq"
      fi
      printf "${SAMPLE}\t${R_1}\t${R_2}\n" >> $S2D_FILE
    done
  fi

  printf "${PROJ_DIR}/${ASSEMBLY}\n" >> $CONTIG_FILE

  printf "%s\n" \
         "{" \
         "  \"contigs\": \"${CONTIG_FILE}\"," \
         "  \"sample_data\": \"${S2D_FILE}\"," \
         "  \"index_size\": \"$(echo ${MEM} / 2.5 | bc)G\"," \
         "  \"minimap_mem\": \"${MEM}gb\"," \
         "  \"minimap_ppn\": \"${THREADS}\"," \
         "  \"vamb_mem\": \"${MEM}gb\"," \
         "  \"vamb_ppn\": \"${THREADS}\"," \
         "  \"checkm_mem\": \"${MEM}gb\"," \
         "  \"checkm_ppn\": \"${THREADS}\"," \
         "  \"vamb_params\": \"-m 2000 --minfasta 500000 --outdir ${OUT_DIR}/vamb\"" \
         "}" > $CONFIG_FILE
  
  # conda activate vamb_sm
  snakemake --cores ${THREADS} --use-conda --conda-prefix $VAMB_ENVS --configfile $CONFIG_FILE --snakefile $SNAKEFILE
  
  if [[ ! -d ${OUT_DIR}/bins ]]; then
     mkdir -p ${OUT_DIR}/bins
   fi

  echo 'Restoring original bin as well as contig names.'
  for bin in ${OUT_DIR}/vamb/bins/*.fna; do
    bin_n=$(basename ${bin} .fna)
    cp $bin ${OUT_DIR}/bins/${bin_n}.fasta
    sed -i -e "s/>S1C/>/g" ${OUT_DIR}/bins/${bin_n}.fasta
  done

  cd -

  if [[ $DATA_SET == 'CAMI2_toy_gut' || $DATA_SET == 'CAMI1_toy_high' || $DATA_SET == 'CAMI1_toy_high_pool' || $DATA_SET == 'pool' ]]; then
    # Create amber input file
    echo "Writing Amber format output table."
    conda activate amber_v2
    python3 convert_fasta_bins_to_biobox_format.py $OUT_DIR/bins/*.fasta -o $OUT_DIR/vamb_bins_amber.tsv
    # Name cleanup for cami2 toy gut samples
    if [[ $PREFIX == cami2_toy_${BENCHMARK_DATA_SET##*_} ]]; then
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

    sed -r "s|_SAMPLEID_|gsa_${PREFIX}_${SAMPLE_NR}|g" $OUT_DIR/vamb_bins_amber.tsv > $OUT_DIR/tmp_${SAMPLE} && mv $OUT_DIR/tmp_${SAMPLE} $OUT_DIR/vamb_bins_amber.tsv

    ./get_amber_taxid_v2.py "False" \
                            ${GSA_TAX_MAP} \
                            ${TAX_FILE} \
                            "$OUT_DIR/vamb_bins_amber.tsv"
  fi
  exit 0
else
  # samtools index ${ASSEMBLY_MAP%.*}_name_sort.bam
  if [[ ! -f ${ASSEMBLY_MAP%.*}_name_sort.bam ]]; then
    conda activate mapping
    samtools sort -n --threads $THREADS $ASSEMBLY_MAP > ${ASSEMBLY_MAP%.*}_name_sort.bam
    conda deactivate
  fi
  # Run VAMB
  conda activate vamb
  # -t 64 for small samples
  vamb --fasta $ASSEMBLY \
       --bamfiles ${ASSEMBLY_MAP%.*}_name_sort.bam \
       --outdir $OUT_DIR \
       -p ${THREADS} \
       --minfasta 200000 \
       -z 0.95 \
       -s 30 \
       -o "C" \
       -t 64
fi

# Move bin fastas and rename.
for i in $OUT_DIR/bins/*.fna; do
  mv $i ${i%.*}.fasta
done

# Exit if no gsa data is processed.
if [[ $PREFIX == 'terabase_soil' ]]
then
  exit 0
fi

if $SNAKEMAKE; then
 exit 0
fi

if [[ $DATA_SET == 'CAMI2_toy_gut' || $DATA_SET == 'CAMI1_toy_high' || $DATA_SET == 'CAMI1_toy_high_pool' ]]; then
  # Create amber input file
  conda activate amber_v2
  python3 convert_fasta_bins_to_biobox_format.py $OUT_DIR/bins/*.fasta -o $OUT_DIR/vamb_bins_amber.tsv
  # Name cleanup for cami2 toy gut samples
  if [[ $PREFIX == cami2_toy_${BENCHMARK_DATA_SET##*_} ]]
  then
    SAMPLE_NR=$(echo $SAMPLE | cut -d "_" -f 3-4)
  elif [[ $PREFIX == cami1_toy_${BENCHMARK_DATA_SET##*_} ]]; then
    SAMPLE_NR="sample_$(echo $SAMPLE | cut -d "0" -f 3)"
  else
    SAMPLE_NR=$SAMPLE
  fi
  sed -r "s|_SAMPLEID_|gsa_${PREFIX}_${SAMPLE_NR}|g" $OUT_DIR/vamb_bins_amber.tsv > $OUT_DIR/tmp_${SAMPLE} && mv $OUT_DIR/tmp_${SAMPLE} $OUT_DIR/vamb_bins_amber.tsv
  
  ./get_amber_taxid_v2.py "False" \
                         "cami_data/${BENCHMARK_DATA_SET}/short_read/${SAMPLE}/contigs/gsa_mapping_amber_hr_tax.tsv"\
                         "cami_data/${BENCHMARK_DATA_SET}/short_read/taxonomic_profile_${SAMPLE_NR##*_}.txt" \
                         "$OUT_DIR/vamb_bins_amber.tsv"
fi

