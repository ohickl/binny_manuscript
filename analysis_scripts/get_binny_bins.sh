#!/bin/bash -l
#SBATCH -J get_binny_bins
#SBATCH --mail-type=end,fail
#SBATCH --mail-user=oskar.hickl@uni.lu
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --time=1-00:00:00
#SBATCH -p batch
#SBATCH --qos=qos-batch

echo "== Starting run at $(date)"
echo "== Job ID: ${SLURM_JOBID}"
echo "== Node list: ${SLURM_NODELIST}"
echo "== Submit dir. : ${SLURM_SUBMIT_DIR}"

ASSEMBLY="binny_out/2017.12.04_18.45.54_sample_0/assembly.fa"
BIN_DIR="binny_out/2017.12.04_18.45.54_sample_0/bins"
BINNY_OUT="binny_out/2017.12.04_18.45.54_sample_0/contigs2clusters.10.4.tsv"


while IFS=$'\t' read -r col1 col2
do
    sed -n '/>$col1/,/>S0C/p' $ASSEMBLY | sed -e '$d' >> $BIN_DIR/$col2.fa
done  < $BINNY_OUT

date
