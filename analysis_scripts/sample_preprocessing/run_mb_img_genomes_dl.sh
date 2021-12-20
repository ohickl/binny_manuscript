#!/bin/bash -l
#SBATCH -J download_metabat_img_genomes
#SBATCH --mail-type=fail
#SBATCH --mail-user=oskar.hickl@uni.lu
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --time=1-00:00:00
#SBATCH -p batch
#SBATCH --qos=normal
#SBATCH -o "slurm_logs/%A-%a-%x-%j.out"

echo "== Starting run at $(date)"
echo "== Job ID: ${SLURM_JOBID}"
echo "== Node list: ${SLURM_NODELIST}"
echo "== Submit dir. : ${SLURM_SUBMIT_DIR}"


curl 'https://signon.jgi.doe.gov/signon/create' --data-urlencode 'login=oskar.hickl@uni.lu' --data-urlencode 'password=eivqb2ACershxLP' -c img_data/cookies > /dev/null
curl -L https://genome.jgi.doe.gov/portal/ext-api/downloads/img/83642-1/tar -b img_data/cookies --output img_data/mb_img_genomes.tar
rm img_data/cookies
