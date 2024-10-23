#!/bin/bash
#SBATCH --job-name=nextflow
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mem=10G
#SBATCH --mail-user=emily.trybulec@uconn.edu

module load nextflow/23.10.1
module load singularity/vcell-3.10.0

export TMPDIR=$PWD/tmp
nextflow pull emilytrybulec/repeat_curation
nextflow run emilytrybulec/repeat_curation -params-file my_params.yaml -c my_config -profile singularity,xanadu -r main
