#!/bin/bash

#SBATCH -A JNATHAN-SL2-CPU
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -p cclake
#SBATCH -D /home/nw416/rds/rds-jan_1-tpuFdqHBAEk/gps_ont
#SBATCH -o ont_qc.log
#SBATCH -c 20
#SBATCH -t 06:00:00
#SBATCH --mem=30G
#SBATCH -J ont_qc

# Initialize Conda for script usage
source "/home/nw416/miniforge3/etc/profile.d/conda.sh"
conda activate ont

DATA_DIR=/home/nw416/rds/rds-jan_1-tpuFdqHBAEk/gps_ont
RAW_FASTQ=$DATA_DIR/284LHC_1_GPS-S1S6-L1L6.fastq.gz

# QC of raw data
# ---------------------------------------
mkdir -p fastqc_raw
mkdir -p nanoplot_raw

#fastqc -t 20 $RAW_FASTQ -o fastqc_raw 2>&1 | tee fastqc_raw.log
NanoPlot --fastq $RAW_FASTQ -o nanoplot_raw --threads 20 2>&1 | tee nanoplot_raw.log
