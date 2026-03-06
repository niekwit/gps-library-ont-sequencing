#!/bin/bash

#SBATCH -A JNATHAN-SL2-CPU
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -p cclake
#SBATCH -D /home/nw416/rds/rds-jan_1-tpuFdqHBAEk/gps_ont
#SBATCH -o barcode_fq_%j.log
#SBATCH -c 4
#SBATCH -t 02:00:00
#SBATCH -J barcode_fq
#SBATCH --mem=8G

# Initialize Conda for script usage
source "/home/nw416/miniforge3/etc/profile.d/conda.sh"
conda activate ont

DATA_DIR=/home/nw416/rds/rds-jan_1-tpuFdqHBAEk/gps_ont

echo "Starting job $SLURM_JOB_ID at $(date)" > barcode_fq.log

# 24 bp barcode library
# ---------------------------
TRIMMED_FASTQ=$DATA_DIR/cutadapt_24.fastq
BARCODES=$DATA_DIR/barcode_24.txt
BARCODE_FQ=$DATA_DIR/barcode_24.fastq

python extract_barcodes.py $BARCODES $TRIMMED_FASTQ $BARCODE_FQ

# 30 bp barcode library
# ---------------------------
TRIMMED_FASTQ=$DATA_DIR/cutadapt_30.fastq
BARCODES=$DATA_DIR/barcode_30.txt
BARCODE_FQ=$DATA_DIR/barcode_30.fastq

python extract_barcodes.py $BARCODES $TRIMMED_FASTQ $BARCODE_FQ

echo "Finished job at $(date)" >> barcode_fq.log