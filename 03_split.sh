#!/bin/bash

#SBATCH -A JNATHAN-SL2-CPU
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -p cclake
#SBATCH -D /home/nw416/rds/rds-jan_1-tpuFdqHBAEk/gps_ont
#SBATCH -o split_%j.log
#SBATCH -c 8
#SBATCH -t 01:30:00
#SBATCH -J split
#SBATCH --mem=16G

# Initialize Conda for script usage
source "/home/nw416/miniforge3/etc/profile.d/conda.sh"
conda activate ont

DATA_DIR="/home/nw416/rds/rds-jan_1-tpuFdqHBAEk/gps_ont"

echo "Starting job $SLURM_JOB_ID at $(date)" > split.log

# 24 bp barcode library splitting
#-----------------------------------
echo "Splitting 24 bp barcode library and ORF sequences with cutadapt..." >> split.log

TRIMMED_FASTQ="$DATA_DIR/cutadapt_24.fastq"
UNSPLIT_FASTQ=$DATA_DIR/trimmed_unsplit_24.fastq
BARCODES=$DATA_DIR/barcode_24.txt
ORFS=$DATA_DIR/orf_24.fastq.gz
LEN_DIST=$DATA_DIR/barcode_24_length_distribution.txt

cutadapt \
	-a AACCCAGCTTTCTTGTACAAAGTGGTTTAATGAGTTTAAACCTCGAGCGTAACTATAACGGTCCTAAGGTAGCGAACCAGTAGGTCCACTATGAGT \
	-e 0.05 \
	-O 96 \
	-j 8 \
	--untrimmed-output $UNSPLIT_FASTQ \
	--rest-file $BARCODES \
	-o $ORFS \
	$TRIMMED_FASTQ >> split.log

echo "Getting length distribution of the 24 bp barcode sequences" >> split.log
awk '{print length($1)}' $BARCODES | sort | uniq -c > $LEN_DIST


# 30 bp barcode library splitting
#-----------------------------------
echo "Splitting 30 bp barcode library and ORF sequences with cutadapt..." >> split.log

TRIMMED_FASTQ="$DATA_DIR/cutadapt_30.fastq"
UNSPLIT_FASTQ=$DATA_DIR/trimmed_unsplit_30.fastq
BARCODES=$DATA_DIR/barcode_30.txt
ORFS=$DATA_DIR/orf_30.fastq.gz
LEN_DIST=$DATA_DIR/barcode_30_length_distribution.txt

cutadapt \
	-a CCAACTTTCTTGTACAAAGTGGTTTAATGAGTTTAAACCTCGAGCGTAACTATAACGGTCCTAAGGTAGCGAACCAGTAGGTCCACTATGAGT \
	-e 0.05 \
	-O 96 \
	-j 8 \
	--untrimmed-output $UNSPLIT_FASTQ \
	--rest-file $BARCODES \
	-o $ORFS \
	$TRIMMED_FASTQ >> split.log

echo "Getting length distribution of the 30 bp barcode sequences" >> split.log
awk '{print length($1)}' $BARCODES | sort | uniq -c > $LEN_DIST


echo "Finished job $SLURM_JOB_ID at $(date)" >> split.log