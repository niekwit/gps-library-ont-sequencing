#!/bin/bash

#SBATCH -A JNATHAN-SL2-CPU
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -p cclake
#SBATCH -D /home/nw416/rds/rds-jan_1-tpuFdqHBAEk/gps_ont
#SBATCH -o ont_trim_%j.log
#SBATCH -c 12
#SBATCH -t 3:00:00
#SBATCH -J ont_trim
#SBATCH --mem=24G

# Initialize Conda for script usage
source "/home/nw416/miniforge3/etc/profile.d/conda.sh"
conda activate ont

DATA_DIR=/home/nw416/rds/rds-jan_1-tpuFdqHBAEk/gps_ont
RAW_FASTQ="$DATA_DIR/284LHC_1_GPS-S1S6-L1L6.fastq.gz"

echo "Starting job $SLURM_JOB_ID at $(date)" > trim.log

{
	zcat "$RAW_FASTQ" | \
	chopper -q 15 | \
	cutadapt -g "CAAGTTTGTACAAAAAAGCAGGCACC...AGCTGTGTAAGCGGAACTAGAGTTAC" \
		--untrimmed-output "$DATA_DIR/untrimmed_24.fastq" \
		--revcomp \
		-e 0.05 \
		-O 26 \
        -j 10 \
		-o "$DATA_DIR/cutadapt_24.fastq" - 
} &> cutadapt_24bp_${SLURM_JOB_ID}.log

{
	zcat "$RAW_FASTQ" | \
	chopper -q 15 | \
	cutadapt -g "AAACAAGTTTGTACAAAAAAGTTGGC...AGCTGTGTAAGCGGAACTAGAGTTAC" \
		--untrimmed-output "$DATA_DIR/untrimmed_30.fastq" \
		--revcomp \
		-e 0.05 \
		-O 26 \
        -j 10 \
		-o "$DATA_DIR/cutadapt_30.fastq" - 
} &> cutadapt_30bp_${SLURM_JOB_ID}.log

echo "Finished job at $(date)" >> trim.log
