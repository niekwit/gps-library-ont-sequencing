#!/bin/bash

#!/bin/bash

#SBATCH -A JNATHAN-SL2-CPU
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -p cclake
#SBATCH -D /home/nw416/rds/rds-jan_1-tpuFdqHBAEk/gps_ont
#SBATCH -o barcode_lengths.log
#SBATCH -c 2
#SBATCH -t 0:30:00
#SBATCH -J barcode_lengths
#SBATCH --mem=8G

# Initialize Conda for script usage
source "/home/nw416/miniforge3/etc/profile.d/conda.sh"
conda activate ont

DATA_DIR=/home/nw416/rds/rds-jan_1-tpuFdqHBAEk/gps_ont
BARCODES=$DATA_DIR/barcodes.txt

echo "Starting job $SLURM_JOB_ID at $(date)" >> barcode_lengths.log
awk '{print length($1)}' $BARCODES | sort | uniq -c | sed 's/^[[:space:]]*//' > $DATA_DIR/barcode_length_distribution.txt
echo "Finished job at $(date)" >> barcode_lengths.log
