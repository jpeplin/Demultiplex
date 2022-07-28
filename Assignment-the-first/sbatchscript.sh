#!/bin/bash

#SBATCH --partition=bgmp        ### Partition (like a queue in PBS)
#SBATCH --job-name=demuxR4      ### Job Name
#SBATCH --output=result-R4-%j.out  ### File in which to store job output
#SBATCH --error=result-R4-%j.err   ### File in which to store job error messages
#SBATCH --nodes=1               ### Number of nodes needed for the job
#SBATCH --cpus-per-task=1       ### CPUs
#SBATCH --account=bgmp          ### Account used for job submission

conda activate bgmp_py310

cd /projects/bgmp/jpeplin5/bioinfo/Bi622/Demultiplex/Assignment-the-first

/usr/bin/time -v ./bioscript.py \
    -f /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz \
    -n 363246735 -l 101 -o R4.png