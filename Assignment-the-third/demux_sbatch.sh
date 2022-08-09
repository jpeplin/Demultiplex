#!/bin/bash

#SBATCH --partition=bgmp        ### Partition (like a queue in PBS)
#SBATCH --job-name=demux      ### Job Name
#SBATCH --output=result-%j.out  ### File in which to store job output
#SBATCH --error=result-%j.err   ### File in which to store job error messages
#SBATCH --nodes=1               ### Number of nodes needed for the job
#SBATCH --cpus-per-task=12       ### CPUs
#SBATCH --account=bgmp          ### Account used for job submission

conda activate bgmp_py310

cd /projects/bgmp/jpeplin5/bioinfo/Bi622/Demultiplex/Assignment-the-third

/usr/bin/time -v ./demux.py \
    -r1 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz \
    -r2 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz \
    -i1 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz \
    -i2 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz \
    -k /projects/bgmp/shared/2017_sequencing/indexes.txt