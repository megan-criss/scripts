#!/bin/bash
#SBATCH --partition=bgmp        ### Partition (like a queue in PBS)
#SBATCH --job-name=demultiplex     ### Job Name
#SBATCH --output=demultiplex_%j.out       ### File in which to store job output
#SBATCH --error=demultiplex_%j.err       ### File in which to store job error messages
#SBATCH --time=10:00:00       ### Wall clock time limit in Days-HH:MM:SS
#SBATCH --nodes=1               ### Number of nodes needed for the job
#SBATCH --ntasks-per-node=1     ### Number of tasks to be launched per Node
#SBATCH --cpus-per-task=1       ### Number of tasks to be launched
#SBATCH --account=bgmp          ### Account used for job submission

conda activate bgmp_py37

./demultiplexing_algorithm.py -idx /projects/bgmp/shared/2017_sequencing/indexes.txt -r1 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz -i1 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz -i2 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz -r2 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz -q 25