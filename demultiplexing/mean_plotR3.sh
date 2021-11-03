#!/usr/bin/bash
#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH --output=R3_%j.out       ### File in which to store job output
#SBATCH --error=R3_%j.err       ### File in which to store job error messages
#SBATCH --job-name=R3                 
#SBATCH --ntasks-per-node=1                                                  
#SBATCH --time=0-05:00:00
#SBATCH --cpus-per-task=1       ### Number of tasks to be launched
#SBATCH --nodes=1               ### Number of nodes needed for the job

conda activate bgmp_py37

./part1.py -l 101 -p r3_avg_phredperposition -f /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz
