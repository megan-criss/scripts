#!/usr/bin/bash
#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH --output=mean_scores_%j.out       ### File in which to store job output
#SBATCH --error=mean_scores_%j.err       ### File in which to store job error messages
#SBATCH --job-name=algorithmpt1                                                              #SBATCH --output=mean+plot_%j.out                                                            #SBATCH --error=mean_plot_%j.err                                                             #SBATCH --nodes=1
#SBATCH --ntasks-per-node=1                                                                  #SBATCH --cpus-per-task=1
#SBATCH --time=0-05:00:00 
#SBATCH --cpus-per-task=1       ### Number of tasks to be launched
#SBATCH --nodes=1               ### Number of nodes needed for the job


conda activate bgmp_py37

./plot_script.py -l 101 -p r1_avg_phredperposition -f /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz 
