#!/usr/bin/env
#SBATCH --partition=bgmp        ### Partition (like a queue in PBS)
#SBATCH --job-name=deduper     ### Job Name
#SBATCH --output=deduper_%j.out       ### File in which to store job output
#SBATCH --error=deduper_j%.err       ### File in which to store job error messages
#SBATCH --nodes=1               ### Number of nodes needed for the job
#SBATCH --ntasks-per-node=1     ### Number of tasks to be launched per Node
#SBATCH --cpus-per-task=7       ### Number of tasks to be launched
#SBATCH --account=bgmp          ### Account used for job submission
#SBATCH --time=00:02:00

conda activate bgmp_py37

/usr/bin/time -v ./deduper.py -f Dataset1.sam -p False -u STL96.txt

