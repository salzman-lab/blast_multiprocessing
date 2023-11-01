#!/bin/bash
#Single commented out SBATCH are active, more than that is inactive
#SBATCH -p horence,normal
#SBATCH --time=10:00:00                     # how much time to run
#SBATCH --mem=20000                          # how much mem in MBs
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --job-name=blast                 # name the job jupyter host
#SBATCH --output=blast.out                 # slurm.out file
#SBATCH --error=blast.err                  # slurm.err file


source /oak/stanford/groups/horence/george/blast_internal_parallelization/environment/bin/activate

python3 /oak/stanford/groups/horence/george/blast_internal_parallelization/blast_internal_parallelization.py $1 $2
