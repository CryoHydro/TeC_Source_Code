#!/bin/bash
#SBATCH --job-name=KYZ_OM_test
#SBATCH -o /nfs/scistore18/pelligrp/ajoubert/Jobs/logs/TC_distributed/TC_Dist_KYZ_%j.out
#SBATCH -e /nfs/scistore18/pelligrp/ajoubert/Jobs/logs/TC_distributed/TC_Dist_KYZ_%j.err
#SBATCH --chdir=/nfs/scistore18/pelligrp/ajoubert/TeC_Source_Code   # launcher paths are relative to this
#SBATCH --mail-user=achille.jouberton@ist.ac.at
#SBATCH --mail-type=BEGIN,END,FAIL,TIME_LIMIT_90
#SBATCH --mem=85G
#SBATCH --time=10-00:00
#SBATCH -c 40
#SBATCH --constraint=matlab
#SBATCH --exclude=eta[406]
#SBATCH --no-requeue
#SBATCH --export=NONE
unset SLURM_EXPORT_ENV

echo "Job $SLURM_JOB_ID on $(hostname), $SLURM_CPUS_PER_TASK cores, started $(date)"

ref=$(sed -n 's/^ref: //p' .git/HEAD)
echo "Code: ${ref#refs/heads/}, commit $(cat .git/$ref 2>/dev/null || grep " $ref\$" .git/packed-refs | cut -c1-40)"

module load matlab/R2023b
srun matlab -batch "delete(gcp('nocreate')); run('Launcher_Dev_Distributed.m')"
echo "Finished $(date), exit code $?"