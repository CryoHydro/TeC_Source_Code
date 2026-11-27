#!/bin/bash
#SBATCH --job-name=KYZ_OM_test
#SBATCH -o /nfs/scistore18/pelligrp/ajoubert/Jobs/logs/TC_distributed/TC_Dist_KYZ_%j.out
#SBATCH -e /nfs/scistore18/pelligrp/ajoubert/Jobs/logs/TC_distributed/TC_Dist_KYZ_%j.err
#SBATCH --chdir=/nfs/scistore18/pelligrp/ajoubert/TeC_Source_Code   # repository: for the commit info and to start the launcher
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

# Settings read by Launcher_Dev_Distributed.m (they replace the values set in the launcher and the config)
export TC_REPO_DIR=/nfs/scistore18/pelligrp/ajoubert/TeC_Source_Code   # T&C code (replaces CONFIG_vals.REPO_Dir)
export TC_ROOT_DIR=$TC_REPO_DIR/Case_study                              # Working directory: folder of the site folders
export TC_SITE_NAME=Kyzylsu                                             # Site folder: $TC_ROOT_DIR/$TC_SITE_NAME
export TC_OUTPUT_DIR=/nfs/scistore18/pelligrp/ajoubert/TC_outputs/Kyzylsu # Outputs: $TC_OUTPUT_DIR/<Run_Folder> (replaces CONFIG_vals.Output_Dir)
# export TC_CONFIG=Site_Config_Run03.m   # Optional: another config in $TC_ROOT_DIR/$TC_SITE_NAME/Parameters/
#                                          (default Site_Config_TC.m)

echo "Job $SLURM_JOB_ID on $(hostname), $SLURM_CPUS_PER_TASK cores, started $(date)"
echo "Site: $TC_ROOT_DIR/$TC_SITE_NAME, config: ${TC_CONFIG:-Site_Config_TC.m}, outputs: $TC_OUTPUT_DIR"

ref=$(sed -n 's/^ref: //p' .git/HEAD)
echo "Code: ${ref#refs/heads/}, commit $(cat .git/$ref 2>/dev/null || grep " $ref\$" .git/packed-refs | cut -c1-40)"

module load matlab/R2023b
srun matlab -batch "delete(gcp('nocreate')); run('Launcher_Dev_Distributed.m')"
rc=$?   # exit code of matlab (read before any other command replaces it)
echo "Finished $(date), exit code $rc"
