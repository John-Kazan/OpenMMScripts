#!/usr/bin/env bash
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=6
#SBATCH --partition=general
#SBATCH --qos=grp_sozkan
#SBATCH --gres=gpu:1
#SBATCH --output=slurm.%j.out
#SBATCH --error=slurm.%j.err
#SBATCH --mail-type=ALL

export PATH="/scratch/ikazan/ztest:${PATH}"
openmm_npt.sh -pdb 1btl.pdb

seff $SLURM_JOB_ID