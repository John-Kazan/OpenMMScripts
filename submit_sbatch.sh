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

(
export PATH="/path/to/your/cloned/repository/OpenMMScripts:${PATH}"
openmm_npt.sh -pdb ./pdb/1btl.pdb
)

sbatch submit_sbatch.sh
