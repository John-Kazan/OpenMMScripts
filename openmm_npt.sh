#!/usr/bin/env bash

# __author__ = 'Can Kazan'
# __version__ = 1.0

code_name=${0##*/}

echo "--${code_name} start: $(date)"

openmm_production_npt_gpu () {
echo "--${FUNCNAME[0]} start: $(date)"
source ~/.bashrc
module load cuda-12.4.0-ld
conda activate OpenMM
python -m openmm.testInstallation

openmm_npt.py -pdb ${pdb_file}

echo "--${FUNCNAME[0]} end: $(date)"
}

while [[ "$#" > 0 ]]; do case $1 in
-pdb|--input_pdb_file) pdb_file="$2"; shift;;
*) echo "Unknown parameter passed: $1"; exit 1;;
esac; shift; done

pdb_file=${pdb_file}; echo "Protein file: ${pdb_file}"
pdb_file_name=${pdb_file/.pdb/}

initial_directory=$(pwd)
work_directory=${initial_directory}/${code_name/.sh/}_run
mkdir -pv ${work_directory}
cp -v ${pdb_file} ${work_directory}/
cd ${work_directory}; pwd

(openmm_production_npt_gpu)

echo "--${code_name} end: $(date)"