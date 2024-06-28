## OpenMM NPT Simulation Scripts

This repository contains scripts for running constant pressure, constant temperature (NPT) molecular dynamics simulations using OpenMM. Some scripts are designed for use on the ASU PHX cluster, but can be adapted for other environments. 

The repository includes:

- `openmm_npt.py`: A Python script for running NPT simulations using OpenMM.
- `openmm_npt.sh`: A wrapper script for automatically loading the correct CUDA version and OpenMM environment on the PHX cluster.
- `submit_sbatch.sh`: A script for submitting simulation jobs to a scheduler (e.g., SLURM).
- `pdb/1btl.pdb`: A PDB file for testing.
- `setup.sh`: A helper script to set up the environment.

### Installation

1. **Install OpenMM with conda:**

    ```bash
    conda create -y -n OpenMM -c conda-forge &&
    conda activate OpenMM && 
    conda install -y openmm pdbfixer -c conda-forge
    ```

    make sure you have conda initialized in your environment:

    ```bash
    conda init bash
    ```

2. **Clone this repository:**

   ```bash
   git clone https://github.com/John-Kazan/OpenMMScripts.git
   ```

3. **Run the setup script:**

    ```bash
    setup.sh
    ```

### Running Simulations

#### Using `openmm_npt.py`

This script runs the NPT simulation directly using OpenMM. It utilizes GPU (CUDA) for acceleration.

**Usage:**

```bash
python openmm_npt.py -pdb ./pdb/1btl.pdb
```

Replace `pdb/1btl.pdb` with the path to your PDB file.

#### Using `openmm_npt.sh` (for PHX cluster)

This wrapper script automatically loads the correct CUDA version and OpenMM environment on the PHX cluster.

**Usage:**

```bash
openmm_npt.sh -pdb ./pdb/1btl.pdb
```

Replace `1btl.pdb` with your PDB file or use the one provided for testing.

### Submitting to Scheduler

To submit the simulation to a scheduler (e.g., SLURM), use the `submit_sbatch.sh` script. 

**Before submission:**

Edit `submit_sbatch.sh`:
- Modify the `openmm_npt.sh -pdb ./pdb/1btl.pdb` line to reflect your own PDB.

**Submission:**

```bash
sbatch submit_sbatch.sh
```

To continue the simulation from the last point simply submit the job again. It will automatially continue.

### Notes

- The simulation is set to run for 100ns. On RTX 2080 the performance is ~200ns/day.
- This repository is designed for running simulations on the PHX cluster, but can be adapted for other environments.
- The provided scripts are examples and may need modifications to suit your specific needs.
- For more detailed information about OpenMM, refer to the official documentation: [https://openmm.org/](https://openmm.org/)
