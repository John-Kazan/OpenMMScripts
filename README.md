## OpenMM NPT Simulation Scripts

This repository provides scripts for running constant pressure, constant temperature (NPT) molecular dynamics simulations using OpenMM. 

### Installation

1. **Install OpenMM with conda:**

    ```bash
    conda create -y -n OpenMM -c conda-forge &&
    conda activate OpenMM && 
    conda install -y openmm
    ```

2. **Clone this repository:**

   ```bash
   git clone <repository_url>
   ```

### Running Simulations

#### Using `openmm_npt.py`

This script runs the NPT simulation directly using OpenMM. It utilizes GPU (CUDA) for acceleration.

**Usage:**

```bash
python openmm_npt.py -pdb 1btl.pdb
```

Replace `1btl.pdb` with the path to your PDB file.

#### Using `openmm_npt_gpu.sh` (for PHX cluster)

This wrapper script automatically loads the correct CUDA version and OpenMM environment on the PHX cluster.

**Usage:**

```bash
./openmm_npt_gpu.sh -pdb 1btl.pdb
```

Replace `1btl.pdb` with your PDB file or use the one provided for testing.

### Submitting to Scheduler

To submit the simulation to a scheduler (e.g., SLURM), use the `submit_sbatch.sh` script. 

**Before submission:**

1. **Edit `submit_sbatch.sh`:**  
   - Modify the `export PATH="/scratch/ikazan/ztest:${PATH}"` line to reflect your own path.

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