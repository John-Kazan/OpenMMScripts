#!/usr/bin/env python
'''
openmm npt 

'''

import sys
import os
import glob
import re
import argparse
import time
from pathlib import Path

from pdbfixer.pdbfixer import PDBFixer
from openmm import app
import openmm as mm
from openmm import unit

__author__ = 'I. Can Kazan'
__version__ = 1.0

parser = argparse.ArgumentParser(
    description='openmm npt gpu version: {}'.format(__version__),
    epilog='Made by: {}'.format(__author__)
)

file_in_out_group = parser.add_argument_group('Input / output files')
file_in_out_group.add_argument(
    '-pdb', '--input_pdb_file',
    help='Input PDB file', type=str, nargs='?', default='N/A')

if len(sys.argv) == 1:
    parser.print_help()
    exit()

args = parser.parse_args()
print('Arguments: {}'.format(args))

def wrap_function():
    separator = '-' * 76
    def decorator(func):
        def wrapper(*args, **kwargs):
            print(separator)
            print(f'Calling function: {func.__name__}')
            start_time = time.time()
            result = func(*args, **kwargs)
            end_time = time.time()
            print(f'Function {func.__name__} completed')
            print(f'Time taken: {end_time - start_time:.4f} seconds')
            print(separator)
            return result
        return wrapper
    return decorator

def find_last_file(pattern):    
    # Use glob to find all files matching the pattern
    files = glob.glob(pattern)
    
    if not files:
        print("No files found.")
        return None
    
    # Extract the numeric part from filenames
    file_numbers = []
    for file in files:
        match = re.search(r'_(\d+)\.xml$', file)
        if match:
            file_numbers.append((int(match.group(1)), file))
    
    # Sort the files by the numeric part
    if not file_numbers:
        print("No files with the expected numeric pattern found.")
        return None
    
    file_numbers.sort()
    
    # Return the file with the highest number
    last_file = file_numbers[-1][1]
    return last_file

@wrap_function()
def fix_pdb(input_pdb_file):
    output_pdb_file = 'fixed.pdb'

    if Path(output_pdb_file).is_file():
        print(f'Output PDB file {output_pdb_file} exists')
        return
    else:
        print(f'Output PDB file {output_pdb_file} does not exist !')

    if Path(input_pdb_file).is_file():
        print(f'{input_pdb_file} exists')
    else:
        print(f'{input_pdb_file} does not exist !')
        exit()
    print('Fixing PDB...')
    pdb = PDBFixer(filename=input_pdb_file)
    pdb.findMissingResidues()
    pdb.findNonstandardResidues()
    pdb.replaceNonstandardResidues()
    pdb.removeHeterogens(True)
    pdb.findMissingAtoms()
    pdb.addMissingAtoms()
    pdb.addMissingHydrogens(7.0)
    topology = pdb.topology
    positions = pdb.positions
    app.PDBFile.writeFile(topology, positions, open(output_pdb_file, 'w'))
    print(f'PDB fixing completed')

@wrap_function()
def prepare_toppos(input_pdb_file):
    output_pdb_file = 'system.pdb'

    if Path(output_pdb_file).is_file():
        print(f'Output PDB file {output_pdb_file} exists')
        return
    else:
        print(f'Output PDB file {output_pdb_file} does not exist !')

    if Path(input_pdb_file).is_file():
        print(f'{input_pdb_file} exists')
    else:
        print(f'{input_pdb_file} does not exist !')
        exit()
    pdb = app.PDBFile(input_pdb_file)
    topology = pdb.getTopology()
    positions = pdb.getPositions()

    force_field = 'amber14-all.xml'
    water_model = 'amber14/tip3p.xml'
    water_box_padding = 16
    water_box_shape = 'cube'
    ions = {
        'positive':'Na+',
        'negative': 'Cl-'
    }
    ion_concentration = 0

    print('Loading hydrogen definitions...')
    app.Modeller.loadHydrogenDefinitions('glycam-hydrogens.xml')

    print('Creating forcefield...')
    forcefield = app.ForceField(force_field, water_model)

    print('Creating modeller...')
    modeller = app.Modeller(topology, positions)

    print('Adding hydrogens...')
    modeller.addHydrogens(forcefield)

    print('Adding solvent...')
    modeller.addSolvent(
        forcefield=forcefield,
        model=os.path.basename(os.path.realpath(water_model)).replace('.xml', ''),
        padding=water_box_padding*unit.angstroms.conversion_factor_to(unit.nanometers),
        boxShape=water_box_shape,
        positiveIon=ions['positive'],
        negativeIon=ions['negative'],
        ionicStrength=ion_concentration*unit.molar,
    )

    topology = modeller.topology
    positions = modeller.positions
    app.PDBFile.writeFile(topology, positions, open(output_pdb_file, 'w'))

    print(f'Topology and positions preparation completed')

def create_system(topology):
    force_field = 'amber14-all.xml'
    water_model = 'amber14/tip3p.xml'

    forcefield = app.ForceField(force_field, water_model)

    nonbonded_method = app.PME
    nonbonded_cutoff = 12.0*unit.angstroms.conversion_factor_to(unit.nanometers)
    constraints = app.HBonds
    rigid_water = True
    ewald_error_tolerance = 0.0005 # suggested 0.0005 i want 0.000001

    print('Creating system...')
    system = forcefield.createSystem(
        topology,
        nonbondedMethod=nonbonded_method,
        nonbondedCutoff=nonbonded_cutoff,
        constraints=constraints,
        rigidWater=rigid_water,
        ewaldErrorTolerance=ewald_error_tolerance
    )
    return system

def manage_restraints(simulation, restraint_mask, **kwargs):
    restraint_force_constant = kwargs.get(
        'restraint_force_constant', 1000.0 * unit.kilojoules_per_mole / unit.nanometer**2
    )
    action = kwargs.get('action', 'apply')

    topology = simulation.topology
    system = simulation.system

    state = simulation.context.getState(getPositions=True)
    positions = state.getPositions()

    if action == 'apply' and restraint_mask:
        print('Applying restraints...')
        force = mm.CustomExternalForce('k*periodicdistance(x, y, z, x0, y0, z0)^2')
        force.addGlobalParameter('k', restraint_force_constant)
        force.addPerParticleParameter('x0')
        force.addPerParticleParameter('y0')
        force.addPerParticleParameter('z0')
        for atom in topology.atoms():
            if atom.residue.name not in restraint_mask:
                # print(atom, atom.index, positions[atom.index])
                force.addParticle(atom.index, positions[atom.index])
        system.addForce(force)
    elif action == 'remove':
        print('Removing restraints...')
        forces_to_remove = []
        for i, force in enumerate(system.getForces()):
            if isinstance(force, mm.CustomExternalForce):
                if force.getEnergyFunction() == 'k*periodicdistance(x, y, z, x0, y0, z0)^2':
                    forces_to_remove.append(i)
        for i in reversed(forces_to_remove):
            system.removeForce(i)
    
    for i, force in enumerate(system.getForces()):
        print(i, force)

def platform_select(platform_name, **kwargs):
    num_cores = kwargs.get('num_cores', 6)
    print('Getting platform...')
    platform = mm.Platform.getPlatformByName(platform_name)
    if platform_name == 'CPU':
        platform.setPropertyDefaultValue('Threads', str(num_cores))
    print(f'Platform: {platform}')
    return platform

def add_reporters(simulation, nsteps, output_prefix, **kwargs):
    interval = 5000
    checkpoint_interval = kwargs.get('checkpoint_interval', interval)
    dcd_interval = kwargs.get('dcd_interval', interval)
    pdb_interval = kwargs.get('pdb_interval', 500000)
    state_interval = kwargs.get('state_interval', interval)

    if nsteps <= pdb_interval:
        pdb_interval = interval

    print('Adding Reporters...')
    print('Adding CheckpointReporter...')
    checkpoint_file = f'{output_prefix}.chk'
    simulation.reporters.append(
        app.CheckpointReporter(
            checkpoint_file,
            checkpoint_interval
        )
    )
    print('Adding DCDReporter...')
    dcd_file = f'{output_prefix}.dcd'
    simulation.reporters.append(
        app.DCDReporter(
            dcd_file,
            dcd_interval
        )
    )
    print('Adding PDBReporter...')
    pdb_file = f'{output_prefix}.pdb'
    simulation.reporters.append(
        app.PDBReporter(
            pdb_file,
            pdb_interval,
            enforcePeriodicBox=True
        )
    )
    print('Adding StateDataReporter...')
    state_file = f'{output_prefix}.csv'
    simulation.reporters.append(
        app.StateDataReporter(
            state_file,
            state_interval,
            step=True, time=True, potentialEnergy=True, kineticEnergy=True, totalEnergy=True, temperature=True,
            volume=True, density=True, progress=True, remainingTime=True, speed=True, elapsedTime=True,
            separator=',', systemMass=None, totalSteps=nsteps
        )
    )

def simulation_prep(**kwargs):
    platform_name = kwargs.get('platform_name', 'OpenCL')
    topology = kwargs.get('topology')
    system = kwargs.get('system')
    integrator = kwargs.get('integrator')
    input_state_file = kwargs.get('input_state_file')
    pdb = kwargs.get('pdb')
    restraint_mask = kwargs.get('restraint_mask')
    nsteps = kwargs.get('nsteps')
    output_prefix = kwargs.get('output_prefix')
    
    print('Creating simulation...')
    platform = platform_select(platform_name)
    simulation = app.Simulation(topology, system, integrator, platform)

    if Path(input_state_file).is_file():
        print(f'State file {input_state_file} exists')
        simulation.loadState(input_state_file)
    else:
        print(f'State file {input_state_file} does not exist !')
        positions = pdb.getPositions()
        simulation.context.setPositions(positions)

    manage_restraints(simulation, restraint_mask)

    add_reporters(simulation, nsteps, output_prefix)
    return simulation

@wrap_function()
def simulation_run(**kwargs):
    input_pdb_file = kwargs.get('input_pdb_file')
    input_state_file = kwargs.get('input_state_file')
    output_prefix = kwargs.get('output_prefix')
    mode = kwargs.get('mode', 'npt')
    platform_name = kwargs.get('platform_name', 'OpenCL')
    temperature = kwargs.get('temperature', 300 * unit.kelvin)
    pressure = kwargs.get('pressure', 1 * unit.bar)
    nsteps = kwargs.get('nsteps', 50000)
    restraint_mask = kwargs.get('restraint_mask', '')
    start_temp = kwargs.get('start_temp', 10 * unit.kelvin)
    end_temp = kwargs.get('end_temp',temperature)

    output_state_file = f'{output_prefix}.xml'

    if Path(output_state_file).is_file():
        print(f'Output state file {output_state_file} exists')
        return
    else:
        print(f'Output state file {output_state_file} does not exist !')

    if Path(input_pdb_file).is_file():
        print(f'PDB file {input_pdb_file} exists')
    else:
        print(f'PDB file {input_pdb_file} does not exist !')
        exit()

    pdb = app.PDBFile(input_pdb_file)
    topology = pdb.getTopology()
    system = create_system(topology)

    if mode == 'npt':
        print('Adding barostat...')
        barostat_frequency = 25
        barostat = mm.MonteCarloBarostat(pressure, temperature, barostat_frequency)
        system.addForce(barostat)
        print('Adding thermostat...')
        friction_coefficient = 2 / unit.picosecond
        timestep = 0.002 * unit.picoseconds
        integrator = mm.LangevinIntegrator(temperature, friction_coefficient, timestep)
    elif mode == 'heat':
        print('Adding thermostat...')
        friction_coefficient = 2 / unit.picosecond
        timestep = 0.002 * unit.picoseconds
        integrator = mm.LangevinIntegrator(temperature, friction_coefficient, timestep)
    elif mode == 'min':
        print('Creating integrator...')
        timestep = 0.002*unit.picoseconds
        integrator = mm.VerletIntegrator(timestep)

    print(f'Creating {mode} simulation...')
    simulation_parameters = {
        'platform_name': platform_name,
        'topology': topology,
        'system': system,
        'integrator': integrator,
        'input_state_file': input_state_file,
        'pdb': pdb,
        'restraint_mask': restraint_mask,
        'nsteps': nsteps,
        'output_prefix': output_prefix
    }
    simulation = simulation_prep(**simulation_parameters)

    print(f'Running {mode} simulation...')

    simulation.currentStep = 0
    if mode == 'min':
        simulation.minimizeEnergy(maxIterations=nsteps)
    elif mode == 'heat':
        current_temp = start_temp
        for step in range(nsteps+1):
            simulation.integrator.setTemperature(current_temp)
            simulation.step(1)
            current_temp = start_temp + ((end_temp - start_temp) / nsteps) * step
    elif mode == 'npt':
        simulation.step(nsteps)
    
    simulation.saveState(f'{output_prefix}.xml')

    with open(f'{output_prefix}_last.pdb', 'w') as output:
        state = simulation.context.getState(getPositions=True)
        app.PDBFile.writeFile(topology, state.getPositions(), output, keepIds=True)

    print(f'{mode} simulation completed')

def main():
    input_pdb_file = args.input_pdb_file

    pdb_file_name = os.path.basename(os.path.realpath(input_pdb_file)).replace('.pdb', '')
    print(pdb_file_name)

    # Fix pdb
    fix_pdb(input_pdb_file)

    # prepare system
    prepare_toppos('fixed.pdb')

    # Solvent minimization
    print('Solvent minimization begin...')
    solvent_minimization_parameters = {
        'input_pdb_file': 'system.pdb',
        'input_state_file': '',
        'output_prefix': f'{pdb_file_name}_solvent_minimized',
        'mode': 'min',
        'nsteps': 50000,
        'restraint_mask': ['HOH','NA','CL']
    }
    simulation_run(**solvent_minimization_parameters)
    print('Solvent minimization end...')
    
    # Solution minimization
    print('Solution minimization begin...')
    solution_minimization_parameters = {
        'input_pdb_file': 'system.pdb',
        'input_state_file': f'{pdb_file_name}_solvent_minimized.xml',
        'output_prefix': f'{pdb_file_name}_solution_minimized',
        'mode': 'min',
        'nsteps': 50000,
    }
    simulation_run(**solution_minimization_parameters)
    print('Solution minimization end...')

    # Heatup with restraints
    print('Heatup begin...')
    heatup_parameters = {
        'input_pdb_file': 'system.pdb',
        'input_state_file': f'{pdb_file_name}_solution_minimized.xml',
        'output_prefix': f'{pdb_file_name}_heated',
        'mode': 'heat',
        'platform_name': 'CPU',
        'nsteps': 50000,
        'restraint_mask': ['HOH','NA','CL']
    }
    simulation_run(**heatup_parameters)
    print('Heatup end...')

    # Production runs (CPU)
    print('Production NPT (CPU) begin...')
    npt_cpu_parameters = {
        'input_pdb_file': 'system.pdb',
        'input_state_file': f'{pdb_file_name}_heated.xml',
        'output_prefix': f'{pdb_file_name}_npt_cpu',
        'platform_name': 'CPU',
        'nsteps': 50000,
    }
    simulation_run(**npt_cpu_parameters)
    print('Production NPT (CPU) end...')

    # Production runs (GPU)
    print('Production NPT (GPU) begin...')
    npt_gpu_parameters = {
        'input_pdb_file': 'system.pdb',
        'input_state_file': f'{pdb_file_name}_npt_cpu.xml',
        'output_prefix': f'{pdb_file_name}_npt_gpu_00000',
        'platform_name': 'CUDA',
        'nsteps': 50000,
    }
    simulation_run(**npt_gpu_parameters)
    print('Production NPT (GPU) end...')

    last_file = find_last_file(f'{pdb_file_name}_npt_gpu_*.xml')

    # Production runs continue (GPU)
    print('Production NPT (GPU) continue...')
    npt_gpu_parameters = {
        'input_pdb_file': 'system.pdb',
        'input_state_file': f'{last_file}',
        'output_prefix': f'{pdb_file_name}_npt_gpu_{int(last_file.replace('.xml', '').split('_')[-1])+1:05}',
        'platform_name': 'CUDA',
        'nsteps': 50000000, # 50000000 * 0.002 = 100ns
    }
    simulation_run(**npt_gpu_parameters)
    print('Production NPT (GPU) pause...')

if __name__ == '__main__':
    main()