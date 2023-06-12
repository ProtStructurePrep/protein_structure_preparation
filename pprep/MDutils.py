import pdbfixer 
from openmm import app, Platform
import mdtraj as md
import nglview as nv
from openff.toolkit.topology import Molecule as OFFMolecule
from openmmforcefields.generators import SystemGenerator
import copy 
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from typing_extensions import Literal
from openmm import unit
import parmed
import openmm
from simtk.openmm.app import PDBReporter, DCDReporter, StateDataReporter
from pathlib import Path
import sys, time


def load_prepared_receptor(receptor_file):
    """
    Load the file of the corrected receptor by using OpenMM.
    
    Parameters
    ----------
    receptor_file: str
        path of the file of the fixed receptor
    
    Returns
    -------
    receptor: openmm.app.pdbfile.PDBFile
        OpenMM object containing the receptor file
    """
    receptor = app.PDBFile(receptor_file)
    return receptor
    
def load_prepared_ligand(rdkit_ligand):
    """
    Load the file of the corrected receptor by using OpenMM.
    
    Parameters
    ----------
    receptor_file: str
        path of the file of the fixed receptor
    
    Returns
    -------
    ligand: openff
        OpenMM object containing the receptor file
    """
    ligand = OFFMolecule.from_rdkit(rdkit_ligand, allow_undefined_stereo=True)
    return ligand

def prepare_system(receptor, ligand, solvate=False):
    """ 
    Function inpired from: https://github.com/cole-group/FEgrow/blob/master/fegrow/receptor.py
    
    Use openmm to prepare the system to simulate with both the ligand and the receptor.

    Parameters
    ----------
    receptor: openmm.app.pdbfile.PDBFile
        Prepared protein (without the ligand) to simulate
    ligand: rdkit.Chem.rdchem.Mol
        Prepared ligand to simulate
    
    Returns
    -------
    system: openmm.openmm.System
        Prepared complex system for further simulations
    modeller: openmm.app.modeller.Modeller
        Topology structure of the complex (with both the receptor and the ligand)
        (with solvernt if solvate = True)
    parmed_structure: parmed.structure.Structure
        Topology structure of the complex (with both the receptor and the ligand) 
        (without solvent)
    """
    
    
    # load the molecule into openff
    openff_mol = OFFMolecule.from_rdkit(ligand, allow_undefined_stereo=True)
    
    # Initialize a SystemGenerator using the GAFF for the ligand
    forcefield_kwargs = { 'constraints' : app.HBonds, 'rigidWater' : True, 'removeCMMotion' : False, 'hydrogenMass' : 4*unit.amu }
    system_generator = SystemGenerator(
        forcefields=['amber/ff14SB.xml', 'amber/tip3p_standard.xml'], 
        small_molecule_forcefield='openff-1.2.0', 
        forcefield_kwargs=forcefield_kwargs, 
        cache='db.json',
        molecules=openff_mol,
    )
    

    # now make a combined receptor and ligand topology
    parmed_receptor = parmed.openmm.load_topology(
        receptor.topology, xyz=receptor.positions
    )
    parmed_ligand = parmed.openmm.load_topology(
        openff_mol.to_topology().to_openmm(), xyz=openff_mol.conformers[0]
    )

    parmed_structure = parmed_receptor + parmed_ligand
    
    
    modeller = app.Modeller(parmed_structure.topology, parmed_structure.positions)
    if solvate:
        # Solvate
        modeller.addSolvent(system_generator.forcefield, model='tip3p', padding=5.0*unit.angstroms, ionicStrength=0.15*unit.molar)
        system = system_generator.create_system(modeller.topology)
    else:
        system = system_generator.create_system(modeller.topology)

    return system, modeller, parmed_structure
    
def setup_openmm_simulation(system, modeller, time_step=1*unit.femtoseconds, temperature=300*unit.kelvin, friction=1/unit.picosecond, platform_name="CPU"):
    """
    Propagate the System with Langevin dynamics and set up the openmm simulation.
    
    Parameters
    ----------
    system: openmm.openmm.System
        Prepared complex system for further simulations
    modeller: 
        Topology structure of the complex (with both the receptor and the ligand)
    time_step: openmm.unit.quantity.Quantity
        Simulaion timestep with default value set to 1 * unit.femtoseconds 
    temperature: openmm.unit.quantity.Quantity
        Simulation temperature with default value set to 300 * unit.kelvin 
    friction: openmm.unit.quantity.Quantity
        Collision rate with default value set to 1 / unit.picosecond 
    platform_name: str
         OpenMM platform name, with the 'cpu' used by default
         
    Returns
    -------
    simulation: openmm.app.simulation.Simulation
        OpenMM initial simulation object
    """
    
    integrator_min = openmm.LangevinIntegrator(temperature, friction, time_step)

    platform = Platform.getPlatformByName(platform_name.upper())

    # set up an openmm simulation
    simulation = app.Simulation(
        modeller.topology, system, integrator_min, platform=platform
    )
    
    simulation.context.setPositions(modeller.positions)
    simulation.context.setVelocitiesToTemperature(temperature)
    
    return simulation
    
def minimize_and_equilibrate(simulation, equilibration_steps):
    """
    Minimize and equilibrate the OpenMM simulation object.
    
    Parameters
    ----------
    simulation: openmm.app.simulation.Simulation
        OpenMM initial simulation object
    equilibration_steps: int
        steps to equilibrate the system.
    
    Returns
    -------
    simulation: openmm.app.simulation.Simulation
        OpenMM prepared simulation object, with the energy minimization and equilibration of the system
    """
    # minimize 
    simulation.minimizeEnergy()
    
    # equilibrate
    simulation.step(equilibration_steps)
    
    return simulation
    
def run_simulation(simulation, num_steps, reporting_interval, output_traj_pdb, output_traj_dcd, temperature=300*unit.kelvin):
    """
    Function inspired from: https://github.com/tdudgeon/simple-simulate-complex/blob/master/simulateComplex.py
    
    Run the OpenMM simulation and add reporters.
    
    Parameters
    ----------
    simulation: openmm.app.simulation.Simulation
        OpenMM initial simulation object
    num_steps: int
        Number of steps to run the simulation
    reporting_interval: int
        Interval at which the reporters are added
    output_traj_pdb: str
        file path to store the trajectory output in .pdb format
    output_traj_dcd: str
        file path to store the trajectory output in .dcd format
    """
    # Run the simulation.
    # The enforcePeriodicBox arg to the reporters is important.
    # It's a bit counter-intuitive that the value needs to be False, but this is needed to ensure that
    # all parts of the simulation end up in the same periodic box when being output.
    simulation.reporters.append(PDBReporter(output_traj_pdb, reporting_interval, enforcePeriodicBox=False))
    simulation.reporters.append(DCDReporter(output_traj_dcd, reporting_interval, enforcePeriodicBox=False))
    simulation.reporters.append(StateDataReporter(sys.stdout, reporting_interval, step=True, potentialEnergy=True, temperature=True))
    print('Starting simulation with', num_steps, 'steps.')
    t0 = time.time()
    simulation.step(num_steps)
    t1 = time.time()
    print('Simulation complete in', t1 - t0, 'seconds at', temperature)
