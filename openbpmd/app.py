import pprep.utils_rdkit as urk
import pprep.MDutils as umd
import sys
from pathlib import Path
import json
from copy import deepcopy

# OpenMM
import openmm
from openmm.app import Modeller, PDBFile, Simulation,  PDBReporter, StateDataReporter, DCDReporter
from openmm import Platform, LangevinIntegrator, app, HarmonicBondForce, NonbondedForce
from simtk import unit
from openmmforcefields.generators import SystemGenerator
from openmm.app.metadynamics import BiasVariable, Metadynamics
from openff.units.openmm import to_openmm
from openff.toolkit.topology import Molecule as OFFMolecule

from openmm.openmm import XmlSerializer

# Others
import os
import MDAnalysis as mda
from MDAnalysis.analysis import rms, contacts
import mdtraj as md
import numpy as np
import pandas as pd
import parmed as pmd
import argparse

import glob
import os

import logging
logger = logging.getLogger(__name__)
logging.basicConfig(level=os.environ.get("LOGLEVEL", "INFO"))


OUTPUT_DIR = Path('Outputs')
NREPS = 10
DEFAULT_LIG_RESNAME = 'LIG'

RECEPTOR_PDB_FILE = 'prepared_receptor.pdb'
LIGAND_PDB_FILE = 'ligand.pdb'
LIG_RESNAME = '03P'
SMILES = 'O=C(CC(O)(C)C)NCCn1ccc2c1c(ncn2)Nc1ccc(c(c1)Cl)Oc1cccc(c1)C(F)(F)F'


def _system_generator(ligand, restart=False):
    if restart or not hasattr(_system_generator, 'generator'):
        forcefield_kwargs = {
            'constraints': app.HBonds, 'rigidWater': True,
            'removeCMMotion': False, 'hydrogenMass': 4*unit.amu }
        generator = SystemGenerator(
            forcefields=['amber/ff14SB.xml', 'amber/tip3p_standard.xml'],
            small_molecule_forcefield='gaff-2.11',
            molecules=[ligand],
            forcefield_kwargs=forcefield_kwargs)
        _system_generator.generator = generator
    return _system_generator.generator

def _get_platform():
    """Get platform."""
    speed = 0
    for i in range(Platform.getNumPlatforms()):
        p = Platform.getPlatform(i)
        # print(p.getName(), p.getSpeed())
        if p.getSpeed() > speed:
            platform = p
            speed = p.getSpeed()
            logger.info(f'Platform: {platform.getName()}, speed: {speed}')
    if platform.getName() == 'CUDA' or platform.getName() == 'OpenCL':
        platform.setPropertyDefaultValue('Precision', 'mixed')
    logger.info(f'Set precision for platform {platform.getName()} to mixed')
    return platform


def _setup_system(receptor_pdb_file, ligand_pdb_file, smiles,
                  lig_resname='LIG'):
    """Return topology, system."""

    "RDKit prepared ligand --> OpenMM ligand"
    rdkit_ligand = urk.prepare_ligand(ligand_pdb_file, lig_resname, smiles)
    
    ligand = umd.load_prepared_ligand(rdkit_ligand)
    ligand.name = DEFAULT_LIG_RESNAME

    "OpenMM receptor"
    receptor = umd.load_prepared_receptor(receptor_pdb_file)

    modeller = Modeller(receptor.topology, receptor.positions)
    modeller.add(ligand.to_topology().to_openmm(), ligand.conformers[0])
    #modeller.add(ligand.to_topology().to_openmm(), to_openmm(ligand.conformers[0]))
    
    system_generator = _system_generator(ligand)
    modeller.addSolvent(system_generator.forcefield, model='tip3p', padding=10.0*unit.angstroms)


    with open('off_ligand.json', 'w') as fp:
        print(ligand.to_json(), file=fp)

    with open('solvated_complex.pdb', 'w') as outfile:
        PDBFile.writeFile(modeller.topology, modeller.positions, outfile)

    topology = modeller.topology
    system = system_generator.create_system(topology, molecules=ligand)

    with open('solvated_complex.xml', 'w') as output:
        output.write(XmlSerializer.serialize(system))
    return modeller, topology, system, receptor, ligand


"""Minimize"""
def minimize(parm, input_positions, system, platform, out_dir, min_file_name):
    """An energy minimization function down with an energy tolerance
    of 10 kJ/mol.

    Parameters
    ----------
    parm : Parmed or OpenMM parameter file object
        Used to create the OpenMM System object.
    input_positions : OpenMM Quantity
        3D coordinates of the equilibrated system.
    system : 
        
    platform : 
        
    out_dir : str
        Directory to write the outputs.
    min_file_name : str
        Name of the minimized PDB file to write.
    """

    # Set up the simulation parameters
    # Langevin integrator at 300 K w/ 1 ps^-1 friction coefficient
    # and a 2-fs timestep
    # NOTE - no dynamics performed, but required for setting up
    # the OpenMM system.
    integrator = LangevinIntegrator(300*unit.kelvin, 1/unit.picosecond,
                                    0.002*unit.picoseconds)
    simulation = Simulation(parm.topology, system, integrator, platform)
    simulation.context.setPositions(input_positions)

    # Minimize the system - no predefined number of steps
    simulation.minimizeEnergy()

    # Write out the minimized system to use w/ MDAnalysis
    positions = simulation.context.getState(getPositions=True).getPositions()
    out_file = os.path.join(out_dir,min_file_name)
    PDBFile.writeFile(simulation.topology, positions,
                      open(out_file, 'w'))

    return None


"""Equilibrate"""
def equilibrate(min_pdb, parm, system, platform, out_dir, eq_file_name):
    """A function that does a 500 ps NVT equilibration with position
    restraints, with a 5 kcal/mol/A**2 harmonic constant on solute heavy
    atoms, using a 2 fs timestep.

    Parameters
    ----------
    min_pdb : str
        Name of the minimized PDB file.
    parm : Parmed or OpenMM parameter file object
        Used to create the OpenMM System object.
    system : 
    
    platform :
    
    out_dir : str
        Directory to write the outputs to.
    eq_file_name : str
        Name of the equilibrated PDB file to write.
    """
    NSTEPS=250 # 500 fs
    # Get the solute heavy atom indices to use
    # for defining position restraints during equilibration
    universe = mda.Universe(min_pdb,
                            format='XPDB', in_memory=True)
    solute_heavy_atom_idx = universe.select_atoms('not resname WAT and\
                                                   not resname SOL and\
                                                   not resname HOH and\
                                                   not resname CL and \
                                                   not resname NA and \
                                                   not name H*').indices
    # Necessary conversion to int from numpy.int64,
    # b/c it breaks OpenMM C++ function
    solute_heavy_atom_idx = [int(idx) for idx in solute_heavy_atom_idx]

    # Add the restraints.
    # We add a dummy atoms with no mass, which are therefore unaffected by
    # any kind of scaling done by barostat (if used). And the atoms are
    # harmonically restrained to the dummy atom. We have to redefine the
    # system, b/c we're adding new particles and this would clash with
    # modeller.topology.
    
    # Add the harmonic restraints on the positions
    # of specified atoms
    system = deepcopy(system)  # copy system before modifications
    restraint = HarmonicBondForce()
    restraint.setUsesPeriodicBoundaryConditions(True)
    system.addForce(restraint)
    nonbonded = [force for force in system.getForces()
                 if isinstance(force, NonbondedForce)][0]
    dummyIndex = []
    input_positions = PDBFile(min_pdb).getPositions()
    positions = input_positions
    # Go through the indices of all atoms that will be restrained
    for i in solute_heavy_atom_idx:
        j = system.addParticle(0)
        # ... and add a dummy/ghost atom next to it
        nonbonded.addParticle(0, 1, 0)
        # ... that won't interact with the restrained atom 
        nonbonded.addException(i, j, 0, 1, 0)
        # ... but will be have a harmonic restraint ('bond')
        # between the two atoms
        restraint.addBond(i, j, 0*unit.nanometers,
                          5*unit.kilocalories_per_mole/unit.angstrom**2)
        dummyIndex.append(j)
        input_positions.append(positions[i])

    integrator = LangevinIntegrator(300*unit.kelvin, 1/unit.picosecond,
                                    0.002*unit.picoseconds)
   
    sim = Simulation(parm.topology, system, integrator,
                     platform)
    sim.context.setPositions(input_positions)
    integrator.step(NSTEPS)  # run 500 ps of equilibration
    all_positions = sim.context.getState(
        getPositions=True, enforcePeriodicBox=True).getPositions()
    # we don't want to write the dummy atoms, so we only
    # write the positions of atoms up to the first dummy atom index
    relevant_positions = all_positions[:dummyIndex[0]]
    out_file = os.path.join(out_dir,eq_file_name)
    PDBFile.writeFile(sim.topology, relevant_positions,
                      open(out_file, 'w'))

    return None


"""Produce"""
def produce(out_dir, idx, eq_pdb, parm, mdtraj_top, set_hill_height, platform,
            ligand, system, lig_resname='LIG'):
    """An OpenBPMD production simulation function. Ligand RMSD is biased with
    metadynamics. The integrator uses a 4 fs time step and
    runs for 10 ns, writing a frame every 100 ps.

    Writes a 'trj.dcd', 'COLVAR.npy', 'bias_*.npy' and 'sim_log.csv' files
    during the metadynamics simulation in the '{out_dir}/rep_{idx}' directory.
    After the simulation is done, it analyses the trajectories and writes a
    'bpm_results.csv' file with time-resolved PoseScore and ContactScore.

    Parameters
    ----------
    out_dir : str
        Directory where your equilibration PDBs and 'rep_*' dirs are at.
    idx : int
        Current replica index.
    lig_resname : str
        Residue name of the ligand.
    eq_pdb : str
        Name of the PDB for equilibrated system.
    parm : Parmed or OpenMM parameter file object
        Used to create the OpenMM System object.
    parm_file : str
        The name of the parameter or topology file of the system.
    coords_file : str
        The name of the coordinate file of the system.
    set_hill_height : float
        Metadynamic hill height, in kcal/mol.
    """
    # the system needs to be created again, as the system has been modified when doing the equilibration
    # I figured this out opening a GitHub issue: https://github.com/openmm/openmm/issues/4091#issuecomment-1570916065

    #system_generator = _system_generator(ligand)
    #system = system_generator.create_system(modeller.topology, molecules=ligand)
    
    # First, assign the replica directory to which we'll write the files
    write_dir = os.path.join(out_dir,f'rep_{idx}')
    # Get the anchor atoms by ...
    universe = mda.Universe(eq_pdb,
                            format='XPDB', in_memory=True)
    # ... finding the protein's COM ...
    prot_com = universe.select_atoms('protein').center_of_mass()
    x, y, z = prot_com[0], prot_com[1], prot_com[2]
    # ... and taking the heavy backbone atoms within 5A of the COM
    sel_str = f'point {x} {y} {z} 5 and backbone and not name H*'
    anchor_atoms = universe.select_atoms(sel_str)
    # ... or 10 angstrom
    if len(anchor_atoms) == 0:
        sel_str = f'point {x} {y} {z} 10 and backbone and not name H*'
        anchor_atoms = universe.select_atoms(sel_str)

    anchor_atom_idx = anchor_atoms.indices.tolist()
    # Get indices of ligand heavy atoms
    lig = universe.select_atoms('not protein and not resname HOH')

    lig_ha_idx = lig.indices.tolist()
    # get the atom positions for the system from the equilibrated
    # system
    input_positions = PDBFile(eq_pdb).getPositions()

    # Add an 'empty' flat-bottom restraint to fix the issue with PBC.
    # Without one, RMSDForce object fails to account for PBC.
    k = 0*unit.kilojoules_per_mole  # NOTE - 0 kJ/mol constant
    upper_wall = 10.00*unit.nanometer
    fb_eq = '(k/2)*max(distance(g1,g2) - upper_wall, 0)^2'
    upper_wall_rest = openmm.CustomCentroidBondForce(2, fb_eq)
    upper_wall_rest.addGroup(lig_ha_idx)
    upper_wall_rest.addGroup(anchor_atom_idx)
    upper_wall_rest.addBond([0, 1])
    upper_wall_rest.addGlobalParameter('k', k)
    upper_wall_rest.addGlobalParameter('upper_wall', upper_wall)
    upper_wall_rest.setUsesPeriodicBoundaryConditions(True)
    system.addForce(upper_wall_rest)

    alignment_indices = lig_ha_idx + anchor_atom_idx

    rmsd = openmm.RMSDForce(input_positions, alignment_indices)
    # Set up the typical metadynamics parameters
    grid_min, grid_max = 0.0, 1.0  # nm
    hill_height = set_hill_height*unit.kilocalories_per_mole
    hill_width = 0.002  # nm, also known as sigma

    grid_width = hill_width / 5
    # 'grid' here refers to the number of grid points
    grid = int(abs(grid_min - grid_max) / grid_width)

    rmsd_cv = BiasVariable(rmsd, grid_min, grid_max, hill_width,
                           False, gridWidth=grid)

    # define the metadynamics object
    # deposit bias every 1 ps, BF = 4, write bias every ns
    meta = Metadynamics(system, [rmsd_cv], 300.0*unit.kelvin, 4.0, hill_height,
                        250, biasDir=write_dir,
                        saveFrequency=250000)

    # Set up and run metadynamics
    integrator = LangevinIntegrator(300*unit.kelvin, 1.0/unit.picosecond,
                                    0.004*unit.picoseconds)

    simulation = Simulation(parm.topology, system, integrator, platform)
    simulation.context.setPositions(input_positions)

    trj_name = os.path.join(write_dir,'trj.dcd')

    sim_time = 10  # ns
    steps = 250000 * sim_time

    simulation.reporters.append(DCDReporter(trj_name, 25000))  # every 100 ps
    simulation.reporters.append(StateDataReporter(
                                os.path.join(write_dir,'sim_log.csv'), 250000,
                                step=True, temperature=True, progress=True,
                                remainingTime=True, speed=True,
                                totalSteps=steps, separator=','))  # every 1 ns

    colvar_array = np.array([meta.getCollectiveVariables(simulation)])
    for i in range(0, int(steps), 500):
        if i % 25000 == 0:
            # log the stored COLVAR every 100ps
            np.save(os.path.join(write_dir,'COLVAR.npy'), colvar_array)
        meta.step(simulation, 500)
        current_cvs = meta.getCollectiveVariables(simulation)
        # record the CVs every 2 ps
        colvar_array = np.append(colvar_array, [current_cvs], axis=0)
    np.save(os.path.join(write_dir,'COLVAR.npy'), colvar_array)

    # center everything using MDTraj, to fix any PBC imaging issues
    
    mdu = md.load(trj_name, top=mdtraj_top)
    mdu.image_molecules()
    mdu.save(trj_name)

    return None
    
"""Pose score"""
def get_pose_score(structure_file, trajectory_file):
    """A function the gets the PoseScore (ligand RMSD) from an OpenBPMD
    trajectory.

    Parameters
    ----------
    'structure_file : str
        The name of the centred equilibrated system
        PDB file that was used to start the OpenBPMD simulation.
    trajectory_file : str
        The name of the OpenBPMD trajectory file.
    lig_resname : str
        Residue name of the ligand that was biased.

    Returns
    -------
    pose_scores : np.array 
        PoseScore for every frame of the trajectory.
    """
    # Load a MDA universe with the trajectory
    u = mda.Universe(structure_file, trajectory_file)
    # Align each frame using the backbone as reference
    # Calculate the RMSD of ligand heavy atoms
    r = rms.RMSD(u, select='backbone',
                 groupselections=['not protein and not resname HOH'],
                 ref_frame=0).run()
    # Get the PoseScores as np.array
    pose_scores = r.rmsd[1:, -1]

    return pose_scores

"""Contact score"""
def get_contact_score(structure_file, trajectory_file, lig_resname='LIG'):
    """A function the gets the ContactScore from an OpenBPMD trajectory.

    Parameters
    ----------
    structure_file : str
        The name of the centred equilibrated system PDB file that 
        was used to start the OpenBPMD simulation.
    trajectory_file : str
        The name of the OpenBPMD trajectory file.
    lig_resname : str
        Residue name of the ligand that was biased.

    Returns
    -------
    contact_scores : np.array 
        ContactScore for every frame of the trajectory.
    """
    u = mda.Universe(structure_file, trajectory_file)

    sel_donor = f"resname {lig_resname} and not name *H*"
    sel_acceptor = f"protein and not name H* and \
                     around 5 resname {lig_resname}"

    # reference groups (first frame of the trajectory, but you could also use
    # a separate PDB, eg crystal structure)
    a_donors = u.select_atoms(sel_donor)
    a_acceptors = u.select_atoms(sel_acceptor)

    cont_analysis = contacts.Contacts(u, select=(sel_donor, sel_acceptor),
                                      refgroup=(a_donors, a_acceptors),
                                      radius=3.5)

    cont_analysis.run()
    # print number of average contacts in the first ns
    # NOTE - hard coded number of frames (100 per traj)
    frame_idx_first_ns = int(len(cont_analysis.timeseries)/10)
    first_ns_mean = np.mean(cont_analysis.timeseries[1:frame_idx_first_ns, 1])
    if first_ns_mean == 0:
        normed_contacts = cont_analysis.timeseries[1:, 1]
    else:
        normed_contacts = cont_analysis.timeseries[1:, 1]/first_ns_mean
    contact_scores = np.where(normed_contacts > 1, 1, normed_contacts)

    return contact_scores

"""Collect results"""
def collect_results(in_dir, out_dir):
    """A function that collects the time-resolved BPM results,
    takes the scores from last 2 ns of the simulation, averages them
    and writes that average as the final score for a given pose.

    Writes a 'results.csv' file in 'out_dir' directory.
    
    Parameters
    ----------
    in_dir : str
        Directory with 'rep_*' directories.
    out_dir : str
        Directory where the 'results.csv' file will be written
    """
    compList = []
    contactList = []
    poseList = []
    # find how many repeats have been run
    glob_str = os.path.join(in_dir,'rep_*')
    nreps = len(glob.glob(glob_str))
    for idx in range(0, nreps):
        f = os.path.join(in_dir,f'rep_{idx}','bpm_results.csv')
        df = pd.read_csv(f)
        # Since we only want last 2 ns, get the index of
        # the last 20% of the data points
        last_2ns_idx = round(len(df['CompScore'].values)/5)  # round up
        compList.append(df['CompScore'].values[-last_2ns_idx:])
        contactList.append(df['ContactScore'].values[-last_2ns_idx:])
        poseList.append(df['PoseScore'].values[-last_2ns_idx:])

    # Get the means of the last 2 ns
    meanCompScore = np.mean(compList)
    meanPoseScore = np.mean(poseList)
    meanContact = np.mean(contactList)
    # Get the standard deviation of the final 2 ns
    meanCompScore_std = np.std(compList)
    meanPoseScore_std = np.std(poseList)
    meanContact_std = np.std(contactList)
    # Format it the Pandas way
    d = {'CompScore': [meanCompScore], 'CompScoreSD': [meanCompScore_std],
         'PoseScore': [meanPoseScore], 'PoseScoreSD': [meanPoseScore_std],
         'ContactScore': [meanContact], 'ContactScoreSD': [meanContact_std]}

    results_df = pd.DataFrame(data=d)
    results_df = results_df.round(3)
    results_df.to_csv(os.path.join(out_dir,'results.csv'), index=False)


def run_replica(idx, cent_eq_pdb, modeller, platform, ligand, system,
                lig_resname='LIG'):
    """Run a replica simulation."""
    mdtraj_topology = md.Topology.from_openmm(modeller.topology)
    rep_dir = os.path.join('Outputs', f'rep_{idx}')
    if not os.path.isdir(rep_dir):
        os.mkdir(rep_dir)

    if os.path.isfile(os.path.join(rep_dir,'bpm_results.csv')):
        return
        
    produce('Outputs/', idx, 'Outputs/equilibrated.pdb', modeller, mdtraj_topology, 0.3, platform, ligand, system)
                
    trj_name = os.path.join(rep_dir,'trj.dcd')
                
    PoseScoreArr = get_pose_score(cent_eq_pdb, trj_name)

    ContactScoreArr = get_contact_score(cent_eq_pdb, trj_name, lig_resname)
    
    # Calculate the CompScore at every frame
    CompScoreArr = np.zeros(99)
    for index in range(ContactScoreArr.shape[0]):
        ContactScore, PoseScore = ContactScoreArr[index], PoseScoreArr[index]
        CompScore = PoseScore - 5 * ContactScore
        CompScoreArr[index] = CompScore

    Scores = np.stack((CompScoreArr, PoseScoreArr, ContactScoreArr), axis=-1)

    # Save a DataFrame to CSV
    df = pd.DataFrame(Scores, columns=['CompScore', 'PoseScore',
                                           'ContactScore'])
    df.to_csv(os.path.join(rep_dir,'bpm_results.csv'), index=False)


platform = _get_platform()

if os.path.isfile('solvated_complex.xml'):
    logger.info('xml file exists')
    with open('solvated_complex.xml') as input:
        system = XmlSerializer.deserialize(input.read())
    modeller = PDBFile('solvated_complex.pdb')
    topology = modeller.topology
    with open('off_ligand.json', 'r') as fp:
        ligand = OFFMolecule.from_json(fp.read())    
else:
    print('xml file does not exist')
    modeller, topology, system, receptor, ligand = _setup_system(
        RECEPTOR_PDB_FILE, LIGAND_PDB_FILE, SMILES,
        lig_resname=LIG_RESNAME)



# create output directory
OUTPUT_DIR.mkdir(exist_ok=True)

minimized_pdb = OUTPUT_DIR / 'minimized.pdb'
if not minimized_pdb.is_file():
    logger.info("Minimizing...")
    minimize(modeller, modeller.positions, system, platform,
             OUTPUT_DIR, minimized_pdb.name)
else:
    logger.info("Minimized structure exists")

equilibrated_pdb = OUTPUT_DIR / 'equilibrated.pdb'
if not equilibrated_pdb.is_file():
    logger.info("Equilibrating...")
    equilibrate(str(minimized_pdb), modeller, system, platform, OUTPUT_DIR,
                equilibrated_pdb.name)
else:
    logger.info("Equilibrated structure exists")

cent_eq_pdb = OUTPUT_DIR / 'cent_equilibrated.pdb'
if not cent_eq_pdb.is_file():
    logger.info("Centering")
    mdtraj_topology = md.Topology.from_openmm(modeller.topology)
    mdu = md.load(equilibrated_pdb, top=mdtraj_topology)
    mdu.image_molecules()
    mdu.save_pdb(cent_eq_pdb)
else:
    logger.info("Centered equilibrated structure exists")

"""OpenBPMD loop"""
logger.info('Starting OpenBPMD loop')
# Run NREPS number of production simulations (nrep=0 has already been finished above)
for idx in range(0, NREPS):
    run_replica(idx, cent_eq_pdb, modeller, platform, ligand, system)
    
logger.info('Collecting results')
collect_results('Outputs', 'Outputs')
