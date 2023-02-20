from openmm.app import PDBFile, ForceField, PME, HBonds, Simulation, PDBReporter, StateDataReporter, Modeller, DCDReporter
from openmm.unit import nanometer, kelvin, picoseconds, picosecond
from openmm import LangevinMiddleIntegrator
from sys import stdout
import copy
from simtk.openmm import app, Platform, LangevinIntegrator
from openmmforcefields.generators import SystemGenerator
from openff.toolkit.topology import Molecule, Topology
from openmm import unit
from pprep.utils_rdkit import prepare_ligand


def simple_protein_simulation(pdb_file, reportInterval, steps):
    """
    Use OpenMM to do a protein simulation.

    Parameters
    ----------
    pdb_file: string
        pdb file path
    reportInterval: int
        The interval (in time steps) at which to write frames
    steps: int
        Number of time steps
    """
    # read protein and force field files
    pdb = PDBFile(pdb_file)
    forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')

    # create the system
    system = forcefield.createSystem(pdb.topology, nonbondedMethod=PME,
            nonbondedCutoff=1*nanometer, constraints=HBonds)
    integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.004*picoseconds)

    # create the simulation context
    simulation = Simulation(pdb.topology, system, integrator)
    simulation.context.setPositions(pdb.positions)

    # minimize the energy
    simulation.minimizeEnergy()

    # add reporteres
    simulation.reporters.append(PDBReporter('output_simple_simulation.pdb', reportInterval))
    simulation.reporters.append(StateDataReporter(stdout, reportInterval, step=True,
            potentialEnergy=True, temperature=True))

    # run the simulation
    simulation.step(steps)

def create_modeller(chain_file, ligand_file, ligand_name, smiles):
    """ 
    Use OpenMM to create the Modeller object, having the protein chain and the ligand.
    
    Parameters
    ----------
    chain_file: str
        PDB file corrected previously by PDBFixer
    ligand_file: str
        PDB file of the ligand
    ligand_name: str
        Name of the ligand
    Smiles: str
        Smiles code of the ligand

    Returns
    -------
    modeller: openmm.app.modeller.Modeller
        Modeller object with the protein chain and the ligand.
    """
    protein = PDBFile(chain_file)
    rdkit_ligand = prepare_ligand(ligand_file, ligand_name, smiles, depict=False)
    ligand_mol = Molecule(rdkit_ligand)
    
    modeller = Modeller(protein.topology, protein.positions)
    modeller.add(ligand_mol.to_topology().to_openmm(), ligand_mol.conformers[0])
    return modeller

def protein_ligand_simulation(modeller, ligand_mol):
    # check whether we have a GPU platform and if so set the precision to mixed
    speed = 0
    for i in range(Platform.getNumPlatforms()):
        p = Platform.getPlatform(i)
        # print(p.getName(), p.getSpeed())
        if p.getSpeed() > speed:
            platform = p
            speed = p.getSpeed()

    if platform.getName() == 'CUDA' or platform.getName() == 'OpenCL':
        platform.setPropertyDefaultValue('Precision', 'mixed')
        print('Set precision for platform', platform.getName(), 'to mixed')


    # Prepare the system
    print('Preparing system')
    forcefield_kwargs = { 'constraints': app.HBonds, 'rigidWater': True, 'removeCMMotion': False, 'hydrogenMass': 4*unit.amu }
    system_generator = SystemGenerator(
        forcefields=['amber/ff14SB.xml'],
        small_molecule_forcefield='gaff-2.11',
        forcefield_kwargs=forcefield_kwargs)

    system = system_generator.create_system(modeller.topology, molecules=ligand_mol)

    integrator = LangevinIntegrator(temperature, 1 / unit.picosecond, 0.002 * unit.picoseconds)
    # system.addForce(openmm.MonteCarloBarostat(1*unit.atmospheres, temperature, 25))
    print('Uses Periodic box:', system.usesPeriodicBoundaryConditions(),
        ', Default Periodic box:', system.getDefaultPeriodicBoxVectors())

    simulation = Simulation(modeller.topology, system, integrator, platform=platform)
    simulation.context.setPositions(modeller.positions)
    print('Minimising ...')
    simulation.minimizeEnergy()

    # write out the minimised PDB
    with open("minimized.pdb", 'w') as outfile:
        PDBFile.writeFile(modeller.topology, simulation.context.getState(getPositions=True, enforcePeriodicBox=False).getPositions(), file=outfile, keepIds=True)

    # equilibrate
    simulation.context.setVelocitiesToTemperature(temperature)
    print('Equilibrating ...')
    simulation.step(equilibration_steps)

    # Run the simulation.
    # The enforcePeriodicBox arg to the reporters is important.
    # It's a bit counter-intuitive that the value needs to be False, but this is needed to ensure that
    # all parts of the simulation end up in the same periodic box when being output.
    # simulation.reporters.append(PDBReporter(output_traj_pdb, reporting_interval, enforcePeriodicBox=False))
    simulation.reporters.append(DCDReporter(output_traj_dcd, reporting_interval, enforcePeriodicBox=False))
    simulation.reporters.append(StateDataReporter(sys.stdout, reporting_interval * 5, step=True, potentialEnergy=True, temperature=True))
    print('Starting simulation with', num_steps, 'steps ...')
    t0 = time.time()
    simulation.step(num_steps)
    t1 = time.time()
    print('Simulation complete in', t1 - t0, 'seconds at', temperature)