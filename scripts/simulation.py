from openmm.app import *
from openmm import *
from openmm.unit import *
from sys import stdout

# reading protein and force field files
pdb = PDBFile('outputs/out_3poz/chain_lig_0/corrected_chain_0.pdb')
forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')


# creating the system
system = forcefield.createSystem(pdb.topology, nonbondedMethod=PME,
        nonbondedCutoff=1*nanometer, constraints=HBonds)
integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.004*picoseconds)

# create the simulation context
simulation = Simulation(pdb.topology, system, integrator)
simulation.context.setPositions(pdb.positions)

# minimize the energy
simulation.minimizeEnergy()

# adding reporteres
simulation.reporters.append(PDBReporter('output.pdb', 1000))
simulation.reporters.append(StateDataReporter(stdout, 1000, step=True,
        potentialEnergy=True, temperature=True))

# running simulation
simulation.step(10000)
