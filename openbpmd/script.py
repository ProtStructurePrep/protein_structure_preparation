import pprep.utils_rdkit as urk
import pprep.MDutils as umd
import sys

# OpenMM
import openmm
from openmm.app import Modeller, PDBFile, Simulation,  PDBReporter, StateDataReporter, DCDReporter
from openmm import Platform, LangevinIntegrator, app, HarmonicBondForce, NonbondedForce
from simtk import unit
from openmmforcefields.generators import SystemGenerator
from openmm.app.metadynamics import BiasVariable, Metadynamics
from openff.units.openmm import to_openmm

from openmm.openmm import XmlSerializer

# Others
import os
import MDAnalysis as mda
import mdtraj as md
import numpy as np
import pandas as pd

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
    
"""RDKit prepared ligand --> OpenMM ligand"""
pdb_file = 'ligand.pdb'
resname = '03P'
smiles = 'O=C(CC(O)(C)C)NCCn1ccc2c1c(ncn2)Nc1ccc(c(c1)Cl)Oc1cccc(c1)C(F)(F)F'
rdkit_ligand = urk.prepare_ligand(pdb_file, resname, smiles)

ligand = umd.load_prepared_ligand(rdkit_ligand)

"""OpenMM receptor"""
receptor = umd.load_prepared_receptor('prepared_receptor.pdb')

modeller = Modeller(receptor.topology, receptor.positions)
modeller.add(ligand.to_topology().to_openmm(), ligand.conformers[0])
#modeller.add(ligand.to_topology().to_openmm(), to_openmm(ligand.conformers[0]))


forcefield_kwargs = {'constraints': app.HBonds, 'rigidWater': True, 'removeCMMotion': False, 'hydrogenMass': 4*unit.amu }
system_generator = SystemGenerator(
    forcefields=['amber/ff14SB.xml', 'amber/tip3p_standard.xml'],
    small_molecule_forcefield='gaff-2.11',
    molecules=[ligand],
    forcefield_kwargs=forcefield_kwargs)
modeller.addSolvent(system_generator.forcefield, model='tip3p', padding=10.0*unit.angstroms)

with open('solvated_complex.pdb', 'w') as outfile:
    PDBFile.writeFile(modeller.topology, modeller.positions, outfile)

system = system_generator.create_system(modeller.topology, molecules=ligand)   

with open('solvated_complex.xml', 'w') as output:
    output.write(XmlSerializer.serialize(system))
