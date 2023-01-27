import pdbfixer
from pdbfixer import PDBFixer
from openmm.app import PDBFile

def chains_id(pdbid):
    l = []
    fixer = pdbfixer.PDBFixer(pdbid + '.pdb')
    for c in fixer.topology.chains():
        l.append(c.id)
    return(l)

def quick_fix(pdbid):
    fixer = PDBFixer(filename=pdbid + '.pdb')
    fixer.findMissingResidues()
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.removeHeterogens(True)
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.0)
    fixer.addSolvent(fixer.topology.getUnitCellDimensions())
    PDBFile.writeFile(fixer.topology, fixer.positions, open(pdbid + 'corrected' + '.pdb', 'w'))
    return fixer
