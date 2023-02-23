import pdbfixer
from pdbfixer import PDBFixer
from openmm.app import PDBFile
import os


def chains_id(pdb_file):
    """
    Use pdbfixer to get a list of the ids of the chains in a protein.

    Parameters
    ----------
    pdb_file: string
        pdb file path

    Returns
    -------
    ids: list
        contains all the chains ids.
    """
    
    ids = []
    fixer = pdbfixer.PDBFixer(pdb_file)
    
    for c in fixer.topology.chains():
        ids.append(c.id)
        
    return ids

def prepare_protein(
    pdb_file, output_file, ignore_missing_residues=True, ignore_terminal_missing_residues=True, ph=7.0):
    """
    Use pdbfixer to prepare the protein from a PDB file. Hetero atoms such as ligands are
    removed and non-standard residues replaced. Missing atoms to existing residues are added.
    Missing residues are ignored by default, but can be included.

    Parameters
    ----------
    pdb_file: pathlib.Path or str
        PDB file containing the system to simulate.
    ignore_missing_residues: bool, optional
        If missing residues should be ignored or built.
    ignore_terminal_missing_residues: bool, optional
        If missing residues at the beginning and the end of a chain should be ignored or built.
    ph: float, optional
        pH value used to determine protonation state of residues

    Returns
    -------
    fixer: pdbfixer.pdbfixer.PDBFixer
        Prepared protein system.
    """
    fixer = pdbfixer.PDBFixer(str(pdb_file))
    fixer.removeHeterogens()  # co-crystallized ligands are unknown to PDBFixer
    fixer.findMissingResidues()  # identify missing residues, needed for identification of missing atoms

    # if missing terminal residues shall be ignored, remove them from the dictionary
    if ignore_terminal_missing_residues:
        chains = list(fixer.topology.chains())
        keys = fixer.missingResidues.keys()
        for key in list(keys):
            chain = chains[key[0]]
            if key[1] == 0 or key[1] == len(list(chain.residues())):
                del fixer.missingResidues[key]

    # if all missing residues shall be ignored ignored, clear the dictionary
    if ignore_missing_residues:
        fixer.missingResidues = {}

    fixer.findNonstandardResidues()  # find non-standard residue
    fixer.replaceNonstandardResidues()  # replace non-standard residues with standard one
    fixer.findMissingAtoms()  # find missing heavy atoms
    fixer.addMissingAtoms()  # add missing atoms and residues
    fixer.addMissingHydrogens(ph)  # add missing hydrogens
    PDBFile.writeFile(fixer.topology, fixer.positions, open(output_file, 'w'))
    return fixer

def apply_pdbfixer():
    """ Iterate over all the protein chains in the `output` directory to create the corrected pdb file. """
    
    for directory1 in os.listdir("outputs/"): # for each different protein
        for directory2 in os.listdir(f"outputs/{directory1}"): # for each chain-ligand folder
            for file in os.listdir(f"outputs/{directory1}/{directory2}"): # for each protein chain
                if "chain" in file:
                    prepare_protein(f"outputs/{directory1}/{directory2}/{file}",f"outputs/{directory1}/{directory2}/corrected_{file}")
                    
    print("PDBFixer applied")
