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

def prepare_protein(pdb_file, output_file, ph=7.0):
    """
    inspired from: https://projects.volkamerlab.org/teachopencadd/talktorials/T019_md_simulation.html
    
    Use pdbfixer to prepare the protein from a PDB file. It removes the heterogens, replace the non standard residues,
    and ignores missing residues.

    Parameters
    ----------
    pdb_file: pathlib.Path or str
        PDB file containing the system to simulate.
    ph: float, optional
        pH value used to determine protonation state of residues

    Returns
    -------
    fixer: pdbfixer.pdbfixer.PDBFixer
        Prepared protein system.
    """
    fixer = pdbfixer.PDBFixer(pdb_file)
    fixer.removeHeterogens()   
    fixer.findMissingResidues()
    fixer.missingResidues = {}
    fixer.findNonstandardResidues()  
    fixer.replaceNonstandardResidues()  
    fixer.findMissingAtoms()  
    fixer.addMissingAtoms()  
    fixer.addMissingHydrogens(ph)  
    PDBFile.writeFile(fixer.topology, fixer.positions, open(output_file, 'w'))
    return fixer

def apply_pdbfixer(output_directory):
    """ 
    Iterate over all the protein chains in the `output` directory to create the corrected pdb file. 
    
    Parameters
    ----------
    output_directory: string
        directory where all the chain files are found (is the output directory created by save_chain_ligand_pdb()
        in utils_mdtraj)
    
    """
    
    for directory1 in os.listdir(output_directory): # for each different protein
        for directory2 in os.listdir(f"{output_directory}/{directory1}"): # for each chain-ligand folder
            for file in os.listdir(f"{output_directory}/{directory1}/{directory2}"): # for each protein chain
                if "chain" in file:
                    prepare_protein(f"{output_directory}/{directory1}/{directory2}/{file}",f"{output_directory}/{directory1}/{directory2}/corrected_{file}")
                    
    print("PDBFixer applied")
