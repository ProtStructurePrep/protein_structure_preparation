import pprep.utils_mdtraj as md
import pprep.utils_pdbfixer as pdbfix
import sys
import os 


def get_files(pdb_name):
    """
    Create a folder (named <pdb_name>) that will contain:
    - a folder for each of the chains
        - inside each folder the pdb of the protein chain and the pdb of the ligand
        
    Parameters
    ----------
    pdb_name: string
        it can be either the pdb file path or the pdb id of the protein  
    """
    pdb = md.load_pdb(pdb_name=pdb_name)
    protein_chains = md.select_protein_chains(pdb)
    ligands = md.select_ligands(pdb)[1]
    dist = md.compute_distance_chain_ligand(pdb, protein_chains, ligands)
    ligand_chain = md.associate_ligand_to_chain(dist)
    md.save_chain_ligand_pdb(pdb_name, pdb, ligands,protein_chains, ligand_chain, "outputs")
    pdbfix.apply_pdbfixer("outputs")
    return("Files created")

for i in sys.argv[1:]:
    print(i)
    print(get_files(i))

