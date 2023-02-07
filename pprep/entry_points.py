import pprep.utils_mdtraj as md
import nglview as nv
import numpy as np 
import pprep.utils_pdbfixer as pdbfix
import sys


def get_files(pdb_name):
    pdb = md.load_pdb(pdb_name=pdb_name)
    protein_chains = md.select_protein_chains(pdb)
    all_ligands = md.select_ligands(pdb)
    common_ligands = md.select_common_ligands(pdb, all_ligands)
    ligands = md.select_definitive_ligands(all_ligands,common_ligands) # exclude common ligands
    dist = md.compute_distance_chain_ligand(pdb, protein_chains, ligands)
    ligand_chain = md.associate_ligand_to_chain(dist)
    md.save_chain_ligand_pdb(pdb, ligands,protein_chains, ligand_chain)
    return("Files created")


def cli():
    pdb = sys.argv[1]
    get_files(pdb)

