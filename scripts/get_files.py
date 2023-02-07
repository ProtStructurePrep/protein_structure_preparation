import pprep.utils_mdtraj as md
import sys

def get_files(pdb_name):
    pdb = md.load_pdb(pdb_name=pdb_name)
    protein_chains = md.select_protein_chains(pdb)
    ligands = md.select_ligands(pdb)
    dist = md.compute_distance_chain_ligand(pdb, protein_chains, ligands)
    ligand_chain = md.associate_ligand_to_chain(dist)
    md.save_chain_ligand_pdb(pdb_name, pdb, ligands,protein_chains, ligand_chain)
    return("Files created")

for i in sys.argv[1:]:
    print(i)
    print(get_files(i))

