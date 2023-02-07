import pprep.utils_mdtraj as md

PDB_NAME = '1J4R'
pdb = md.load_pdb(PDB_NAME)

def test_get_ligands():
    assert len(md.get_ligands(pdb)) == 3

def test_read_pdb():
    assert pdb.n_atoms == 2795

def test_associate_ligand():
    protein_chains = md.select_protein_chains(pdb)
    all_ligands = md.select_ligands(pdb)
    common_ligands = md.select_common_ligands(pdb, all_ligands)
    ligands = md.select_definitive_ligands(all_ligands,common_ligands) # exclude common ligands
    dist = md.compute_distance_chain_ligand(pdb, protein_chains, ligands)
    ligand_chain = md.associate_ligand_to_chain(dist)
    assert ligand_chain == {0: 0, 1: 1, 2: 2}


def main():
    test_read_pdb()
    test_get_ligands()
        

if __name__ == '__main__':
    main()
