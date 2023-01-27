from pprep.utils import get_ligands, load_pdb

PDB_NAME = '1J4R'

def test_get_ligands():
    pdb = load_pdb(PDB_NAME)
    assert len(get_ligands(pdb)) == 3

def test_read_pdb():
    pdb = load_pdb(PDB_NAME)
    assert pdb.n_atoms == 2795


def main():
    test_read_pdb()
    test_get_ligands()
        

if __name__ == '__main__':
    main()
