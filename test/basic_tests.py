from pprep.io import read_pdb
from pprep.utils import get_ligands

PDB_FILE='1J4R.pdb'


def test_read_pdb():
    pdb = read_pdb(PDB_FILE)
    assert pdb.n_atoms == 2795


def test_get_ligands():
    pdb = read_pdb(PDB_FILE)
    assert len(get_ligands(pdb)) == 3


def main():
    test_read_pdb()
    test_get_ligands()
        

if __name__ == '__main__':
    main()
