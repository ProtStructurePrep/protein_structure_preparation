from pprep.io import read_pdb

PDB_FILE='1J4R.pdb'


def test_read_pdb():
    pdb = read_pdb(PDB_FILE)
    assert pdb.n_atoms == 2795


def main():
    test_read_pdb()
        

if __name__ == '__main__':
    main()
