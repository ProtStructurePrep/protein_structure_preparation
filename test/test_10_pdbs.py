import pprep.utils_mdtraj as md

PDB_NAMES = ['1j4r', '2fv9', '2jfz', '3fnu', '3qs1', '3kdb', '2bua', '2buc', '3ckr', '3cw8']
EXPECTED_LIGANDS = ['001','INN','003','006','006','006','007','008','009','00A']

def test_select_ligand():
    l = []
    for i in range(len(PDB_NAMES)):
        print(PDB_NAMES[i])
        pdb = md.load_pdb(PDB_NAMES[i])
        ligand = md.select_ligands(pdb)[0][0]
        print(ligand)
        assert ligand == EXPECTED_LIGANDS[i]

def main(): 
    test_select_ligand()


if __name__ == '__main__':
    main()

