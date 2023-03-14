import pprep.utils_mdtraj as md

PDB_NAMES = ['1j4r', '2fv9', '2jfz', '3fnu', '3qs1', '3kdb', '2bua', '2buc', '3ckr', '3cw8', '1a5g', '1a61', '1d9i', 
            '1tbz','1d6w' , '4omd', '4omc', '4ii8', '3n84', '3n7y']
EXPECTED_LIGANDS = ['001','INN','003','006','006','006','007','008','009','00A','00L','00N','00P', '00Q', '00R', 
                   '00S', '00S', '010', '011', '011']

def test_select_ligand():
    l = []
    for i in range(len(PDB_NAMES)):
        print(PDB_NAMES[i])
        pdb = md.load_pdb(PDB_NAMES[i])
        ligands = md.select_ligands(pdb)[0]
        print(ligands)
        assert EXPECTED_LIGANDS[i] in ligands

def main():
    test_select_ligand()


if __name__ == '__main__':
    main()
