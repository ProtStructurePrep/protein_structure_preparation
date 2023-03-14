import pprep.utils_mdtraj as md
import time

PDB_NAMES = ['1j4r', '2fv9', '2jfz', '3fnu', '3qs1', '3kdb', '2bua', '2buc', '3ckr', '3cw8', '1a5g', '1a61', '1d9i',
            '1tbz','1d6w' , '4omd', '4omc', '4ii8', '3n84', '3n7y', '3ckp', '3ouh', '3nls', '4npt', '4dqh', '4dqe', '3jw2', '2hs2', '2hs1', '2idw']
EXPECTED_LIGANDS = ['001','INN','003','006','006','006','007','008','009','00A','00L','00N','00P', '00Q', '00R',
                   '00S', '00S', '010', '011', '011', '012', '014', '016', '017', '017', '017', '017', '017', '017', '017']

def test_select_ligand():
    start = time.time()
    
    l = []
    for i in range(len(PDB_NAMES)):
        print(PDB_NAMES[i])
        pdb = md.load_pdb(PDB_NAMES[i])
        ligands = md.select_ligands(pdb)[0]
        print(ligands)
        assert EXPECTED_LIGANDS[i] in ligands
        
    end = time.time()
    total_time = end - start
    print(total_time)


def main():
    test_select_ligand()


if __name__ == '__main__':
    main()
