data = open("/home/eva/tfg/data/cc-to-pdb.tdd","r")

def common_ligands(data, n=100):
    excluded_ligands = []
    threshold = n  # check threshold

    for line in data: 
        line = line.split()
        ligand = line[0]
        n_prot = len(line[1:])
        # print(ligand, n_prot)
        if n_prot > threshold:
            excluded_ligands.append(ligand)
        
    return(excluded_ligands)

if __name__ == '__main__':
    excluded_ligands = common_ligands(data)
    for ligand in excluded_ligands:
        print(ligand)
