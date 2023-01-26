from pprep import utils

'''
create the 'common_ligands' list
'''

data = open("/home/eva/tfg/data/cc-to-pdb.tdd","r")

len(utils.common_ligands(data))

'''
load a pdb
'''

pdb = utils.load_pdb('1J4R')


'''
get the ligands of the protein (ignoring the ligands found in the 'common_ligands' list)
'''

ligands = utils.get_ligands(pdb)
print(ligands)
