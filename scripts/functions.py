import mdtraj as md

pdb = md.load_pdb('http://www.rcsb.org/pdb/files/1J4R.pdb')

from pprep import utils

'''
create the 'common_ligands' list
'''

data = open("/home/eva/tfg/data/cc-to-pdb.tdd","r")

len(utils.common_ligands(data))


'''
get the ligands of the protein (ignoring the ligands found in the 'common_ligands' list)
'''

utils.get_ligands(pdb)
