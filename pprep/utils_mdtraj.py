import mdtraj as md
import numpy as np
import itertools 
import os
from pathlib import Path

COMMON_LIGANDS = '''
1PE
3DR
5CM
5MC
5MU
7MG
8OG
ABA
ACE
ACO
ACP
ACT
ACY
ADN
ADP
AKG
ALF
ALY
AMP
ANP
ATP
AZI
BA
BCL
BCT
BEF
BEN
BEZ
BGC
BMA
BME
BOG
BR
BRU
BTB
C8E
CA
CAC
CAS
CD
CIT
CL
CME
CMO
CO
CO3
COA
CSD
CSO
CSW
CU
CU1
CYN
DAL
DGT
DIO
DMS
DOC
DTP
DTT
EDO
EOH
EPE
F3S
FAD
FE
FE2
FES
FLC
FME
FMN
FMT
FUC
FUL
GAL
GDP
GLA
GLC
GNP
GOL
GSH
GTP
H4B
HEC
HED
HEM
HG
HYP
IMD
IOD
IPA
K
KCX
LDA
LLP
LMT
M3L
MAL
MAN
MES
MG
MLI
MLY
MN
MPD
MRD
MYR
NA
NAD
NAG
NAI
NAP
NCO
NDG
NDP
NH2
NH4
NI
NO
NO3
OCS
OXY
P6G
PCA
PE4
PEG
PG4
PGE
PLM
PLP
PO4
POP
PSU
PTR
PYR
RET
SAH
SAM
SCN
SEP
SF4
SIA
SIN
SO4
SPM
SR
SUC
THP
TLA
TPO
TPP
TRS
TYS
UDP
UMP
UNX
XE
XYP
ZN
'''.split()

def common_residues(data, n=100):
    """
    
    Get the pdb residues that are common.
    
    Parameters
    ----------
    data: .tdd file
		it contains the residues id in the first column and the proteins 
		in wich it appears in the second column
	n: integer
		threshold to consider a residue as a common residue. The default value is
		n = 100, meaning that if a residue appears in more than 100 proteins,
		this residue will be considered a common residue.
		
    Returns
    ----------
    List of the common residues.

    """
    excluded_ligands = []
    threshold = n # check threshold

    for line in data:
        line = line.split()
        ligand = line[0]
        n_prot = len(line[1:])
        # print(ligand, n_prot)
        if n_prot > threshold:
            excluded_ligands.append(ligand)

    return(excluded_ligands)

def load_pdb(pdb_name):
    """
    
    Parse a pdb file/fetch a structure from the PDB.
    
    Parameters
    ----------
    pdb_name: string
        it can be either: the pdb file path or the pdb id of the protein
    
    Returns
    ----------
    Mdtraj object.

    Example: load_pdb('/home/username/1.pdb')
    
    """
    if Path(pdb_name).is_file():
        return md.load_pdb(pdb_name)
    else:
        # pdb_name is a PDB ID
        link = f'http://www.rcsb.org/pdb/files/{pdb_name}.pdb'
        return md.load_pdb(link)
        
def select_chain(pdb, n):
    """

    Select atomic ids for protein chain with index `n`
    
    Parameters
    ----------
    pdb: mdtraj object
        it is the object that load_pdb() returns
    n: integer
        it indicates the index of the chain we want to access

    Returns
    ----------
    numpy.ndarray containg the atomic indices of the protein chain specified.
    
    
    """
    
    idx = pdb.top.select(f'(chainid {n}) and protein')
    
    if len(idx) == 0:
        raise ValueError('Empty chain')
        
    return idx

def select_protein_chains(pdb):
    """
    Return the atomic indices of each protein chain.
    
    Parameters
    ----------
    pdb: mdtraj object
        it is the object that load_pdb() returns
    
    Returns
    ----------
    List of arrays containing the atomic indices of each protein chain
    
    """
    n_chains = len(list(pdb.top.chains))
    protein_chains = []
    for n in range(n_chains):
        try:
            chain = select_chain(pdb, n)
        except ValueError:
            pass
        else:
            protein_chains.append(chain)
    return protein_chains


def select_ligands(pdb):
    """
    
    Parameters
    ----------
    pdb: mdtraj object
        it is the object that load_pdb() returns
    
    Returns
    ----------
    List of arrays containing the atomic indices of each ligand
    
    
    """
    idx_ligands = []
    
    for res in pdb.top.residues: 
        if not res.is_protein and not res.is_water and res.name not in COMMON_LIGANDS:
            idx_ligands.append(pdb.top.select(f"resid {res.index}"))
    
    return idx_ligands

def compute_distance_chain_ligand(pdb,protein_chains, ligands):
    """

    Parameters
    ----------
    pdb: mdtraj object
        it is the object that load_pdb() returns
    protein_chains: list of arrays
        it contains the atomic indices of the protein chains
    ligands: list of arrays
        it contains the atomic indices of the ligands of interest

    Returns
    ----------
    Dictionary wich key is the number of ligand (0-based) and wich value is a list of
    arrays, each of them containing the distances between the ligand atoms and each protein chain atoms.
    
    """
    distances = {}
    for nlig in range(len(ligands)):
        l = []
        for nchain in range(len(protein_chains)):
            pairs = list(itertools.product(protein_chains[nchain], ligands[nlig]))
            d = md.compute_distances(pdb, pairs)
            l.append(d)
        distances[nlig] = l
    return distances

def associate_ligand_to_chain(dist):
    """

    Parameters
    ----------
    dist: dictionary of lists containing numpy.ndarray
        the key is the index of the ligand and the value is a list containing all the distance numpy.ndarray

    Returns
    ----------
    A dictionary with the ligand index as key and the most closest chain as value

    """

    d = {}
    for i in range(len(dist)): # for each ligand
        l = []
        for j in range(len(dist[i])): # for each array of distances
            min1 = np.amin(dist[i][j])
            l.append(min1)

        min_dist = min(l)
        closest_chain = l.index(min_dist)
        d[i] = closest_chain

    return d

def save_chain_ligand_pdb(pdbid,pdb,ligands,protein_chains, ligand_chain):
    """
    
    Parameters
    ----------
    pdbid: string
        it is the PDB Id of the protein 
    pdb: mdtraj object
        it is the object that load_pdb() returns
    protein_chains: list of arrays
        it contains the atomic indices of the protein chains
    ligands: list of arrays
        it contains the atomic indices of the ligands of interest
    ligand_chain: dictionary
        it contains the ligand index as key and the protein chain index as value
        
    Returns
    ----------
    It saves the protein chain and the ligand into two separated pdb files in a file named <pdbid>. 
    It returns 'Well saved!' if everything went okey.
    
    """
    dir1 = f"out_{pdbid}"
    os.mkdir(dir1)
    for i in range(len(ligands)):
        lig = pdb.atom_slice(ligands[i])
        chain = pdb.atom_slice(protein_chains[ligand_chain[i]])

        dir2 = f"{dir1}/chain_lig_{i}"

        os.mkdir(dir2)

        lig.save_pdb(f"{dir2}/ligand_{i}.pdb")
        chain.save_pdb(f"{dir2}/chain_{i}.pdb")

    return("Well saved!")

