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

def common_ligands(data, n=100):
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

    Example: load_pdb('/home/username/1.pdb')
    """
    if Path(pdb_name).is_file():
        return md.load_pdb(pdb_name)
    else:
        # pdb_name is a PDB ID
        link = f'http://www.rcsb.org/pdb/files/{pdb_name}.pdb'
        return md.load_pdb(link)

def select_chain(pdb, n, selection=None):
    """select atomic ids for chain with index `n`"""
    idx = pdb.top.select(f'(chainid {n}) and {selection}')
    if len(idx) == 0:
        raise ValueError('Empty chain')
    return idx

def select_protein_chains(pdb):
    """Return the atomic indices of each protein chain"""
    n_chains = len(list(pdb.top.chains))
    protein_chains = []
    for n in range(n_chains):
        try:
            chain = select_chain(pdb, n, selection='protein')
        except ValueError:
            pass
        else:
            protein_chains.append(chain)
    return protein_chains


def select_ligands(pdb, common_residues=COMMON_LIGANDS):
    """ get the atomic indices of each ligand (common residues included) """
    n_chains = len(list(pdb.top.chains))
    all_ligands_idx = []
    for n in range(n_chains):
        try:
            chain = select_chain(pdb, n, selection='not protein and not water')
        except ValueError:
            pass
        else:
            all_ligands_idx.append(chain)

    """ get atomic indices of the common residues """
    common_residues_idx = []
    for res in pdb.top.residues:
        if res.name in common_residues:
            common_residues_idx.append(pdb.top.select(f"resname == {res.name}"))

    """ select all ligands but the common residues """
    definitive_ligands_idx = []
    for i in range(len(all_ligands_idx)):
        definitive_ligands_idx.append(np.setdiff1d(all_ligands_idx[i],common_residues_idx[i]))

    return definitive_ligands_idx





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
    List of arrays, each of them containing the distances between each ligand and each protein chain
    The order of the array is as follows: (ligand_1-chain_1), (ligand_1-chain_2), (ligand_1-chain_3) ... ,
    (ligand_2-chain_1), (ligand_2-chain_2), (ligand_2-chain_3) ...

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
    os.mkdir(pdbid)
    for i in range(len(ligands)):
        lig = pdb.atom_slice(ligands[i])
        chain = pdb.atom_slice(protein_chains[ligand_chain[i]])

        directory = f"{pdbid}/chain_lig_{i}"

        os.mkdir(directory)

        lig.save_pdb(f"{directory}/ligand_{i}.pdb")
        chain.save_pdb(f"{directory}/chain_{i}.pdb")

    return("Well saved!")

