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
    Get the pdb residues that are fond in more than `n` proteins in the PDB.
    
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
    excluded_residues = []
    threshold = n # check threshold

    for line in data:
        line = line.split()
        res = line[0]
        n_prot = len(line[1:])
        if n_prot > threshold:
            excluded_residues.append(res)

    return excluded_residues

def load_pdb(pdb_name):
    """
    Parse a pdb file/fetch a structure from the PDB.
    
    Parameters
    ----------
    pdb_name: string
        it can be either the pdb file path or the pdb id of the protein
    
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
    Select atomic ids for protein chain with index `n`.
    
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
    Select the atomic indices of each protein chain.
    
    Parameters
    ----------
    pdb: mdtraj object
        it is the object that load_pdb() returns
    
    Returns
    ----------
    List of numpy.ndarray containing the atomic indices of each protein chain
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

def select_ligands_name(pdb):
    """
    Get the ligand names.
    
    Parameters
    ----------
    pdb: mdtraj object
        it is the object that load_pdb() returns
    
    Returns
    ----------
    List. `ligand_names`, containing the names of the ligands
    """
    ligand_names = []
    
    for res in pdb.top.residues: 
        if not res.is_protein and not res.is_water:
            if res.name not in COMMON_LIGANDS:
                ligand_names.append(res.name)
    
    return ligand_names
    
def select_ligands(pdb):
    """
    Get the atomic indices of the ligands.
    
    Parameters
    ----------
    pdb: mdtraj object
        it is the object that load_pdb() returns
    
    Returns
    ----------
    List. `ligand_idx`, containig an array of the ligand atomic indices for each ligand.
    """
    ligand_idx = []
    
    for res in pdb.top.residues: 
        if not res.is_protein and not res.is_water:
            if res.name not in COMMON_LIGANDS:
                ligand_idx.append(pdb.top.select(f"resid {res.index}"))
    
    return ligand_idx

def compute_distance_chain_ligand(pdb, protein_chains, ligands):
    """
    Compute the distance between the ligand atoms and each protein chain atoms.
    
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
    Dictionary. Its key is the number of ligand (0-based) and its value is a list of
    arrays, each of them containing the distances between the ligand atoms and each protein chain atoms.
    """
    distances = {}
    
    for nlig in range(len(ligands)):
        l = []
        for nchain in range(len(protein_chains)):
            pairs = list(itertools.product(protein_chains[nchain], ligands[nlig]))
            l.append(md.compute_distances(pdb, pairs))
        distances[nlig] = l
        
    return distances

def associate_ligand_to_chain(dist):
    """
    Associate each ligand to a protein chain by considering the distance at which the ligand is 
    from the protein chain.
    
    Parameters
    ----------
    dist: dictionary of lists containing numpy.ndarray
        the key is the index of the ligand and the value is a list containing all the distance numpy.ndarray

    Returns
    ----------
    A dictionary with the ligand index as key and the most closest chain index as value.
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

def save_chain_ligand_pdb(pdbid, pdb, ligands, protein_chains, ligand_chain, output_directory):
    """
    Create a directory for the protein called `out_pdbid` with a directory for each ligand. In each ligand
    directory you can find the pdb of the ligand and the pdb of its associated chain.
    
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
        it contains the ligand index as key and its the index of its associated chain as value
    output_directory: str
        Directory in wich the files will be saved.
    """
    if not os.path.exists(output_directory):
        os.mkdir(output_directory)
  
    dir1 = f"{output_directory}/out_{pdbid}" # create the protein directory
    os.mkdir(dir1)
    
    for i in range(len(ligands)):
        lig = pdb.atom_slice(ligands[i])
        chain = pdb.atom_slice(protein_chains[ligand_chain[i]])

        dir2 = f"{dir1}/chain_lig_{i}" # create the ligand directory

        os.mkdir(dir2)

        lig.save_pdb(f"{dir2}/ligand_{i}.pdb")
        chain.save_pdb(f"{dir2}/chain_{i}.pdb")

    print("Well saved!")
