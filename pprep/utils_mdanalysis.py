import MDAnalysis as mda
import nglview as nv
from pathlib import Path

COMMON_RESIDUES = '''
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


def load_pdb(pdb_name):
    """Parse a pdb file/fetch a structure from the PDB.

    Parameters
    ----------
    pdb_name: string
        it can be either the pdb file path or the pdb id of the protein

    Returns
    ----------
    MDAnalysis object.

    Example: load_pdb('/home/username/1.pdb')
    """

    if Path(pdb_name).is_file():
        return mda.Universe(pdb_name)
    else:
        # pdb_name is a PDB ID
        return mda.fetch_mmtf(pdb_name)


def show_nv(u):
    """Show MDAnalysis object
    
    Parameters
    ----------
    u: MDAnalysis.core.universe.Universe
        MDAnalysis object to show
    """
    view = nv.show_mdanalysis(u)
    return view

def find_potential_ligand(u):
    """Find the name of the potential ligand and remove ions and residues from the complex
    
    Parameters
    ----------
    u: MDAnalysis.core.universe.Universe
        Protein-ligand complex object
        
    Returns
    ----------
    u: MDAnalysis.core.universe.Universe
        Protein-ligand complex object without any unwanted molecule
    potential_ligand: str
        Name of the potential ligand of the protein-ligand complex
    """
    not_protein = set(u.atoms.select_atoms("not protein and not resname HOH").residues.resnames.tolist())
    for res in not_protein:
        if res in COMMON_RESIDUES:
            u = u.atoms.select_atoms(f'not resname {res}')
        else:
            potential_ligand = res
            return u, potential_ligand
    return u, None

def select_ligand(u, ligand_name, output_file_name):
    """Select the ligand to simulate and save it in a pdb file
    
    Parameters
    ----------
    u: MDAnalysis.core.universe.Universe
        Protein-ligand complex object
    ligand_name: string
        Name of the ligand to simulate
    output_file_name: string
        Name of the output file pdb where the ligand structure will be stored
        
    Returns
    ----------
    lig: MDAnalysis.core.universe.Universe
        Ligand object
    """
    lig = u.atoms.select_atoms(f'resname {ligand_name}')
    lig.write(output_file_name)
    
    return lig

def display_chains(u):
    """Display the chain ids of the complex chains
    
    Parameters
    ----------
    u: MDAnalysis.core.universe.Universe
        Protein-ligand complex object
        
    Returns
    ----------
    chains: list
        List of all the chains present in the protein-ligand complex
    """
    chains = list(set(u.segments.segids.tolist()))
    return chains

def select_protein(u, output_file_name, leave_cristalographic_waters=False, chainid=None):
    """Select the receptor to simulate and save it in a pdb file
    
    Parameters
    ----------
    u: MDAnalysis.core.universe.Universe
        Protein-ligand complex object
    output_file_name: string
        Name of the output file pdb where the ligand structure will be stored
    chainid: list (default : None)
        List containing the chains of the protein-ligand complex that you want to simulate
           
    Returns
    ----------
    receptor: MDAnalysis.core.universe.Universe
        Receptor object
    """
    if leave_cristalographic_waters:
        receptor = u.atoms.select_atoms("protein or resname HOH")
    else:
        receptor = u.atoms.select_atoms("protein")
        
    if not chainid:
        receptor.write(output_file_name)
        return receptor
    else:
        sep = ' and segid '
        string = 'segid ' + sep.join(chainid)
        receptor = receptor.atoms.select_atoms(string)
        receptor.write(output_file_name)
        return receptor
