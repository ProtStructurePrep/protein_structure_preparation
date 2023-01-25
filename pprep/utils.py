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

def get_ligands(pdb):
    ligands = []

    top = pdb.top
    idx = top.select("not protein and not water")
    subset = pdb.atom_slice(idx)
    topology_ligands = subset.top

    # exclude the ligands in the 'excluded_ligands' list
    for res in topology_ligands.residues:
        if res.name not in COMMON_LIGANDS:
            ligands.append(res.name)

    return(ligands)

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
