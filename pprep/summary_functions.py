import pprep.utils_mdanalysis as umda
import pprep.utils_rdkit as urk
import pprep.utils_pdbfixer as updbfix
import pprep.MDutils as umd


def parse_complex(pdb_name, output_ligand, output_prot):
    """Select the receptor and the ligand from a protein ligand complex
    
    Parameters
    ----------
    pdb_name: string
        it can be either the pdb file path or the pdb id of the protein
    """
    u = umda.load_pdb(pdb_name)
    u, ligand_name = umda.find_potential_ligand(u)
    ligand = umda.select_ligand(u, ligand_name, output_ligand)
    protein = umda.select_protein(u, output_prot)
    return ligand, protein

def fix_complex(pdb_file, resname, smiles, protein_pdb, output_prepared_prot):
	"""Fix the ligand and the protein for further simulations """
	rdkit_ligand = urk.prepare_ligand(pdb_file, resname, smiles)
	pdbfix_prot = updbfix.prepare_protein(protein_pdb, output_prepared_prot)
	return rdkit_ligand, pdbfix_prot
	


