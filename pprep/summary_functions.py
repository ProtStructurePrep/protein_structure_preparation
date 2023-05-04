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
	
def simulate_complex(rdkit_ligand, receptor_file, 
                     equilibration_steps, num_steps, reporting_interval,
                     output_traj_pdb, output_traj_dcd):
    """Simulate the protein-ligand complex
    
    Parameters
    -----------
    rdkit_ligand:
        prepared ligand to simulate
    receptor_file: string
        prepared receptor file to simulate
    equilibration_steps: int
        steps to perform the equilibration
    num_steps: int
        number of steps to run the simulation
    reporting interval: int
        interval at whcih the reporters are added
    output_traj_pdb: string
        name of the pdb output trajectory file
    output_traj_dcd: string
        name of the dcd output trajectory file
    """
    
    ligand = rdkit_ligand  
    receptor = umd.load_prepared_receptor(receptor_file)
    system, complex_structure = umd.prepare_system(receptor, ligand)
    simulation = umd.setup_openmm_simulation(system, complex_structure)
    simulation = umd.minimize_and_equilibrate(simulation, equilibration_steps)
    umd.run_simulation(simulation, num_steps, reporting_interval, output_traj_pdb, output_traj_dcd)

