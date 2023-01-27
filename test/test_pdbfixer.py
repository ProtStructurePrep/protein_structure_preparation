import pdbfixer

PDBID = '1J4R'
EXPECTED_CHAIN_IDS_REMAINING = ['A']
IDX_CHAINS_TO_REMOVE = [1,2,3,4,5,6,7,8]

def remove_chains_and_verify():
    # Create a PDBFixer instance for the given pdbid
    fixer = pdbfixer.PDBFixer(PDBID + '.pdb')
    # Remove specified chains
    fixer.removeChains(IDX_CHAINS_TO_REMOVE)
    # Check to make sure asserted chains remain
    chain_ids_remaining = [c.id for c in fixer.topology.chains()]
    assert EXPECTED_CHAIN_IDS_REMAINING == chain_ids_remaining


def main():
    remove_chains_and_verify()


if __name__ == '__main__':
    main()
