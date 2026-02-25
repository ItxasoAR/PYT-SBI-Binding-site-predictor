# Import parser of PDB from BioPython
from Bio.PDB import PDBParser

# Import numpy to handle numeric arrays
import numpy as np


def load_structure(pdb_file):
    """
    Load a PDB structure from a file.
    
    Parameters:
        pdb_file (str): Path to PDB file.
    
    Returns:
        structure: object structure of BioPython.
    """

    # Create the parser. QUIET=True avoid innecessary messages in the console.
    parser = PDBParser(QUIET=True)

    # Read the PDB file and convert it in a structure object
    structure = parser.get_structure("protein", pdb_file)

    # Return the load structure
    return structure

def extract_ca_atoms(structure):
    """
    Extract the C-alpha atoms fromm the protein.

    Parameters:
        structure: object structure of Biopython.

    Return:
        residues (list): residues list which contain the C-alpha.
        coordinates (np.array): array with 3D coordinates for each C-alpha.
    """

    # List where we will save the residues
    residues = []

    # List where we will save the 3D coordinates
    coordinates = []

    # Iterate over each module (normally there is only one)
    for model in structure:

        # Iterate over each chain (A, B, etc.)
        for chain in model:

            # Iterate over each residue from the chain 
            for residue in chain:

                # Check if the residue has C-alpha atom
                if "CA" in residue:

                    # Save the C-alpha atom
                    ca = residue["CA"]

                    # Save the full residue
                    residues.append(residue)

                    # Save its 3D coordinates 
                    coordinates.append(ca.get_coord())
    
    # Convert the coordinates list to a numpy array
    return residues, np.array(coordinates)

