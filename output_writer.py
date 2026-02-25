# Import class to write PDB files
from Bio.PDB import PDBIO

def write_pocket_report(residues, pockets, filename="pockets.txt"):
    """
    Write a txt file with the residues of each pocket
    """

    # Open the file in write mode
    with open(filename, "w") as f:

        # Iterate over each detected pocket
        for pocket_id, indices in pockets.items():

            # Write the header
            f.write(f"Pocket {pocket_id}\n")
            f.write("-" * 20 + "\n")

            # Iterate over the pokcet residues
            for i in indices:

                # Recover the residue
                res = residues[i]

                # Obtain the chain
                chain = res.get_parent().id

                # Name of the aminoacid
                resname = res.get_resname()

                # Number of the residue
                resnum = res.get_id()[1]

                # Write the formatted information
                f.write(f"{resname} {chain}{resnum}\n")

            f.write("\n")


def write_pdb_with_bfactor(structure, residues, pockets, output_file="pockets.pdb"):
    """
    Generate a new PDB where pocket residues
    have a high (50) B-factor to visualization.
    """

    # Create a set with all the residues that belongs to pockets
    pocket_residue_indices = set()

    for indices in pockets.values():
        pocket_residue_indices.update(indices)

    # Iterate over each residue
    for i, residue in enumerate(residues):

        # Iterate over each atom of the residue
        for atom in residue:

            # If it belongs to the pocket --> high B-factor
            if i in pocket_residue_indices:
                atom.set_bfactor(50.0)

            # If not --> low B-factor
            else:
                atom.set_bfactor(0.0)

    # Save new PDB file
    io = PDBIO()
    io.set_structure(structure)
    io.save(output_file)

def write_pymol_script(residues, pockets, filename="visualize_pockets.pml"):
    """
    Generate a PyMOL script to color pockets
    """

    with open(filename, "w") as f:

        for pocket_id, indices in pockets.items():

            selection = []

            # Create logic selection of residues
            for i in indices:
                res = residues[i]
                chain = res.get_parent().id
                resnum = res.get_id()[1]

                selection.append(f"(chain {chain} and resi {resnum})")

            # Join condition with OR
            sel_string = "or".join(selection)

            # Create selection in PyMOL
            f.write(f"select pocket_{pocket_id}, {sel_string}\n")

            # Color in red
            f.write(f"color red, pocket_{pocket_id}\n")

            