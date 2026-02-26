
import Bio 
import Bio.PDB
import sys
 

def read_pdb(pdbfile):
    try:
        pdbparser = Bio.PDB.PDBParser(QUIET=True)   # suppress PDBConstructionWarning
        struct = pdbparser.get_structure("protein", pdbfile) #get structure needs 2 arg: ID + FILE
        return struct
    except Exception as err:
        print(str(err), file=sys.stderr)
        return None 
    
def extract_atoms(struct):
    atoms = []

    # Get only first model
    model = struct[0]

    for chain in model:
        for residue in chain:

            # Ignore heteroatoms (waters, ligands, etc.)
            if residue.id[0] != " ":
                continue

            for atom in residue:

                # Ignore hydrogens
                if atom.element == "H":
                    continue

                atom_name = atom.get_name()

                atoms.append({
                    "x": atom.coord[0],
                    "y": atom.coord[1],
                    "z": atom.coord[2],
                    "atom_name": atom_name,
                    "res_name": residue.get_resname(),
                    "res_id": residue.get_id()[1],
                    "chain": chain.get_id(),
                    "alpha": 0 if atom_name == "CA" else 1
                })

    return atoms