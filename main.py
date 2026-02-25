# Import library for arguments of commandline
import argparse

# Import functions of our modules
from pdb_parser import load_structure, extract_ca_atoms
from pocket_detector import detect_pockets
from output_writer import (
    write_pocket_report,
    write_pdb_with_bfactor,
    write_pymol_script
)

def main():
    """
    Main function that executes all the pipeline
    """

    # Create the arguments
    parser = argparse.ArgumentParser(
        description="Ligand Binding Site Prediction"
    )

    # Mandatory argument: PDB file
    parser.add_argument("pdb_file", help="Input PDB file")

    # Optional parameter to DBSCAN
    parser.add_argument(
        "--eps",
        type=float,
        default=8.0,
        help="DBSCAN epsilon distance"
    )

    # Optional parameter: minimum size cluster
    parser.add_argument(
        "--min_samples",
        type=int,
        default=4,
        help="Minimum residues per cluster"
    )

    # Parser arguments
    args = parser.parse_args()

    print("Loading structure...")
    structure = load_structure(args.pdb_file)
    
    print("Extracting C-alpha atoms...")
    residues, coordinates = extract_ca_atoms(structure)

    print("Detecting pockets...")
    pockets = detect_pockets(
        coordinates,
        eps=args.eps,
        min_samples=args.min_samples
    )

    print(f"{len(pockets)} pockets detected.")

    print("Writing ouputs...")
    write_pocket_report(residues, pockets)
    write_pdb_with_bfactor(structure, residues, pockets)
    write_pymol_script(residues, pockets)

    print("Done!")


# Execute only if the file is directly executed
if __name__ == "__main__":
    main()  


#PDBParser = clase de BioPython que sabe leer archivos pdb y convertirlos en un objeto estructurado
#numpy = usado para manejar coordenadas 3D como array numérico
#función load_structure = recibe la ruta de un archivo PDB y devuelve la estructura cargada
#funcion extract_ca_atoms = recorre toda la proteina, extrae solo los atomos calpha y guarda sus coordenadas
#DBSCAN es un algortimo de clustering basado en densidad, agrupa puntos que estan:
#   cerca entre si, forman regiones densas y tienen un minimo de vecinos (y marca como ruido los aislados)