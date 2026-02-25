# Import the DBSCAN algorithm for the spatial clustering
from sklearn.cluster import DBSCAN

# Import numpy
import numpy as np

def detect_pockets(coordinates, eps=0.8, min_samples=4):
    """
    Detect possible pockets using spatial clustering (DBSCAN).

    Parameters:
        coordinates (np.array): 3D coordinates of the residues.
        eps (float): maximum distance to consider neighbour dots.
        min_samples (int): minimum number of dots to consider a cluster.

    Returns:
        pockets (dict): dictionary with detected clusters
    """

    # Create the DBSCAN model with defined parameters
    clustering = DBSCAN(eps=eps, min_samples=min_samples)

    # Apply clustering to the coordinates
    labels = clustering.fit_predict(coordinates)

    # Dictionary where we will save the pockets
    pockets = {}

    # Iterate over each dot and its label of cluster
    for idx, label in enumerate(labels):

        # label = -1 means noise (don't belong to any cluster)
        if label == -1:
            continue

        # If it is a new cluster, we will create it in the dictionary
        if label not in pockets:
            pockets[label] = []

        # Add the index residue to their correspondent cluster
        pockets[label].append(idx)

    # Return the pockets dictionary
    return pockets