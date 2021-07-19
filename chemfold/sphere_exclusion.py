import chemfold as cf

import numpy as np
import scipy.sparse
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Sphere exclusion on the given X")
    parser.add_argument("--x", help="Molecule descriptor file, e.g., ECFPs (matrix market or numpy)", type=str, required=True)
    parser.add_argument("--out", help="Output file for the clusters (.npy)", type=str, required=True)
    parser.add_argument("--dists", nargs="+", help="Distances", type=float, required=True)

    args = parser.parse_args()
    print(args)

    print(f"Loading '{args.x}'.")
    X = cf.load_sparse(args.x).tocsr()

    print("Clustering.")
    # two step clustering, first at 0.5, then at 0.6
    cl_hier = cf.hierarchical_clustering(X, dists=args.dists)

    np.save(args.out, cl_hier)
    print(f"Saved clustering into '{args.out}'.")
