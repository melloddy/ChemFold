import scipy.io
import numpy as np

def load_sparse(filename):
    """Loads sparse from Matrix market or Numpy .npy file."""
    if filename is None:
        return None
    if filename.endswith('.mtx'):
        return scipy.io.mmread(filename).tocsr()
    elif filename.endswith('.npy'):
        return np.load(filename, allow_pickle=True).item().tocsr()
    raise ValueError(f"Loading '{filename}' failed. It must have a suffix '.mtx' or '.npy'.")

