from scipy.io import loadmat
import numpy as np
import sys

if __name__ == "__main__":
    print('input_ants_rigid_transform input_bvecs output_bvec')

    # Loading Ants transform
    trans = loadmat(sys.argv[1])
    matrice = trans['AffineTransform_double_3_3'][:9].reshape((3,3))

    #loading bvecs
    bvecs = np.genfromtxt(sys.argv[2])
    # heuristic to have bvecs shape 3xN
    if bvecs.shape[0] != 3:
        bvecs = bvecs.T

    # Rotating bvecs
    newbvecs = np.dot(matrice, bvecs)

    # saving bvecs
    np.savetxt(sys.argv[3], newbvecs)
