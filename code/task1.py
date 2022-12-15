import numpy as np
from readData import total_mass

scale_length = 0


def hernquist_density(r):
    return ((total_mass()/ (2*np.pi)) * (scale_length / r) * (1/np.pow((r+scale_length),3)))

