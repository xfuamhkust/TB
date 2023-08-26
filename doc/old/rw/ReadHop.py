import numpy as np


def ReadHopping(FileName):
    HopValIn  = np.array([float(line.strip()) for line in open(FileName,"r").readlines()])
    return HopValIn