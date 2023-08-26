import numpy as np

def GetEigenSolution(ParaHmtK):
    
    # Parameters
    HmtKpt = ParaHmtK["HmtKpt"]
    NumKpt, NumStt = HmtKpt[:,:,0].shape
    
    # Calculate eigenvalues and eigenvectors
    EigVal = np.zeros((NumKpt,NumStt))
    EigVct = np.zeros((NumKpt,NumStt,NumStt),complex)
    for iKpt in range(NumKpt):
        EigVal[iKpt], EigVct[iKpt] = np.linalg.eigh(HmtKpt[iKpt])
    
    # Output
    ParaEig = {"EigVal": EigVal,
               "EigVct": EigVct,
               }
    return ParaEig