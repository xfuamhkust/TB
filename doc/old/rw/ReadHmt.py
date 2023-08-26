import numpy as np

def ReadHamiltonianReal(Name):
    
    # Read AtOrbInd
    AtOrbInd = np.load("data/" + Name + "/AtOrbInd.npy")
    
    # Read HmtRealSpace
    HmtRealSpace = np.load("data/" + Name + "/HmtRealSpace.npy")
    
    # Read HmtMatLv & HmtLv
    HmtMatLv = np.load("data/" + Name + "/HmtMatLv.npy")
    HmtLv    = np.load("data/" + Name + "/HmtLv.npy"   )
    
    # Output
    ParaHmtR = {"AtOrbInd":     AtOrbInd,
                "HmtRealSpace": HmtRealSpace,
                "HmtMatLv":     HmtMatLv,
                "HmtLv":        HmtLv,
                }
    
    return ParaHmtR