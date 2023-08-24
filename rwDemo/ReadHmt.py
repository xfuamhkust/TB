import numpy as np

def ReadHamiltonianReal(ParaIn):
    folderName=ParaIn["Folder"]
    # Read AtOrbInd
    AtOrbInd = np.load(folderName + "/AtOrbInd.npy")
    
    # Read HmtRealSpace
    HmtRealSpace = np.load(folderName + "/HmtRealSpace.npy")
    
    # Read HmtMatLv & HmtLv
    HmtMatLv = np.load(folderName+ "/HmtMatLv.npy")
    HmtLv    = np.load(folderName + "/HmtLv.npy"   )
    
    # Output
    ParaHmtR = {"AtOrbInd":     AtOrbInd,
                "HmtRealSpace": HmtRealSpace,
                "HmtMatLv":     HmtMatLv,
                "HmtLv":        HmtLv,
                }
    
    return ParaHmtR