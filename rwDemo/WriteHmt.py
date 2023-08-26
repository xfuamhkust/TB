import numpy as np

def WriteHamiltonianReal(ParaIn,ParaHmtR):
    
    # Parameters
    Name          = ParaIn["Name"]
    AtOrbInd      = ParaHmtR["AtOrbInd"]
    HmtRealSpace  = ParaHmtR["HmtRealSpace"]
    HmtMatLv      = ParaHmtR["HmtMatLv"]
    HmtLv         = ParaHmtR["HmtLv"]
    folderName=ParaIn["Folder"]
    NumLv, NumStt = HmtMatLv[:,:,0].shape
    
    # Write AtOrbInd
    np.save(folderName+ "/AtOrbInd", AtOrbInd)
    
    # Write HmtRealSpace
    np.save(folderName+ "/HmtRealSpace", HmtRealSpace)
    
    # Write HmtMatLv & HmtLv
    np.save(folderName + "/HmtMatLv", HmtMatLv)
    np.save(folderName+ "/HmtLv"   , HmtLv)
    
    