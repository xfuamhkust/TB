import numpy as np

def WriteHamiltonianReal(ParaIn,ParaHmtR):
    
    # Parameters
    Name          = ParaIn["Name"]
    AtOrbInd      = ParaHmtR["AtOrbInd"]
    HmtRealSpace  = ParaHmtR["HmtRealSpace"]
    HmtMatLv      = ParaHmtR["HmtMatLv"]
    HmtLv         = ParaHmtR["HmtLv"]
    NumLv, NumStt = HmtMatLv[:,:,0].shape
    
    # Write AtOrbInd
    np.save("data/" + Name + "/AtOrbInd", AtOrbInd)
    
    # Write HmtRealSpace
    np.save("data/" + Name + "/HmtRealSpace", HmtRealSpace)
    
    # Write HmtMatLv & HmtLv
    np.save("data/" + Name + "/HmtMatLv", HmtMatLv)
    np.save("data/" + Name + "/HmtLv"   , HmtLv)
    
    