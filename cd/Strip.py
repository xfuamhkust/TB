import sys
import numpy as np
pi = np.pi

'''
This is a code to transfer a 2D material to quasi 1D strip.

x↓ y→
  
 □   □   □   □   □   □   □   □   □   □   □ ↑ 
                                           |
 □   □   □   □   □   □   □   □   □   □   □ W 
                                           |
 □   □   □   □   □   □   □   □   □   □   □ ↓

                 □                           □                   □
T0: □        Tx: ↓       Ty: □ → □      Txy:   ↘        Tyx:   ↗
                 □                               □           □

 H10 H00 H01
 --- --- ---
| □ | □ | □ |          | T0   Tx       |          | Ty   Txy      |
|   |   |   |          |               |          |               |
| □ | □ | □ |    H00 = | Tx+  T0   Tx  |    H01 = | Tyx  Ty   Txy |
|   |   |   |          |               |          |               |
| □ | □ | □ |          |      Tx+  T0  |          |      Tyx  Ty  |
 --- --- ---

'''

def GetStrip(ParaIn,ParaHmtR,ParaStripIn):
    
    # Parameters
    Dim        = ParaIn["Dimension"]
    Lv         = ParaIn["LatticeVector"]
    AtName      = ParaIn["AtomName"]
    AtNum      = ParaIn["AtomNumber"]
    AtLv       = ParaIn["AtomSite"]
    HmtMatLv   = ParaHmtR["HmtMatLv"]
    HmtLv      = ParaHmtR["HmtLv"]
    Cell       = ParaStripIn["Cell"]
    Dir        = ParaStripIn["Dir"]
    Wid        = ParaStripIn["Wid"]
    
    # Check dimension
    if Dim != 2:
        print("Dimension can only be two!")
        sys.exit(0)
    
    # Expand cell
    c1, c2 = Cell
    LvNew = np.copy(Lv)
    LvNew[0] *= c1
    LvNew[1] *= c2
    AtNumNew = (np.array(AtNum)*c1*c2).tolist()
    AtLvNew = []
    for AtLvi in AtLv:
        for ic1 in range(c1):
            for ic2 in range(c2):
                AtLvNew.append(Lv[0]*ic1 + Lv[1]*ic2 + AtLvi)
    AtLvNew = np.array(AtLvNew)
    AtLvNew[0] /= c1
    AtLvNew[1] /= c2
    
    # Expand Hamiltonian
    NumLv, NumStt, _ = HmtMatLv.shape
    NumSttNew = NumStt * c1 * c2
    HmtLvNew = np.unique(HmtLv // np.array([c1,c2,1]), axis = 0)
    NumLvNew = len(HmtLvNew)
    HmtMatLvNew = np.zeros((NumLvNew, NumSttNew, NumSttNew), complex)
    for iLvNew in range(NumLvNew):
        n1, n2, _ = HmtLvNew[iLvNew]
        counti = 0
        for ic1 in range(c1):
            for ic2 in range(c2):
                countj = 0
                for jc1 in range(c1):
                    for jc2 in range(c2):
                        iLv = FindIndex(np.array([n1*c1+jc1-ic1, n2*c2+jc2-ic2, 0]),HmtLv)
                        if iLv == -1: continue
                        HmtMatLvNew[iLvNew,counti*NumStt:(counti+1)*NumStt,\
                                           countj*NumStt:(countj+1)*NumStt] = HmtMatLv[iLv]
                        countj += 1
                counti += 1
    
    # Calculate T0, Tx, Ty, Txy, Tyx
    if Dir == 1:
        # x = a1, y = a2
        # Lv1D = np.linalg.norm(LvNew[1])
        IndT0  = FindIndex(np.array([ 0, 0, 0]), HmtLvNew)
        IndTx  = FindIndex(np.array([ 1, 0, 0]), HmtLvNew)
        IndTy  = FindIndex(np.array([ 0, 1, 0]), HmtLvNew)
        IndTxy = FindIndex(np.array([ 1, 1, 0]), HmtLvNew)
        IndTyx = FindIndex(np.array([-1, 1, 0]), HmtLvNew)
    elif Dir == 2:
        # x = -a2, y = a1
        # Lv1D = np.linalg.norm(LvNew[0])
        IndT0  = FindIndex(np.array([ 0, 0, 0]), HmtLvNew)
        IndTx  = FindIndex(np.array([ 0,-1, 0]), HmtLvNew)
        IndTy  = FindIndex(np.array([ 1, 0, 0]), HmtLvNew)
        IndTxy = FindIndex(np.array([ 1,-1, 0]), HmtLvNew)
        IndTyx = FindIndex(np.array([ 1, 1, 0]), HmtLvNew)
    T0  = HmtMatLvNew[IndT0]
    Tx  = HmtMatLvNew[IndTx]
    Ty  = HmtMatLvNew[IndTy]
    Txy = HmtMatLvNew[IndTxy]
    Tyx = HmtMatLvNew[IndTyx]
    
    # Calculate H00 H01
    H00 = np.kron(np.eye(Wid),T0) + np.kron(np.eye(Wid, k=1),Tx)  + np.kron(np.eye(Wid, k=-1),Tx.T.conj())
    H01 = np.kron(np.eye(Wid),Ty) + np.kron(np.eye(Wid, k=1),Txy) + np.kron(np.eye(Wid, k=-1),Tyx)
    
    # Output
    ParaStrip = {"Dir":    Dir,
                 "Wid":    Wid,
                 "Lv":     LvNew,
                 "AtName": AtName,
                 "AtNum":  AtNumNew,
                 "AtLv":   AtLvNew,
                 "Hmt00":  H00,
                 "Hmt01":  H01,
                 }
    
    return ParaStrip

def FindIndex(x,X,tol = 1e-3):
    dX = np.sum(abs(X - x), axis=1)
    min_dX = min(dX)
    if min_dX < tol:
        ind = np.where(dX == min_dX)[0][0]
        return ind
    else:
        return -1

