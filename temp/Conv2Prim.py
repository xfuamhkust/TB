import numpy as np

def Conv2Prim(LvC, AtC, AtTypeInd, Dim, error=1e-3):
    '''
    This function is to transfer conventional cell gauge to primitive cell gauge.
    Input:
        LvC: Lattice vectors in conventional cell.
        AtC: Atom sites in conventional cell.
        AtTypeInd: The type of each atom. [0,1,2,2,2] for ABO3.
        Dim: dimension, 2 or 3
    Output:
        LatSys: Lattice systems (7)
        BraLat: Bravais lattices (14)
        LvP: Lattice vectors in primitive cell.
    '''
    # Check dimension
    if Dim == 2:
        pass
    elif Dim == 3:
        # Get Lattice System
        LatSys = GetLatSys3D(LvC, error)
        BraLat = GetBraLat3D(LatSys, AtC, AtTypeInd)
        LvP    = GetLvP(LvC,BraLat)
    return LatSys, BraLat, LvP
        
def GetLatSys3D(LvC, error=1e-3):
    '''
    This function is to get LatSys from LvC.
    Input: LvC
    Output: LatSys
    '''
    a = np.linalg.norm(LvC[0])
    b = np.linalg.norm(LvC[1])
    c = np.linalg.norm(LvC[2])
    alpha = np.arccos(LvC[1] @ LvC[2] /b/c)*180/np.pi
    beta  = np.arccos(LvC[2] @ LvC[0] /c/a)*180/np.pi
    gamma = np.arccos(LvC[0] @ LvC[1] /a/b)*180/np.pi
    print(a,b,c)
    print(alpha,beta,gamma)
    if CkEq([alpha,beta,gamma,90], error): #c/t/o
        if CkEq([a,b,c], error): #c
            LatSys = "c"
        elif CkEq([a,b], error): #t
            LatSys = "t"
        else:             #o
            LatSys = "o"
    elif CkEq([alpha,beta,90], error) and CkEq([gamma,120], error): #hP
        if CkEq([a,b], error):   #hP
            LatSys = "hP"
        else:
            print("Wrong LatSys!")
    elif CkEq([alpha,beta,gamma], error) and 1-CkEq([alpha,beta,gamma,90], error): #hR
        if CkEq([a,b,c], error): #hR
            LatSys = "hR"
        else:
            print("Wrong LatSys!")
    elif CkEq([alpha,gamma,90], error) and 1-CkEq([beta,90], error): #m
        if 1-CkEq([a,b], error) and 1-CkEq([c,a], error) and 1-CkEq([b,c], error): #m
            LatSys = "m"
        else:
            print("Wrong LatSys!")
    elif 1-CkEq([alpha,beta], error) and 1-CkEq([gamma,alpha], error) and 1-CkEq([beta,gamma], error) and 1-CkEq([alpha,90], error) and 1-CkEq([beta,90], error) and 1-CkEq([gamma,90], error): #a
        print(1)
        if 1-CkEq([a,b], error) and 1-CkEq([c,a], error) and 1-CkEq([b,c], error): #a
            LatSys = "a"
        else:
            print("Wrong LatSys!")
    else:
        print("Wrong LatSys!")
    return LatSys
     
def GetBraLat3D(LatSys, AtC, AtTypeInd):
    '''
    This function is to get BraLat.
    Input: LatSys, AtC, AtTypeInd
    Output: BraLat
    '''
    if LatSys == "hP":
        BraLat = "hP"
    elif LatSys == "hR":
        BraLat = "hR"
    elif LatSys == "a":
        BraLat = "aP"
    elif LatSys == "c":
        if CheckBraLat(AtC,AtTypeInd,"F"):
            BraLat = "cF"
        elif CheckBraLat(AtC,AtTypeInd,"I"):
            BraLat = "cI"
        else:
            BraLat = "cP"
    elif LatSys == "t":
        if CheckBraLat(AtC,AtTypeInd,"I"):
            BraLat = "tI"
        else:
            BraLat = "tP"
    elif LatSys == "o":
        if CheckBraLat(AtC,AtTypeInd,"F"):
            BraLat = "oF"
        elif CheckBraLat(AtC,AtTypeInd,"I"):
            BraLat = "oI"
        elif CheckBraLat(AtC,AtTypeInd,"S"):
            BraLat = "oS"
        else:
            BraLat = "oP"
    elif LatSys == "m":
        if CheckBraLat(AtC,AtTypeInd,"S"):
            BraLat = "mS"
        else:
            BraLat = "mP"
    return BraLat
    
def CheckBraLat(AtC,AtTypeInd,BraLat,error=1e-3):
    '''
    This function is to check face-centered (F) / body-centered (I) / based-centered (S)
    Input: AtC,AtTypeInd,BraLat
    Output: True or False
    '''
    if BraLat == "F":
        LvP0 = np.array([[  0,1/2,1/2],[1/2,  0,1/2],[1/2,1/2,  0]])
    elif BraLat == "I":
        LvP0 = np.array([[1/2,1/2,1/2]])
    elif BraLat == "S":
        LvP0 = np.array([[1/2,1/2,  0]])
    NumAt = len(AtC); NumLv = len(LvP0)
    for iAt in range(NumAt):
        for iLv in range(NumLv):
            flag = 0
            for jAt in range(NumAt):
                if AtTypeInd[iAt] == AtTypeInd[jAt]:
                    dLv = AtC[iAt] + LvP0[iLv] - AtC[jAt]
                    if np.linalg.norm(dLv - np.round(dLv)) < error:
                        flag = 1
                        break
            if flag == 0:
                return False
    return True
                    
def GetLvP(LvC,BraLat):
    '''
    This function is to get LvP from LvC and BraLat
    Input: LvC,BraLat
    Output: LvP
    '''
    if BraLat[1] == "P" or BraLat == "hR":
        LvP = np.copy(LvC)
    elif BraLat[1] == "F":
        LvP1 = 1/2 * LvC[1] + 1/2 * LvC[2]
        LvP2 = 1/2 * LvC[2] + 1/2 * LvC[0]
        LvP3 = 1/2 * LvC[0] + 1/2 * LvC[1]
        LvP = np.array([LvP1,LvP2,LvP3])
    elif BraLat[1] == "I":
        LvP1 = (-1/2)*LvC[0] + ( 1/2) * LvC[1] + ( 1/2) * LvC[2]
        LvP2 = ( 1/2)*LvC[0] + (-1/2) * LvC[1] + ( 1/2) * LvC[2]
        LvP3 = ( 1/2)*LvC[0] + ( 1/2) * LvC[1] + (-1/2) * LvC[2]
        LvP = np.array([LvP1,LvP2,LvP3])
    elif BraLat[1] == "S":
        LvP1 = ( 1/2)*LvC[0] + ( 1/2) * LvC[1]
        LvP2 = ( 1/2)*LvC[0] + (-1/2) * LvC[1]
        LvP3 = LvC[2]
        LvP = np.array([LvP1,LvP2,LvP3])    
    return LvP
                
def CkEq(x,error=1e-3): 
    # Check Equal
    return np.std(x)/np.mean(x) < error
