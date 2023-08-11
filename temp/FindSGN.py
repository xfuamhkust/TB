import numpy as np

'''
This code is to find the space group number (SGN) of an input system.
'''

def FindSGN(AtC,AtTypeInd):
    '''
    This function is to find the SGN by comparing the matices (not symmetry itself!) in a smart way.
    Input:
        AtC: Atom sites in conventional cell.
        AtTypeInd: The type of each atom. [0,1,2,2,2] for ABO3.
    Output:
        SGN: The index of space group, 1~230.
    Comments:
        The core idea is to exclude the unsatisfied groups in each checking.
        Maxium checking times:   20
        Average checking times: ~10
    '''
    # Load SGID (1010...) & SGU (Space group matrices to check)
    SGID  = np.load("SGID.npy")
    SGU   = np.load("SGU.npy")
    NumSG,NumID = SGID.shape
    Order = np.argsort(np.sum(SGID,axis=0))[::-1]
    Ind = np.ones(NumSG,bool)
    for i in Order:
        # If the remaining groups can not be distinguished by ith matrix, continue.
        if np.var(SGID[Ind,i]) < 1e-3:
            continue
        # Use function CheckSG to check if symetry SGU[i] is in the system.
        flag = CheckSG(SGU[i],AtC,AtTypeInd)
        # Assign "False" to unsatisfied groups.
        IndInvj = np.where(1-(SGID[:,i] == flag))
        Ind[IndInvj] = False
        # If only one space group remains, return its SGN.
        if sum(Ind) == 1:
            return np.where(Ind)[0][0] + 1

def CheckSG(SymLv,AtLv,AtTypeInd,error=1e-3):
    '''
    This function is to check if one symmetry matrix is satisfied in system.
    The choice of origin is important.
    '''
    NumAt = len(AtLv)
    for iAt in range(1,NumAt+1):
        SymAt0ii = SymLv[0:3,0:3] @ AtLv[iAt-1] + SymLv[0:3,3]
        flag = 0
        for jAt in range(1,NumAt+1):
            if AtTypeInd[iAt-1] != AtTypeInd[jAt-1]:
                continue
            dLv = SymAt0ii - AtLv[jAt-1]
            if np.linalg.norm(dLv - np.round(dLv)) < error:
                flag = 1
                break
        if flag == 0:
            return False
    return True

'''
The following part is to use several material to check the validity of this code.
'''
# ABO3
AtLv = np.array([[ 0.00000000,     0.00000000,     0.00000000],
                 [ 0.50000000,     0.50000000,     0.50000000],
                 [ 0.00000000,     0.50000000,     0.50000000],
                 [ 0.50000000,     0.00000000,     0.50000000],
                 [ 0.50000000,     0.50000000,     0.00000000]])
AtTypeInd = np.array([0,1,2,2,2])
SGN = FindSGN(AtLv,AtTypeInd)
print(SGN)

# NaCl
AtLv = np.array([[ 0.00000000,     0.00000000,     0.00000000],
                 [ 0.00000000,     0.50000000,     0.50000000],
                 [ 0.50000000,     0.00000000,     0.50000000],
                 [ 0.50000000,     0.50000000,     0.00000000],
                 [ 0.50000000,     0.50000000,     0.50000000],
                 [ 0.50000000,     0.00000000,     0.00000000],
                 [ 0.00000000,     0.50000000,     0.00000000],
                 [ 0.00000000,     0.00000000,     0.50000000]])
AtTypeInd = np.array([0,0,0,0,1,1,1,1])
SGN = FindSGN(AtLv,AtTypeInd)
print(SGN)

# Graphene
AtLv = np.array([[ 0.33333333,     0.66666666,     0.00000000],
                 [ 0.66666666,     0.33333333,     0.00000000]])
AtTypeInd = np.array([0,0])
SGN = FindSGN(AtLv,AtTypeInd)
print(SGN)

# h-BN
AtLv = np.array([[ 0.33333333,     0.66666666,     0.00000000],
                 [ 0.66666666,     0.33333333,     0.00000000]])
AtTypeInd = np.array([0,1])
SGN = FindSGN(AtLv,AtTypeInd)
print(SGN)

# Si
# AtLv = np.array([[ 0.00000000,     0.00000000,     0.00000000],
#                  [ 0.00000000,     0.50000000,     0.50000000],
#                  [ 0.50000000,     0.00000000,     0.50000000],
#                  [ 0.50000000,     0.50000000,     0.00000000],
#                  [ 0.25000000,     0.25000000,     0.25000000],
#                  [ 0.75000000,     0.75000000,     0.25000000],
#                  [ 0.75000000,     0.25000000,     0.75000000],
#                  [ 0.25000000,     0.75000000,     0.75000000]])
AtLv = np.array([[ 0.00000000,     0.00000000,     0.00000000],
                 [ 0.00000000,     0.50000000,     0.50000000],
                 [ 0.50000000,     0.00000000,     0.50000000],
                 [ 0.50000000,     0.50000000,     0.00000000],
                 [ 0.25000000,     0.25000000,     0.25000000],
                 [ 0.75000000,     0.75000000,     0.25000000],
                 [ 0.75000000,     0.25000000,     0.75000000],
                 [ 0.25000000,     0.75000000,     0.75000000]])-0.125
AtTypeInd = np.array([0,0,0,0,0,0,0,0])
SGN = FindSGN(AtLv,AtTypeInd)
print(SGN)
