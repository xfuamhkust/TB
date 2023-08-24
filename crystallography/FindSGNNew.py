import numpy as np
# import sys


'''
This code is to find the space group number (SGN) of an input system.
'''

"""
       @Author: FXZ
       """

def FindSGN(AtC,AtTypeInd):
    '''
    Based on FindSGN, This function is added with an origin-finding process.
    '''
    # Load SGID (1010...) & SGM (Space group matrices to check)
    SGID  = np.load("crystallography/SGID.npy")
    
    # Get number of elements for each space group
    NumElm = np.sum(SGID,axis=1) #   1   2   3   4   6   8   9  12  16  18  24  32  36  48  96 192
    
    # Get Origin, containing [0,0,0], AtC, and centers of bonds AtC--(Lv+Atc).
    # ([n1,n2,n3] + Ati + Atj)/2 = [n1/2,n2/2,n3/2] + (Ati+Atj)/2
    # ni/2 = -1/2 0 1/2 1; Total number is 4^3 * NumAt^2
    error = 1e-3
    NumAt = len(AtC)
    n0 = [-1/2, 0, 1/2, 1]
    AtC0 = np.mod(AtC,1)
    NumOrg = 4**3 * NumAt**2 + 1
    # Get all origins
    Origin0 = np.zeros((NumOrg,3))
    count = 1
    for x in n0:
        for y in n0:
            for z in n0:
                Lv_2 = np.array([x,y,z])/2
                for iAt in range(1,NumAt+1):
                    for jAt in range(1,NumAt+1):
                        Origin0[count] = 0.5*(AtC0[iAt-1] + AtC0[jAt-1]) + Lv_2
                        count += 1
    Origin0 = np.mod(Origin0,1)
    # Reduce repeat origins
    for i in range(NumOrg):
        if np.sum(np.linalg.norm(Origin0 - Origin0[i],axis=1) < error) > 1:
            Origin0[i] = [0.,0.,0.]
    Ind0 = np.where(np.linalg.norm(Origin0,axis=1) != 0.)[0]
    Origin = np.zeros((len(Ind0) + 1,3))
    Origin[1:] = Origin0[Ind0]
    NumOrg = len(Origin)
    # Check all origins
    SGN       = np.zeros(NumOrg,int)
    NumElmOrg = np.zeros(NumOrg,int)
    for i in range(NumOrg):
        SGN[i] = FindSGN0(AtC-Origin[i],AtTypeInd)
        if SGN[i]:
            NumElmOrg[i] = NumElm[SGN[i]-1]
        else:
            NumElmOrg[i] = 0
    # Choose the space group with highest symmetry (shortest norm)
    MaxNumElm = np.max(NumElmOrg)
    IndMax = np.where(NumElmOrg == MaxNumElm)[0]
    OrgNorm = np.linalg.norm(Origin[IndMax],axis=1)
    IndMin = np.where(OrgNorm == np.min(OrgNorm))[0][0]
    IndOrg = IndMax[IndMin]
    return SGN[IndOrg], Origin[IndOrg]
        
    


def FindSGN0(AtC,AtTypeInd):
    '''
    This function is to find the SGN by comparing the matices (not symmetry itself!) in a smart way.
    Input:
        AtC: Atom sites in conventional cell.
        AtTypeInd: The type of each atom. [0,1,2,2,2] for ABO3.
    Output:
        SGN: The index of space group, 1~230.
    Comments:
        ???
    '''
    # Load SGID (1010...) & SGM (Space group matrices to check)
    SGID  = np.load("crystallography/SGID.npy")
    SGM   = np.load("crystallography/SGM.npy")
    NumSG,NumID = SGID.shape
    ID = np.zeros(NumID,int)
    for i in range(NumID):
        # Use function CheckSG to check if symetry SGU[i] is in the system.
        ID[i] = CheckSG(SGM[i],AtC,AtTypeInd)
    Ind = np.where(np.sum(abs(SGID - ID),axis=1)==0)[0]
    if len(Ind) == 1:
        return Ind[0] + 1
    else:
        # print("Wrong Space Group!!!")
        return 0
    

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

def swapVec(vec,i,j):
    """
    i: ind
    j: ind
    """
    vec[i],vec[j]=vec[j],vec[i]
    return vec

def swapList(vecList,i,j):
    """
    swap ith and jth coordinates for each vec in vecList
    """
    retList=[swapVec(vec,i,j) for vec in vecList]
    return np.array(retList)
def inverseVec(vec,i):
    """
    inverse element at ith position of vec
    """
    vec[i]=-vec[i]
    return vec
def inverseList(vecList,i):
    """
    inverse each vector
    """
    retList=[inverseVec(vec,i) for vec in vecList]
    return np.array(retList)

# # ABO3
# AtLv = np.array([[ 0.00000000,     0.00000000,     0.00000000],
#                  [ 0.50000000,     0.50000000,     0.50000000],
#                  [ 0.00000000,     0.50000000,     0.50000000],
#                  [ 0.50000000,     0.00000000,     0.50000000],
#                  [ 0.50000000,     0.50000000,     0.00000000]])
# AtTypeInd = np.array([0,1,2,2,2])
# #swap x,y
# # AtLv=swapList(AtLv,0,1)
# #swap x,z
# # AtLv=swapList(AtLv,0,2)
# #swap y,z
# # AtLv=swapList(AtLv,1,2)
#
# # SGN, Org = FindSGN(AtLv,AtTypeInd)
# # print("ABO3 group number="+str(SGN)+", origin="+str(Org))
#
# # NaCl
# AtLv = np.array([[ 0.00000000,     0.00000000,     0.00000000],
#                  [ 0.00000000,     0.50000000,     0.50000000],
#                  [ 0.50000000,     0.00000000,     0.50000000],
#                  [ 0.50000000,     0.50000000,     0.00000000],
#                  [ 0.50000000,     0.50000000,     0.50000000],
#                  [ 0.50000000,     0.00000000,     0.00000000],
#                  [ 0.00000000,     0.50000000,     0.00000000],
#                  [ 0.00000000,     0.00000000,     0.50000000]])
# AtTypeInd = np.array([0,0,0,0,1,1,1,1])
# # SGN, Org = FindSGN(AtLv,AtTypeInd)
# # print("NaCl group number="+str(SGN)+", origin="+str(Org))
#
# # Graphene
# AtLv = np.array([[ 0.33333333,     0.66666666,     0.00000000],
#                  [ 0.66666666,     0.33333333,     0.00000000]])
# AtTypeInd = np.array([0,0])
# SGN, Org = FindSGN(AtLv,AtTypeInd)
# print("graphene group number="+str(SGN)+", origin="+str(Org))
#
# # h-BN
# AtLv = np.array([[ 0.33333333,     0.66666666,     0.00000000],
#                  [ 0.66666666,     0.33333333,     0.00000000]])
# AtTypeInd = np.array([0,1])
# # SGN, Org = FindSGN(AtLv,AtTypeInd)
# # print("h-BN group number="+str(SGN)+", origin="+str(Org))
#
# # Si
# # AtLv = np.array([[ 0.00000000,     0.00000000,     0.00000000],
# #                  [ 0.00000000,     0.50000000,     0.50000000],
# #                  [ 0.50000000,     0.00000000,     0.50000000],
# #                  [ 0.50000000,     0.50000000,     0.00000000],
# #                  [ 0.25000000,     0.25000000,     0.25000000],
# #                  [ 0.75000000,     0.75000000,     0.25000000],
# #                  [ 0.75000000,     0.25000000,     0.75000000],
# #                  [ 0.25000000,     0.75000000,     0.75000000]])
# AtLv = np.array([[ 0.00000000,     0.00000000,     0.00000000],
#                  [ 0.00000000,     0.50000000,     0.50000000],
#                  [ 0.50000000,     0.00000000,     0.50000000],
#                  [ 0.50000000,     0.50000000,     0.00000000],
#                  [ 0.25000000,     0.25000000,     0.25000000],
#                  [ 0.75000000,     0.75000000,     0.25000000],
#                  [ 0.75000000,     0.25000000,     0.75000000],
#                  [ 0.25000000,     0.75000000,     0.75000000]])#-0.125
# AtTypeInd = np.array([0,0,0,0,0,0,0,0])
# #swap x,y
# # AtLv=swapList(AtLv,0,1)
# #swap x,z
# # AtLv=swapList(AtLv,0,2)
# #swap y,z
# # AtLv=swapList(AtLv,1,2)
# #inverse x
# # AtLv=inverseList(AtLv,0)
# #inverse y
# # AtLv=inverseList(AtLv,1)
# #inverse z
# # AtLv=inverseList(AtLv,2)
# # SGN, Org = FindSGN(AtLv,AtTypeInd)
# # print("Si group number="+str(SGN)+", origin="+str(Org))
