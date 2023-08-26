import numpy as np

def GetSpaceGroup(ParaIn):
    SGN  = ParaIn["SpaceGroupNumber"]
    LvSG = ParaIn["SpaceGroupLatticeVector"]
    Lv   = ParaIn["LatticeVector"]
    SymLvSG   = GetSymLvSG(SGN)
    SymXyz, SymXyzt = GetSymXyz(SymLvSG,LvSG)
    SymLv   = GetSymLv(SymXyzt,Lv)
    SymOrb  = GetSymOrb(SymXyz)
    SymSpn  = np.array([GetSymSpin(SymXyz[i]) for i in range(len(SymLvSG))],"complex")
    ParaSym = {"SymLvSG": SymLvSG,
               "SymXyz":  SymXyz,
               "SymXyzt": SymXyzt,
               "SymLv":   SymLv,
               "SymOrb":  SymOrb,
               "SymSpn":  SymSpn,
               }
    
    return ParaSym

def GetSymLvSG(SGN):
    FileName = "data/sys/SpGrpMat.txt"
    f = [line.strip() for line in open(FileName,"r").readlines()]
    for iff in range(len(f)):
        fispl = f[iff].split(" ")
        if fispl[0] == "_" + "%d"%SGN + "_":
            NumMat = int(fispl[1])
            SGM = np.zeros((NumMat,3,4))
            for iMat in range(NumMat):
                fiispl = f[iff+iMat+1].split(" ")
                MatLine = []
                for iNum in range(12):
                    if fiispl[iNum].find("/") + 1:
                        [stri1,stri2] = fiispl[iNum].split("/")
                        MatLine.append(float(stri1)/float(stri2))
                    else:
                        MatLine.append(int(fiispl[iNum]))
                SGM[iMat] = np.array(MatLine).reshape((3,4))
            break
    return SGM

def GetSymXyz(SymLvSG,LvSG):
    LvSGT = LvSG.T; LvSGTI = np.linalg.inv(LvSGT); NumSym = len(SymLvSG)
    SymXyzt = np.zeros((NumSym,3,4))
    for i in range(NumSym):
        SymXyzt[i,:,:3] = LvSGT @ SymLvSG[i,:,:3] @ LvSGTI
        SymXyzt[i,:, 3] =         SymLvSG[i,:, 3] @ LvSG
    SymXyz  = SymXyzt[:,:,:3]
    return SymXyz, SymXyzt

def GetSymLv(SymXyzt,Lv):
    LvT = Lv.T;LvI = np.linalg.inv(Lv);LvTI = np.linalg.inv(LvT); NumSym = len(SymXyzt)
    SymLv = np.zeros((NumSym,3,4))
    for i in range(NumSym):
        SymLv[i,:,:3] = LvTI @ SymXyzt[i,:,:3] @ LvT
        SymLv[i,:, 3] =        SymXyzt[i,:, 3] @ LvI
    return SymLv

def GetSymOrb(SymXyz):
    NumSym = len(SymXyz)
    # S: s
    SymS = np.ones((NumSym,1,1))
    # P: px,py,pz
    SymP = np.copy(SymXyz)
    # D: d_xy,d_yz,d_zx,d_(x^2-y^2 ),d_(3z^2-r^2 )
    SymD = np.array([GetSymD(SymXyz[i]) for i in range(NumSym)])
    # F: fz3,fxz2,fyz2,fxyz,fz(x2-y2),fx(x2-3y2),fy(3x2-y2)
    SymF = np.array([GetSymF(SymXyz[i]) for i in range(NumSym)])
    return [SymS,SymP,SymD,SymF]

def GetSymD(R):
    [[R_11, R_12, R_13], [R_21, R_22, R_23], [R_31, R_32, R_33]] = R
    RD = np.zeros((5,5))
    sr3 = np.sqrt(3)
    #
    RD[0,0] = R_11*R_22+R_12*R_21
    RD[0,1] = R_21*R_32+R_22*R_31
    RD[0,2] = R_11*R_32+R_12*R_31
    RD[0,3] = 2*R_11*R_12+R_31*R_32
    RD[0,4] = sr3*R_31*R_32
    #
    RD[1,0] = R_12*R_23+R_13*R_22
    RD[1,1] = R_22*R_33+R_23*R_32
    RD[1,2] = R_12*R_33+R_13*R_32
    RD[1,3] = 2*R_12*R_13+R_32*R_33
    RD[1,4] = sr3*R_32*R_33
    #
    RD[2,0] = R_11*R_23+R_13*R_21
    RD[2,1] = R_21*R_33+R_23*R_31
    RD[2,2] = R_11*R_33+R_13*R_31
    RD[2,3] = 2*R_11*R_13+R_31*R_33
    RD[2,4] = sr3*R_31*R_33
    #
    RD[3,0] = R_11*R_21-R_12*R_22
    RD[3,1] = R_21*R_31-R_22*R_32
    RD[3,2] = R_11*R_31-R_12*R_32
    RD[3,3] = (R_11**2-R_12**2 )+1/2*(R_31**2-R_32**2 )
    RD[3,4] = sr3/2*(R_31**2-R_32**2 )
    #
    RD[4,0] = 1/sr3*(2*R_13*R_23-R_11*R_21-R_12*R_22)
    RD[4,1] = 1/sr3*(2*R_23*R_33-R_21*R_31-R_22*R_32)
    RD[4,2] = 1/sr3*(2*R_13*R_33-R_11*R_31-R_12*R_32)
    RD[4,3] = 1/sr3*(2*R_13**2-R_11**2-R_12**2 )+1/sr3/2*(2*R_33**2-R_31**2-R_32**2 )
    RD[4,4] = 1/2*(2*R_33**2-R_31**2-R_32**2 )
    
    return RD.T

def GetSymSpin(R,T=0):
    # This part refers to Ziyu's TB code.
    RS = np.zeros((2,2),"complex")
    # R -> axis & theta
    tol = 1e-6
    if (np.linalg.det(R)+1) < tol:# inversion
        R = -R
    theta = np.arccos((np.trace(R)-1)/2)
    if abs(theta) < tol:
        vec_x = 0
        vec_y = 0
        vec_z = 1
    elif abs(theta-np.pi) < tol:
        vec_x = np.sqrt((R[0, 0] + 1) / 2)
        vec_y = np.sqrt((R[1, 1] + 1) / 2)
        vec_z = np.sqrt((R[2, 2] + 1) / 2)
    else:
        vec_x = (R[2, 1] - R[1, 2]) / (2*np.sin(theta))
        vec_y = (R[0, 2] - R[2, 0]) / (2*np.sin(theta))
        vec_z = (R[1, 0] - R[0, 1]) / (2*np.sin(theta))
    # axis & theta -> RS
    sigma_0 = np.eye(2)
    sigma_x = np.array([[0,   1], [1 ,  0]])
    sigma_y = np.array([[0, -1j], [1j,  0]])
    sigma_z = np.array([[1,   0], [0 , -1]])
    RS = np.cos(theta / 2) * sigma_0 +\
         np.sin(theta / 2) * 1j * \
         (vec_x * sigma_x + vec_y * sigma_y + vec_z * sigma_z)
    if T:
        RS = -1j * np.matmul(RS,sigma_y)
    
    return RS

def GetSymF(R):
    sr3 = np.sqrt(3); sr5 = np.sqrt(5); sr15 = np.sqrt(15)
    x1x2x3 = np.array([[ 1, 1, 1], # x3
                       [ 2, 2, 2], # y3
                       [ 3, 3, 3], # z3
                       [ 1, 1, 2], # x2y
                       [ 1, 2, 2], # xy2
                       [ 1, 1, 3], # x2z
                       [ 1, 3, 3], # xz2
                       [ 2, 2, 3], # y2z
                       [ 2, 3, 3], # yz2
                       [ 1, 2, 3]  # xyz
                       ],int)
    # 10*10, x3~xyz, a~j
    Rx1x2x3 = np.zeros((10,10))
    for i in range(10):
        n1,n2,n3 = x1x2x3[i]
        Rx1x2x3[i,0] = R[ 1- 1,n1- 1] * R[ 1- 1,n2- 1] * R[ 1- 1,n3- 1] # a
        Rx1x2x3[i,1] = R[ 2- 1,n1- 1] * R[ 2- 1,n2- 1] * R[ 2- 1,n3- 1] # b
        Rx1x2x3[i,2] = R[ 3- 1,n1- 1] * R[ 3- 1,n2- 1] * R[ 3- 1,n3- 1] # c
        Rx1x2x3[i,3] = R[ 1- 1,n1- 1] * R[ 1- 1,n2- 1] * R[ 2- 1,n3- 1]\
                     + R[ 1- 1,n1- 1] * R[ 2- 1,n2- 1] * R[ 1- 1,n3- 1]\
                     + R[ 2- 1,n1- 1] * R[ 1- 1,n2- 1] * R[ 1- 1,n3- 1] # d
        Rx1x2x3[i,4] = R[ 1- 1,n1- 1] * R[ 2- 1,n2- 1] * R[ 2- 1,n3- 1]\
                     + R[ 2- 1,n1- 1] * R[ 2- 1,n2- 1] * R[ 1- 1,n3- 1]\
                     + R[ 2- 1,n1- 1] * R[ 1- 1,n2- 1] * R[ 2- 1,n3- 1] # e
        Rx1x2x3[i,5] = R[ 1- 1,n1- 1] * R[ 1- 1,n2- 1] * R[ 3- 1,n3- 1]\
                     + R[ 1- 1,n1- 1] * R[ 3- 1,n2- 1] * R[ 1- 1,n3- 1]\
                     + R[ 3- 1,n1- 1] * R[ 1- 1,n2- 1] * R[ 1- 1,n3- 1] # f
        Rx1x2x3[i,6] = R[ 1- 1,n1- 1] * R[ 3- 1,n2- 1] * R[ 3- 1,n3- 1]\
                     + R[ 3- 1,n1- 1] * R[ 3- 1,n2- 1] * R[ 1- 1,n3- 1]\
                     + R[ 3- 1,n1- 1] * R[ 1- 1,n2- 1] * R[ 3- 1,n3- 1] # g
        Rx1x2x3[i,7] = R[ 2- 1,n1- 1] * R[ 2- 1,n2- 1] * R[ 3- 1,n3- 1]\
                     + R[ 2- 1,n1- 1] * R[ 3- 1,n2- 1] * R[ 2- 1,n3- 1]\
                     + R[ 3- 1,n1- 1] * R[ 2- 1,n2- 1] * R[ 2- 1,n3- 1] # h
        Rx1x2x3[i,8] = R[ 2- 1,n1- 1] * R[ 3- 1,n2- 1] * R[ 3- 1,n3- 1]\
                     + R[ 3- 1,n1- 1] * R[ 3- 1,n2- 1] * R[ 2- 1,n3- 1]\
                     + R[ 3- 1,n1- 1] * R[ 2- 1,n2- 1] * R[ 3- 1,n3- 1] # i
        Rx1x2x3[i,9] = R[ 1- 1,n1- 1] * R[ 2- 1,n2- 1] * R[ 3- 1,n3- 1]\
                     + R[ 1- 1,n1- 1] * R[ 3- 1,n2- 1] * R[ 2- 1,n3- 1]\
                     + R[ 2- 1,n1- 1] * R[ 1- 1,n2- 1] * R[ 3- 1,n3- 1]\
                     + R[ 2- 1,n1- 1] * R[ 3- 1,n2- 1] * R[ 1- 1,n3- 1]\
                     + R[ 3- 1,n1- 1] * R[ 1- 1,n2- 1] * R[ 2- 1,n3- 1]\
                     + R[ 3- 1,n1- 1] * R[ 2- 1,n2- 1] * R[ 1- 1,n3- 1] # j
    # 7*10, fz3~fy(3x2-y2), x3~xyz
    '''           [     "x3",     "y3",     "z3",    "x2y",    "xy2",    "x2z",    "xz2",    "y2z",    "yz2",    "xyz"] '''
    F = np.array([[        0,        0,   1/sr15,        0,        0,-3/2/sr15,        0,-3/2/sr15,        0,        0], # fz3
                  [ -1/2/sr5,        0,        0,        0, -1/2/sr5,        0,    2/sr5,        0,        0,        0], # fxz2
                  [        0, -1/2/sr5,        0, -1/2/sr5,        0,        0,        0,        0,    2/sr5,        0], # fyz2
                  [        0,        0,        0,        0,        0,        0,        0,        0,        0,        1], # fxyz
                  [        0,        0,        0,        0,        0,      1/2,        0,     -1/2,        0,        0], # fz(x2-y2)
                  [  1/2/sr3,        0,        0,        0,   -sr3/2,        0,        0,        0,        0,        0], # fx(x2-3y2)
                  [        0, -1/2/sr3,        0,    sr3/2,        0,        0,        0,        0,        0,        0]  # fy(3x2-y2)
                  ])
    FR = F @ Rx1x2x3 # 7*10, fz3~fy(3x2-y2), a~j
    # 7*10, fz3~fy(3x2-y2), a~j
    '''            [    "a",    "b",    "c",    "d",    "e",    "f",    "g",    "h",    "i",    "j"] '''
    CF = np.array([[      0,      0,   sr15,      0,      0,      0,      0,      0,      0,      0], # fz3
                   [      0,      0,      0,      0,      0,      0,  sr5/2,      0,      0,      0], # fxz2
                   [      0,      0,      0,      0,      0,      0,      0,      0,  sr5/2,      0], # fyz2
                   [      0,      0,      0,      0,      0,      0,      0,      0,      0,      1], # fxyz
                   [      0,      0,      3,      0,      0,      2,      0,      0,      0,      0], # fz(x2-y2)
                   [  2*sr3,      0,      0,      0,      0,      0,  sr3/2,      0,      0,      0], # fx(x2-3y2)
                   [      0, -2*sr3,      0,      0,      0,      0,      0,      0, -sr3/2,      0]  # fy(3x2-y2)
                   ]) 
    RF = FR @ CF.T
    
    return RF.T
    
def GetSymD_new(R):
    sr3 = np.sqrt(3)
    x1x2 = np.array([[ 1, 1], # x2
                     [ 2, 2], # y2
                     [ 3, 3], # z2
                     [ 1, 2], # xy
                     [ 1, 3], # xz
                     [ 2, 3], # yz
                     ],int)
    # 6*6, x2~yz, a~f
    Rx1x2 = np.zeros((6,6))
    for i in range(6):
        n1,n2 = x1x2[i]
        Rx1x2[i,0] = R[ 1- 1,n1- 1] * R[ 1- 1,n2- 1] # a
        Rx1x2[i,1] = R[ 2- 1,n1- 1] * R[ 2- 1,n2- 1] # b
        Rx1x2[i,2] = R[ 3- 1,n1- 1] * R[ 3- 1,n2- 1] # c
        Rx1x2[i,3] = R[ 1- 1,n1- 1] * R[ 2- 1,n2- 1]\
                   + R[ 2- 1,n1- 1] * R[ 1- 1,n2- 1] # d
        Rx1x2[i,4] = R[ 1- 1,n1- 1] * R[ 3- 1,n2- 1]\
                   + R[ 3- 1,n1- 1] * R[ 1- 1,n2- 1] # e
        Rx1x2[i,5] = R[ 2- 1,n1- 1] * R[ 3- 1,n2- 1]\
                   + R[ 3- 1,n1- 1] * R[ 2- 1,n2- 1] # f
    # 5*6, dxy~dz2, x2~yz
    '''           [     "x2",     "y2",     "z2",     "xy",     "xz",    "yz"] '''
    D = np.array([[        0,        0,        0,        2,        0,       0], # dxy
                  [        0,        0,        0,        0,        0,       2], # dyz
                  [        0,        0,        0,        0,        2,       0], # dzx
                  [        1,       -1,        0,        0,        0,       0], # dx2-y2
                  [   -1/sr3,   -1/sr3,      sr3,        0,        0,       0], # dz2
                  ])
    DR = D @ Rx1x2 # 5*6, dxy~dz2, a~f
    # 5*6, dxy~dz2, a~f
    '''            [    "a",    "b",    "c",    "d",    "e",    "f"] '''
    CD = np.array([[      0,      0,      0,    1/2,      0,      0], # dxy
                   [      0,      0,      0,      0,      0,    1/2], # dyz
                   [      0,      0,      0,      0,    1/2,      0], # dzx
                   [      1,      0,    1/3,      0,      0,      0], # dx2-y2
                   [      0,      0,  sr3/3,      0,      0,      0]  # dz2
                   ]) 
    RD = DR @ CD.T
    
    return RD.T
    