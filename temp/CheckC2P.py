import numpy as np
from Conv2Prim import Conv2Prim

def GenConv():
    BraLat0 = ["aP","mP","mS","oP","oS","oI","oF",\
               "tP","tP","hP","hR","cP","cI","cF"]
    Dim = 3
    NumAt = np.random.poisson(3) + 3
    AtC0 = np.random.random((NumAt,3))
    AtTypeInd0 = np.arange(NumAt)
    BraLat = BraLat0[np.random.randint(len(BraLat0))]
    a,b,c = np.random.random(3)
    th0   = np.random.random()*np.pi
    if BraLat == "aP":
        LvC = np.random.random((3,3))
        AtC = AtC0
        AtTypeInd = AtTypeInd0
    elif BraLat == "mP":
        LvC = np.array([[             a, 0,             0],
                        [             0, b,             0],
                        [ c*np.cos(th0), 0, c*np.sin(th0)]])
        AtC = AtC0
        AtTypeInd = AtTypeInd0
    elif BraLat == "mS":
        LvC = np.array([[             a, 0,             0],
                        [             0, b,             0],
                        [ c*np.cos(th0), 0, c*np.sin(th0)]])
        AtC = np.zeros((2*NumAt,3))
        AtC[    0:  NumAt] = AtC0
        AtC[NumAt:2*NumAt] = AtC0 + np.array([[1/2,1/2,  0]])
        AtTypeInd = np.tile(AtTypeInd0,2)
    elif BraLat == "oP":
        LvC = np.array([[ a, 0, 0],
                        [ 0, b, 0],
                        [ 0, 0, c]])
        AtC = AtC0
        AtTypeInd = AtTypeInd0
    elif BraLat == "oS":
        LvC = np.array([[ a, 0, 0],
                        [ 0, b, 0],
                        [ 0, 0, c]])
        AtC = np.zeros((2*NumAt,3))
        AtC[    0:  NumAt] = AtC0
        AtC[NumAt:2*NumAt] = AtC0 + np.array([[1/2,1/2,  0]])
        AtTypeInd = np.tile(AtTypeInd0,2)
    elif BraLat == "oI":
        LvC = np.array([[ a, 0, 0],
                        [ 0, b, 0],
                        [ 0, 0, c]])
        AtC = np.zeros((2*NumAt,3))
        AtC[    0:  NumAt] = AtC0
        AtC[NumAt:2*NumAt] = AtC0 + np.array([[1/2,1/2,1/2]])
        AtTypeInd = np.tile(AtTypeInd0,2)
    elif BraLat == "oF":
        LvC = np.array([[ a, 0, 0],
                        [ 0, b, 0],
                        [ 0, 0, c]])
        AtC = np.zeros((4*NumAt,3))
        AtC[0*NumAt:1*NumAt] = AtC0
        AtC[1*NumAt:2*NumAt] = AtC0 + np.array([[  0,1/2,1/2]])
        AtC[2*NumAt:3*NumAt] = AtC0 + np.array([[1/2,  0,1/2]])
        AtC[3*NumAt:4*NumAt] = AtC0 + np.array([[1/2,1/2,  0]])
        AtTypeInd = np.tile(AtTypeInd0,4)
    elif BraLat == "tP":
        LvC = np.array([[ a, 0, 0],
                        [ 0, a, 0],
                        [ 0, 0, c]])
        AtC = AtC0
        AtTypeInd = AtTypeInd0
    elif BraLat == "tI":
        LvC = np.array([[ a, 0, 0],
                        [ 0, a, 0],
                        [ 0, 0, c]])
        AtC = np.zeros((2*NumAt,3))
        AtC[    0:  NumAt] = AtC0
        AtC[NumAt:2*NumAt] = AtC0 + np.array([[1/2,1/2,1/2]])
        AtTypeInd = np.tile(AtTypeInd0,2)
    elif BraLat == "hP":
        LvC = np.array([[   a,              0, 0],
                        [-a/2, a*np.sqrt(3)/2, 0],
                        [   0,              0, c]])
        AtC = AtC0
        AtTypeInd = AtTypeInd0
    elif BraLat == "hR":
        LvC = np.array([[ a,                0, c],
                        [-a/2, a*np.sqrt(3)/2, c],
                        [-a/2,-a*np.sqrt(3)/2, c]])
        AtC = AtC0
        AtTypeInd = AtTypeInd0
    elif BraLat == "cP":
        LvC = np.array([[1,0,0],
                        [0,1,0],
                        [0,0,1]])*a
        AtC = AtC0
        AtTypeInd = AtTypeInd0
    elif BraLat == "cI":
        LvC = np.array([[1,0,0],
                        [0,1,0],
                        [0,0,1]])*a
        AtC = np.zeros((2*NumAt,3))
        AtC[    0:  NumAt] = AtC0
        AtC[NumAt:2*NumAt] = AtC0 + np.array([[1/2,1/2,1/2]])
        AtTypeInd = np.tile(AtTypeInd0,2)
    elif BraLat == "cF":
        LvC = np.array([[1,0,0],
                        [0,1,0],
                        [0,0,1]])*a
        AtC = np.zeros((4*NumAt,3))
        AtC[0*NumAt:1*NumAt] = AtC0
        AtC[1*NumAt:2*NumAt] = AtC0 + np.array([[  0,1/2,1/2]])
        AtC[2*NumAt:3*NumAt] = AtC0 + np.array([[1/2,  0,1/2]])
        AtC[3*NumAt:4*NumAt] = AtC0 + np.array([[1/2,1/2,  0]])
        AtTypeInd = np.tile(AtTypeInd0,4)
    r       = np.random.random(3)
    x, y, z = r/np.linalg.norm(r)
    theta   = np.random.random()*2*np.pi
    ct = np.cos(theta); st = np.sin(theta)
    R = np.array([[ct+(1-ct)*x**2,  (1-ct)*x*y-st*z, (1-ct)*x*z+st*y],
                  [(1-ct)*y*x+st*z, ct+(1-ct)*y**2,  (1-ct)*y*z-st*x],
                  [(1-ct)*z*x-st*y, (1-ct)*z*y+st*x, ct+(1-ct)*z**2 ]])
    LvC = LvC @ R.T
    
    return BraLat, LvC, AtC, AtTypeInd, Dim
        
# def CheckC2P(n=10):
#     for i in range(n):
#         BraLatIn, LvC, AtC, AtTypeInd, Dim = GenConv()
#         print(BraLatIn, LvC, AtC, AtTypeInd, Dim)
#         LatSys, BraLatOut, LvP = Conv2Prim(LvC, AtC, AtTypeInd, Dim)
#         print(BraLatIn,BraLatOut)
#         if BraLatIn != BraLatOut:
#             print("Wrong!")
#             break
  
# CheckC2P(1000)      

for i in range(1000):
    BraLatIn, LvC, AtC, AtTypeInd, Dim = GenConv()
    print(BraLatIn, LvC, AtC, AtTypeInd, Dim)
    LatSys, BraLatOut, LvP = Conv2Prim(LvC, AtC, AtTypeInd, Dim, error=1e-6)
    print(BraLatIn,BraLatOut)
    if BraLatIn != BraLatOut:
        print("Wrong!")
        break

# LatSys, BraLatOut, LvP = Conv2Prim(LvC, AtC, AtTypeInd, Dim)




