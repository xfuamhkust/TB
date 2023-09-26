import os
import numpy as np
import re
import pathlib
OrbGrp = np.array(["1S","2S","2P","2P","2P",
                   "3S","3P","3P","3P","3D","3D","3D","3D","3D",
                   "4S","4P","4P","4P","4D","4D","4D","4D","4D",
                   "4F","4F","4F","4F","4F","4F","4F",
                   "5S","5P","5P","5P","5D","5D","5D","5D","5D",
                   "5F","5F","5F","5F","5F","5F","5F",
                   "6S","6P","6P","6P","6D","6D","6D","6D","6D",
                   "6F","6F","6F","6F","6F","6F","6F",
                   "7S","7P","7P","7P","7D","7D","7D","7D","7D",
                   "7F","7F","7F","7F","7F","7F","7F",
                   "S","P","P","P","D","D","D","D","D",
                   "F","F","F","F","F","F","F"])
OrbIdv = np.array(["1s","2s","2px","2py","2pz",
                   "3s","3px","3py","3pz","3dxy","3dyz","3dzx","3dx2-y2","3dz2",
                   "4s","4px","4py","4pz","4dxy","4dyz","4dzx","4dx2-y2","4dz2",
                   "4fz3","4fxz2","4fyz2","4fxyz","4fz(x2-y2)","4fx(x2-3y2)","4fy(3x2-y2)",
                   "5s","5px","5py","5pz","5dxy","5dyz","5dzx","5dx2-y2","5dz2",
                   "5fz3","5fxz2","5fyz2","5fxyz","5fz(x2-y2)","5fx(x2-3y2)","5fy(3x2-y2)",
                   "6s","6px","6py","6pz","6dxy","6dyz","6dzx","6dx2-y2","6dz2",
                   "6fz3","6fxz2","6fyz2","6fxyz","6fz(x2-y2)","6fx(x2-3y2)","6fy(3x2-y2)",
                   "7s","7px","7py","7pz","7dxy","7dyz","7dzx","7dx2-y2","7dz2",
                   "7fz3","7fxz2","7fyz2","7fxyz","7fz(x2-y2)","7fx(x2-3y2)","7fy(3x2-y2)",
                   "s","px","py","pz","dxy","dyz","dzx","dx2-y2","dz2",
                   "fz3","fxz2","fyz2","fxyz","fz(x2-y2)","fx(x2-3y2)","fy(3x2-y2)"])


def ReadInput(FileName):
    """

    :param FileName: configuration containing the info of the lattice
    :return: info of lattice
    """
    f = [line.strip() for line in open(FileName,"r").readlines()]#TODO: check readin error
    Name=CaptInfo("Name",f)
    Dim=int(CaptInfo("Dim",f))
    Spn=int(CaptInfo("Spin",f))
    # SGN=int(CaptInfo("SGN",f));
    NumAtType=int(CaptInfo("AtTpNum",f))
    LatType=CaptInfo("LatType",f)
    # print(LatType)
    if  re.search("[^a-zA-Z]+",LatType) or LatType is None:
        raise ValueError("Invalid lattice type name.")
    LatType=LatType.lower()#lattice type
    # print(LatType)
    Nbr0=CaptInfo("Nbr",f)
    if type(Nbr0) == str:
        n = int(Nbr0); Nbr = [n, n, 0 if Dim == 2 else n]
    else:
        Nbr = [int(Nbr0[0]), int(Nbr0[1]), 0 if Dim == 2 else int(Nbr0[2])]
    # LvSG   = np.array([CaptInfo("LatVecSG",f,3)[i].split() for i in range(3)],float)
    # Lv     = np.array([CaptInfo("LatVec"  ,f,3)[i].split() for i in range(3)],float)
    Bas    = [CaptInfo("Bases",f,NumRow=NumAtType)[i].split() for i in range(NumAtType)]
    AtName = [Bas[i][0] for i in range(NumAtType)]
    AtNum  = [int(Bas[i][1]) for i in range(NumAtType)]
    AtSite = np.array([CaptInfo("AtomSite",f,NumRow=sum(AtNum))[i].split() for i in range(sum(AtNum))],float)
    AtOrb  = np.array([GetAtOrb(Bas[i][2:]) for i in range(NumAtType)])
    # Relate atome index iAt with AtomType
    AtTypeInd = []
    for iAt in range(NumAtType):
        for jAt in range(AtNum[iAt]):
            AtTypeInd.append(iAt)
    AtTypeInd = np.array(AtTypeInd,"int")
    # Creat a new dictory for this material
    # if (1 - os.path.exists("data/" + Name)):
    #     os.mkdir("data/" + Name)
    inConfigFolder=str(pathlib.Path(FileName).parent)
    if LatType=="primitive":
        Lv = np.array([CaptInfo("LatVec", f, 3)[i].split() for i in range(3)], float)
        ParaIn = {"Name":                       Name if Name else "Material",
                  "Dimension":                  Dim  if Dim  else 3,
                  "Spin":                       Spn  if Spn  else 0,
                  "Lattice type":               LatType,
                # "SpaceGroupNumber":           SGN  if SGN  else 1,
                "NeighborNumber":             Nbr  if Nbr  else 1,
                "LatticeVector":              Lv,
                # "SpaceGroupLatticeVector":    LvSG,
                "AtomName":                   AtName,
                "AtomNumber":                 AtNum,
                "AtomSite":                   AtSite,
                "AtomOrbital":                AtOrb,
                "AtomTypeIndex":              AtTypeInd,
                "OrbIdv":                     OrbIdv,
                "Folder":                     inConfigFolder
                }
        return ParaIn
    elif LatType=="conventional":
        LvSG = np.array([CaptInfo("LatVecSG", f, 3)[i].split() for i in range(3)], float)
        ParaIn = {"Name": Name if Name else "Material",
                  "Dimension": Dim if Dim else 3,
                  "Spin": Spn if Spn else 0,
                  "Lattice type": LatType,
                  # "SpaceGroupNumber":           SGN  if SGN  else 1,
                  "NeighborNumber": Nbr if Nbr else 1,
                  # "LatticeVector": Lv,
                  "SpaceGroupLatticeVector":    LvSG,
                  "AtomName": AtName,
                  "AtomNumber": AtNum,
                  "AtomSite": AtSite,
                  "AtomOrbital": AtOrb,
                  "AtomTypeIndex": AtTypeInd,
                  "Folder": inConfigFolder
                  }
        return ParaIn
    else:
        raise ValueError("Wrong lattice name.")

def CaptInfo(Name, FileList, NumRow = 0):
    NumList = len(FileList)
    for iList in range(NumList):
        FileLine = FileList[iList].split()#TODO: more efficient ways of iteration
        if FileLine[0] == Name:
            if NumRow:
                return FileList[iList+1:iList+1+NumRow]
            else:
                return FileLine[1:] if len(FileLine)-2 else FileLine[1]
            
def GetAtOrb(Orb):
    # n=1~7 or unknown
    Orb0 = np.zeros(94,int) #1+4+9+16+16+16+16+16
    for Orbi in Orb:
        IndGrp  = np.where(OrbGrp ==Orbi)[0]
        IndIdv  = np.where(OrbIdv ==Orbi)[0]
        if len(IndGrp):
            Orb0[IndGrp[0]] = 1
        elif len(IndIdv):
            Orb0[IndIdv[0]] = 1
        else:
            print("No orbital named" + Orbi + "!")
    return Orb0