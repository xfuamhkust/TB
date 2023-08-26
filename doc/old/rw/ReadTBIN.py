import os
import numpy as np

def ReadInput(FileName="data//TBIN.txt"):
    f = [line.strip() for line in open(FileName,"r").readlines()]
    Name=CaptInfo("Name",f);Dim=int(CaptInfo("Dim",f));Spn=int(CaptInfo("Spin",f));
    SGN=int(CaptInfo("SGN",f));NumAtType=int(CaptInfo("AtTpNum",f))
    Nbr0=CaptInfo("Nbr",f)
    if type(Nbr0) == str:
        n = int(Nbr0); Nbr = [n, n, 0 if Dim == 2 else n]
    else:
        Nbr = [int(Nbr0[0]), int(Nbr0[1]), 0 if Dim == 2 else int(Nbr0[2])]
    LvSG   = np.array([CaptInfo("LatVecSG",f,3)[i].split() for i in range(3)],float)
    Lv     = np.array([CaptInfo("LatVec"  ,f,3)[i].split() for i in range(3)],float)
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
    if (1 - os.path.exists("data/" + Name)):
        os.mkdir("data/" + Name)
    ParaIn = {"Name":                       Name if Name else "Material",
              "Dimension":                  Dim  if Dim  else 3,
              "Spin":                       Spn  if Spn  else 0,
              "SpaceGroupNumber":           SGN  if SGN  else 1,
              "NeighborNumber":             Nbr  if Nbr  else 1,
              "LatticeVector":              Lv,
              "SpaceGroupLatticeVector":    LvSG,
              "AtomName":                   AtName,
              "AtomNumber":                 AtNum,
              "AtomSite":                   AtSite,
              "AtomOrbital":                AtOrb,
              "AtomTypeIndex":              AtTypeInd,
              }
    return ParaIn

def CaptInfo(Name, FileList, NumRow = 0):
    NumList = len(FileList)
    for iList in range(NumList):
        FileLine = FileList[iList].split()
        if FileLine[0] == Name:
            if NumRow:
                return FileList[iList+1:iList+1+NumRow]
            else:
                return FileLine[1:] if len(FileLine)-2 else FileLine[1]
            
def GetAtOrb(Orb):
    Orb0 = np.zeros(16,int) #1+3+5+7
    OrbGrp = np.array(["S","P","D","F"])
    OrbIdv = np.array(["s","px","py","pz","dxy","dyz","dzx","dx2-y2","dz2",
              "fz3","fxz2","fyz2","fxyz","fz(x2-y2)","fx(x2-3y2)","fy(3x2-y2)"])
    for Orbi in Orb:
        IndGrp = np.where(OrbGrp==Orbi)[0]
        IndIdv = np.where(OrbIdv==Orbi)[0]
        if len(IndGrp):
            Orb0[IndGrp[0]**2:(IndGrp[0]+1)**2] = 1
        elif len(IndIdv):
            Orb0[IndIdv[0]] = 1
        else:
            print("No orbital named" + Orbi + "!")
    return Orb0