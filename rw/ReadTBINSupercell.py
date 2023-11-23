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

keywords=["Name", "Dim", "Spin", "Nbr",
          "LatType", "LatVec", "AtTpNum","Bases",
          "AtomSite", "supercellSize","supercellVacancy","supercellSubstitution",
          "supercellInterstitial"]


def readContents(fileName):
    """

    :param fileName: configuration containing the info of the lattice and supercell
    :return: contents in file
    """
    f = [line.strip() for line in open(fileName, "r").readlines()]
    contents=[]
    for oneline in f:
        lineContents=oneline.split()
        if len(lineContents)==0:#skip empty lines
            continue
        if lineContents[0] in keywords:
            contents.append(lineContents)
        else:
            if len(contents)==0:
                raise ValueError("file not starting with a keyword.")
            contents[-1].append(lineContents)

    return contents

def str2num(numStrList,dtype="float"):
    """

    :param strVec: number string to be converted
    :param dtype: converted data type
    :return:
    """

    if dtype=="float":
        ret=np.array([float(elem) for elem in numStrList])
        return ret
    if dtype=="int":
        ret=np.array([int(elem) for elem in numStrList])
        return ret
    if dtype=="complex":
        ret=np.array([complex(elem) for elem in numStrList])
        return ret


def readInputSupercell(fileName):
    """

    :param fileName:  configuration containing the info of the lattice and supercell
    :return: info of supercell
    """
    contents=readContents(fileName)
    for item in   contents:
        kw=item[0]
        if kw=="Name":
            Name=item[1]
        elif kw=="Dim":
            Dim=int(item[1])
        elif kw=="Spin":
            Spn=int(item[1])
        elif kw=="Nbr":
            n=int(item[1])
            Nbr=[n,n,0 if Dim==2 else n]
        elif kw=="LatType":
            LatType=item[1]
        elif kw=="LatVec":
            Lv=[]
            for vecStr in item[1:]:
                Lv.append(str2num(vecStr,"float"))
            Lv=np.array(Lv)
        elif kw=="AtTpNum":
            NumAtType=int(item[1])
        elif kw=="Bases":
            AtName=[]
            AtNum=[]
            AtOrb=[]
            for strList in item[1:]:
                AtName.append(strList[0])
                AtNum.append(int(strList[1]))
                AtOrb.append(GetAtOrb(strList[2:]))
        elif kw=="AtomSite":
            AtSite=[]
            for strVec in item[1:]:
                AtSite.append(str2num(strVec,"float"))
            AtSite=np.array(AtSite)
        elif kw=="supercellSize":
            supercellSize=str2num(item[1],"int")
            supercellVacList=[]
            supercellSubsList=[]
            supercellInterstitialList=[]
        elif kw=="supercellVacancy":












def GetAtOrb(Orb):
    """

    :param Orb: a list of strings of orbitals
    :return: vector of 1 and 0 to represent orbitals
    """
    # n=1~7 or unknown
    Orb0 = np.zeros(94, int)  # 1+4+9+16+16+16+16+16
    for Orbi in Orb:
        IndGrp = np.where(OrbGrp == Orbi)[0]
        IndIdv = np.where(OrbIdv == Orbi)[0]
        if len(IndGrp):
            Orb0[IndGrp[0]] = 1
        elif len(IndIdv):
            Orb0[IndIdv[0]] = 1
        else:
            print("No orbital named" + Orbi + "!")
    return Orb0









fileName="../data/ABO3/supercell_TBIN_ABO3.txt"
readInputSupercell(fileName)
