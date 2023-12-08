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

keywords=["Name", "Dim", "Spin", "Nbr",
          "LatType", "LatVec", "AtTpNum","Bases",
          "AtomSite", "supercellSize","supercellVacancy","supercellSubstitution",
          "supercellInterstitial","supercellNbr"]
def ReadInput(fileName):
    """

    :param fileName: configuration containing the info of the lattice
    :return: info of lattice
    """
    contents = readContents(fileName)
    ParaIn = {}
    ParaIn["supercell"]={}
    for item in contents:
        kw = item[0]
        if kw == "Name":
            Name = item[1]
            ParaIn["Name"]= Name if Name else "Material"

        elif kw == "Dim":
            Dim = int(item[1])
            ParaIn["Dimension"]=Dim  if Dim  else 3
        elif kw == "Spin":
            Spn = int(item[1])
            ParaIn["Spin"]=Spn
        elif kw == "Nbr":
            n = int(item[1])

        elif kw == "LatType":
            LatType = item[1]
            ParaIn["Lattice type"]=LatType
        elif kw == "LatVec":
            #each row of Lv is a basis vector
            Lv = []
            for vecStr in item[1:]:
                Lv.append(str2num(vecStr, "float"))
            Lv = np.array(Lv)
            ParaIn["LatticeVector"]=Lv
        elif kw=="LatVecSG":
            LvSG=[]
            for vecStr in item[1:]:
                LvSG.append(str2num(vecStr,"float"))
            LvSG=np.array(LvSG)
            ParaIn["SpaceGroupLatticeVector"]=LvSG
        elif kw == "AtTpNum":
            NumAtType = int(item[1])
        elif kw == "Bases":
            AtName = []
            AtNum = []
            AtOrb = []
            AtOrbStr=[]
            for strList in item[1:]:
                AtName.append(strList[0])
                AtNum.append(int(strList[1]))
                AtOrb.append(GetAtOrb(strList[2:],OrbGrp,OrbIdv))
                AtOrbStr.append(strList[2:])
            ParaIn["AtomName"]=AtName
            ParaIn["AtomNumber"]=AtNum
            ParaIn["AtomOrbital"]=np.array(AtOrb)
            ParaIn["AtomOrbitalStr"]=AtOrbStr
        elif kw == "AtomSite":
            AtSite = []
            for strVec in item[1:]:
                AtSite.append(str2num(strVec, "float"))
            AtSite = np.array(AtSite)
            ParaIn["AtomSite"]=AtSite
        elif kw == "supercellSize":
            supercellSize = str2num(item[1], "int")
            if len(supercellSize)!=3:
                raise ValueError("Please give supercell size in 3 directions.")
            ParaIn["supercell"]["supercellSize"]=supercellSize
        elif kw == "supercellVacancy":
            supercellVacList = []
            for row in item[1:]:
                tmpList = [str2num(row[0:3], "int"), int(row[3]), row[4]]
                supercellVacList.append(tmpList)
            ParaIn["supercell"]["supercellVacancy"]=supercellVacList
        elif kw == "supercellSubstitution":
            supercellSubsList = []
            for row in item[1:]:
                tmpList = [str2num(row[0:3], "int"), int(row[3]), row[4], row[5:]]
                supercellSubsList.append(tmpList)
            ParaIn["supercell"]["supercellSubstitution"]=supercellSubsList
        elif kw == "supercellInterstitial":
            supercellInterstitialList = []
            for row in item[1:]:
                tmpList = [str2num(row[0:3], "int"), str2num(row[3:6], "float"), row[6], row[7:]]
                supercellInterstitialList.append(tmpList)
            ParaIn["supercell"]["supercellInterstitial"]=supercellInterstitialList
        elif kw=="supercellNbr":
            ParaIn["supercell"]["supercellNbr"]=int(item[1][0])
    inConfigFolder = str(pathlib.Path(fileName).parent)
    if len(ParaIn["supercell"])>0:
        ParaIn["supercell"]["baseLatticeType"]=LatType
    ParaIn["Folder"]=inConfigFolder
    Nbr = [n, n, 0 if Dim == 2 else n]
    ParaIn["NeighborNumber"]=Nbr
    # Relate atome index iAt with AtomType
    AtTypeInd = []
    for iAt in range(NumAtType):
        for jAt in range(AtNum[iAt]):
            AtTypeInd.append(iAt)
    AtTypeInd = np.array(AtTypeInd, "int")
    ParaIn["AtomTypeIndex"]=AtTypeInd
    ParaIn["OrbIdv"]=OrbIdv


    return ParaIn


def str2num(numStrList, dtype="float"):
    """

    :param numStrList: number string to be converted
    :param dtype: converted data type
    :return: numerical values of the string
    """
    if dtype == "float":
        ret = np.array([float(elem) for elem in numStrList])
        return ret
    if dtype == "int":
        ret = np.array([int(elem) for elem in numStrList])
        return ret
    if dtype == "complex":
        ret = np.array([complex(elem) for elem in numStrList])
        return ret

def readContents(fileName):
    """

    :param fileName: configuration containing the info of the lattice and supercell
    :return: contents in file
    """
    f = [line.strip() for line in open(fileName, "r").readlines()]
    contents = []
    for oneline in f:
        lineContents = oneline.split()
        if len(lineContents) == 0:  # skip empty lines
            continue
        if lineContents[0] in keywords:  # the variable keywords is defined in class lattice
            contents.append(lineContents)
        else:
            if len(contents) == 0:
                raise ValueError("file not starting with a keyword.")
            contents[-1].append(lineContents)
    return contents


def GetAtOrb(Orb,OrbGrpVal,OrbIdvVal):
    # n=1~7 or unknown
    Orb0 = np.zeros(94,int) #1+4+9+16+16+16+16+16
    for Orbi in Orb:
        IndGrp  = np.where(OrbGrpVal ==Orbi)[0]
        IndIdv  = np.where(OrbIdvVal ==Orbi)[0]
        if len(IndGrp):
            Orb0[IndGrp[0]] = 1
        elif len(IndIdv):
            Orb0[IndIdv[0]] = 1
        else:
            print("No orbital named" + Orbi + "!")
    return Orb0


#######################################
#OOP

# class readInput:
#     #mixin class, module for reading config file
#
#     def readContents(self,fileName):
#         """
#
#         :param fileName: configuration containing the info of the lattice and supercell
#         :return: contents in file
#         """
#         f = [line.strip() for line in open(fileName, "r").readlines()]
#         contents = []
#         for oneline in f:
#             lineContents = oneline.split()
#             if len(lineContents) == 0:  # skip empty lines
#                 continue
#             if lineContents[0] in self.keywords:#the variable keywords is defined in class lattice
#                 contents.append(lineContents)
#             else:
#                 if len(contents) == 0:
#                     raise ValueError("file not starting with a keyword.")
#                 contents[-1].append(lineContents)
#
#         return contents
#
#     @staticmethod
#     def str2num(numStrList, dtype="float"):
#         """
#
#         :param numStrList: number string to be converted
#         :param dtype: converted data type
#         :return:
#         """
#
#         if dtype == "float":
#             ret = np.array([float(elem) for elem in numStrList])
#             return ret
#         if dtype == "int":
#             ret = np.array([int(elem) for elem in numStrList])
#             return ret
#         if dtype == "complex":
#             ret = np.array([complex(elem) for elem in numStrList])
#             return ret
#
#     def readInput(self,fileName):
#         """
#
#         :param fileName:
#         :return: configuration containing the info of the lattice and supercell
#         """
#         contents = self.readContents(fileName)
#         for item in contents:
#             kw = item[0]
#             if kw == "Name":
#                 self.Name = item[1]
#             elif kw == "Dim":
#                 self.Dim = int(item[1])
#             elif kw == "Spin":
#                 self.Spn = int(item[1])
#             elif kw == "Nbr":
#                 n = int(item[1])
#                 self.Nbr = [n, n, 0 if self.Dim == 2 else n]
#             elif kw == "LatType":
#                 self.LatType = item[1]
#             elif kw == "LatVec":
#                 Lv = []
#                 for vecStr in item[1:]:
#                     Lv.append(self.str2num(vecStr, "float"))
#                 self.Lv = np.array(Lv)
#             elif kw == "AtTpNum":
#                 self.NumAtType = int(item[1])
#             elif kw == "Bases":
#                 self.AtName = []
#                 self.AtNum = []
#                 self.AtOrb = []
#                 for strList in item[1:]:
#                     self.AtName.append(strList[0])
#                     self.AtNum.append(int(strList[1]))
#                     self.AtOrb.append(self.GetAtOrb(strList[2:]))
#             elif kw == "AtomSite":
#                 AtSite = []
#                 for strVec in item[1:]:
#                     AtSite.append(self.str2num(strVec, "float"))
#                 self.AtSite = np.array(AtSite)
#             elif kw == "supercellSize":
#                 self.supercellSize = self.str2num(item[1], "int")
#             elif kw=="supercellVacancy":
#                 self.supercellVacList=[]
#                 for row in item[1:]:
#                     tmpList=[self.str2num(row[0:3],"int"),int(row[3]),row[4]]
#                     self.supercellVacList.append(tmpList)
#             elif kw=="supercellSubstitution":
#                 self.supercellSubsList=[]
#                 for row in item[1:]:
#                     tmpList=[self.str2num(row[0:3],"int"), int(row[3]),row[4], row[5:]]
#                     self.supercellSubsList.append(tmpList)
#             elif kw=="supercellInterstitial":
#                 self.supercellInterstitialList=[]
#                 for row in item[1:]:
#                     tmpList=[self.str2num(row[0:3],"int"),self.str2num(row[3:6],"float"),row[6],row[7:]]
#                     self.supercellInterstitialList.append(tmpList)
#
#
#
#     def GetAtOrb(self,Orb):
#         """
#
#         :param Orb: a list of strings of orbitals
#         :return: vector of 1 and 0 to represent orbitals
#         """
#         # n=1~7 or unknown
#         Orb0 = np.zeros(94, int)  # 1+4+9+16+16+16+16+16
#         for Orbi in Orb:
#             IndGrp = np.where(self.OrbGrp == Orbi)[0]
#             IndIdv = np.where(self.OrbIdv == Orbi)[0]
#             if len(IndGrp):
#                 Orb0[IndGrp[0]] = 1
#             elif len(IndIdv):
#                 Orb0[IndIdv[0]] = 1
#             else:
#                 print("No orbital named" + Orbi + "!")
#         return Orb0
#
