import numpy as np


def WriteRelation(ParaIn,ParaRel):
    
    # Parameters
    # Name          = ParaIn["Name"]
    OrbIdv        = ParaIn["OrbIdv"]
    LvAtAtOO      = ParaRel["LvAtAtOO"]
    HopRelAltClas = ParaRel["HopRelAltClas"]
    folderName=ParaIn["Folder"]
    NumClas = len(HopRelAltClas)
    
    # Write hopping relations
    NameHopRel = folderName + "/HopRel.txt"
    f = open(NameHopRel,"w")
    f.write("###HoppingTerms###\n\n")
    for iClas in range(NumClas):
        LvAtAtOOi = LvAtAtOO[iClas]
        NumHopi = len(LvAtAtOOi)
        for iHop in range(NumHopi):
            for j in LvAtAtOOi[iHop]:
                f.write("%d"%j + " ")
            f.write("\n")
        f.write("\n")
    f.write("###HoppingRelations###\n\n")
    for iClas in range(NumClas):
        FreIndi = HopRelAltClas[iClas][0]
        UnfIndi = HopRelAltClas[iClas][1]
        UnfVali = HopRelAltClas[iClas][2]
        f.write("[ ")
        for i in FreIndi:
            f.write("%d"%i + " ")
        f.write("]\n")
        f.write("[ ")
        for i in UnfIndi:
            f.write("%d"%i + " ")
        f.write("]\n")
        # print(UnfVali)
        for UnfValii in UnfVali:
            for UnfValiii in UnfValii:
                i,j = UnfValiii
                f.write("%d"%i + " " + "%f"%j + " ")
            f.write("\n")
        f.write("\n")
    f.close()
    
    # Write input instrunction of free hopping terms
    NameHopFree = folderName+ "/HopFree.txt"
    f = open(NameHopFree,"w")
    for iClas in range(NumClas):
        FreIndi = HopRelAltClas[iClas][0]
        LvAtAtOOi = LvAtAtOO[iClas]
        for FreIndii in FreIndi:
            n1, n2, n3, iAt, jAt, iOrb, jOrb = LvAtAtOOi[FreIndii]
            iOrbName = OrbIdv[iOrb]; jOrbName = OrbIdv[jOrb]
            f.write("%d"%n1 + " " + "%d"%n2 + " " + "%d"%n3 + " " + "%d"%iAt +\
                    " " + "%d"%jAt + " " + iOrbName + " " + jOrbName + "\n")
    f.close()