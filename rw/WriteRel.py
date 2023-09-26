import numpy as np
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

def WriteRelation(ParaIn,ParaRel):
    
    # Parameters
    Name          = ParaIn["Name"]
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