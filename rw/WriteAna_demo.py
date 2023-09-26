import numpy as np

def WriteHamiltonianAnalytic(ParaIn,ParaRel):
    
    # Parameters
    Name          = ParaIn["Name"]
    OrbIdv        = ParaIn["OrbIdv"]
    AtName        = ParaIn["AtomName"]
    AtNum         = ParaIn["AtomNumber"]
    AtTypeInd     = ParaIn["AtomTypeIndex"]
    folderName    = ParaIn["Folder"]
    LvAtAtOO      = ParaRel["LvAtAtOO"]
    HopRelAltClas = ParaRel["HopRelAltClas"]
    NumAt = len(AtTypeInd)
    NumClas = len(HopRelAltClas)
    
    # Write Analytic of real space Hamiltonian
    NameHrAna = folderName + "/HrAna.txt"
    f = open(NameHrAna,"w")
    # Step1: write free hopping terms t_{i,j}
    f.write("Free hopping terms\n\n")
    for iClas in range(NumClas):
        FreIndi = HopRelAltClas[iClas][0]
        LvAtAtOOi = LvAtAtOO[iClas]
        for FreIndii in FreIndi:
            tij = "t_{%d"%(iClas+1) + "," + "%d}"%(FreIndii+1)
            f.write(tij + "\n")
    f.write("\n")
    # Step2: define hopping terms t_{i,j} = <Stti|H|Sttj>
    f.write("Definition of hopping terms\n\n")
    for iClas in range(NumClas):
        LvAtAtOOi = LvAtAtOO[iClas]
        NumHopi = len(LvAtAtOOi)
        for iHop in range(NumHopi):
            n1, n2, n3, iAt, jAt, iOrb, jOrb = LvAtAtOOi[iHop]
            Lati = "(0,0,0)";               Latj = "(%d"%n1 + ",%d"%n2 + ",%d"%n3 + ")"
            Ati = AtName[AtTypeInd[iAt-1]]; Atj = AtName[AtTypeInd[jAt-1]]
            Orbi = OrbIdv[iOrb];            Orbj = OrbIdv[jOrb]
            Stti = Ati + "," + Orbi + "," + Lati
            Sttj = Atj + "," + Orbj + "," + Latj
            Hopij = "<" + Stti + "|H|" + Sttj + ">"
            tij = "t_{%d"%(iClas+1) + "," + "%d}"%(iHop+1)
            f.write(tij + " = " + Hopij)
            f.write("\n")
        f.write("\n")
    # Step3: write relations between free ones and unfree ones
    f.write("Relations between hopping terms\n\n")
    for iClas in range(NumClas):
        UnfIndi = HopRelAltClas[iClas][1]
        UnfVali = HopRelAltClas[iClas][2]
        NumUnf  = len(UnfIndi)
        for iUnf in range(NumUnf):
            tUnf = "t_{%d"%(iClas+1) + ",%d}"%(UnfIndi[iUnf]+1)
            f.write(tUnf + " = ")
            for iFre, Val in UnfVali[iUnf][:-1]:
                tFre = "t_{%d"%(iClas+1) + ",%d}"%(iFre+1)
                f.write("%f"%Val + " * " + tFre + " + ")
            if len(UnfVali[iUnf]):
                iFre, Val = UnfVali[iUnf][-1]
                tFre = "t_{%d"%(iClas+1) + ",%d}"%(iFre+1)
                f.write("%f"%Val + " * " + tFre)
            f.write("\n")
        f.write("\n")
    f.close()
    
    '''
    # Write Analytic of k-space Hamiltonian
    NameHkAna = folderName + "/HkAna.txt"
    f = open(NameHkAna,"w")
    # Step1: write free hopping terms t_{i,j}
    f.write("Free hopping terms\n\n")
    for iClas in range(NumClas):
        FreIndi = HopRelAltClas[iClas][0]
        LvAtAtOOi = LvAtAtOO[iClas]
        for FreIndii in FreIndi:
            tij = "t_{%d"%(iClas+1) + "," + "%d}"%(FreIndii+1)
            f.write(tij + "\n")
    f.write("\n")
    # Step2: write Bloch states
    下面从GetHamiltonianReal中copy了部分内容，是为了保证Hk的bases和算能带的bases一致。
    # ----------------------- From GetHamiltonianReal -----------------------#
    LvAtAtOOAll = []
    for iClas in range(NumClas):
        LvAtAtOOAll += LvAtAtOO[iClas].tolist()
    LvAtAtOOAll = np.array(LvAtAtOOAll)
    AtOrbInd = []
    for iAt in range(1,NumAt+1):
        Orbi = np.unique(LvAtAtOOAll[np.where(LvAtAtOOAll[:,3]==iAt)[0],5])
        for Orbii in Orbi:
            AtOrbInd.append([iAt,Orbii])
    AtOrbInd = np.array(AtOrbInd)
    NumStt = len(AtOrbInd)
    # --------------------------------- End ---------------------------------#
    f.write("Bloch states\n\n")
    for iStt in range(NumStt):
        iAt, iOrb = AtOrbInd[iStt]
        Ati = AtName[AtTypeInd[iAt-1]]
        Orbi = OrbIdv[iOrb]
        Stti = "|" + Ati + "," + Orbi + ">"
        f.write(Stti + "\n")
    f.write("\n")
    # Step3: write expressions of k-space Hamiltonian elements by free hopping terms
    f.write("k-space Hamiltonian\n\n")
    for iClas in range(NumClas):
        LvAtAtOOi = LvAtAtOO[iClas]
        NumHopi = len(LvAtAtOOi)
        for iHop in range(NumHopi):
            Ln1, n2, n3, iAt, jAt, iOrb, jOrb = LvAtAtOOi[iHop]
            Lvii       = LvAtAtOOii[0:3]
            AtAtOOii   = LvAtAtOOii[3:7]
            iStt = FindIndex(AtAtOOii,X)
            
            
def FindIndex(x,X,tol = 1e-3):
    dX = np.sum(abs(X - x), axis=1)
    min_dX = min(dX)
    if min_dX < tol:
        ind = np.where(dX == min_dX)[0][0]
        return ind
    else:
        return -1
    '''
    
    
    