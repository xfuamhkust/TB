import numpy as np
import sympy as smp
import copy


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
            tij = "t_{%d"%(iClas) + "," + "%d}"%(FreIndii)
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
            Ati = AtName[AtTypeInd[iAt]]; Atj = AtName[AtTypeInd[jAt]]
            Orbi = OrbIdv[iOrb];            Orbj = OrbIdv[jOrb]
            Stti = Ati + "," + Orbi + "," + Lati
            Sttj = Atj + "," + Orbj + "," + Latj
            Hopij = "<" + Stti + "|H|" + Sttj + ">"
            tij = "t_{%d"%(iClas) + "," + "%d}"%(iHop)
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
            tUnf = "t_{%d"%(iClas) + ",%d}"%(UnfIndi[iUnf])
            f.write(tUnf + " = ")
            for iFre, Val in UnfVali[iUnf][:-1]:
                tFre = "t_{%d"%(iClas) + ",%d}"%(iFre)
                f.write("%f"%Val + " * " + tFre + " + ")
            if len(UnfVali[iUnf]):
                iFre, Val = UnfVali[iUnf][-1]
                tFre = "t_{%d"%(iClas) + ",%d}"%(iFre)
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


def writeHkAnalytic(ParaIn,ParaRel,ParaHmtR):
    """

    :param ParaIn:
    :param ParaRel:
    :param ParaHmtR
    :return: print and return Hk's analytic form
    """

    # Parameters
    Name = ParaIn["Name"]
    OrbIdv = ParaIn["OrbIdv"]
    AtName = ParaIn["AtomName"]
    AtNum = ParaIn["AtomNumber"]
    AtTypeInd = ParaIn["AtomTypeIndex"]
    folderName = ParaIn["Folder"]
    LvAtAtOO = ParaRel["LvAtAtOO"]
    HopRelAltClas = ParaRel["HopRelAltClas"]
    NumAt = len(AtTypeInd)
    NumClas = len(HopRelAltClas)
    Lv=ParaIn["LatticeVector"]




    LvAtAtOOSp=[]#holding matrix element symbols
    for elem in LvAtAtOO:
        LvAtAtOOSp.append(elem.tolist())

    rtnFree=[]
    #free coefficients
    for iClas in range(0,NumClas):
        FreIndi=HopRelAltClas[iClas][0]
        LvAtAtOOSpi=LvAtAtOOSp[iClas]

        for j in range(0,len(LvAtAtOOSpi)):
            row=LvAtAtOOSpi[j]
            row.append(0)# real hopping coefficient
            row.append(0)# empty term
            LvAtAtOOSpi[j]=smp.Matrix(row).T
        for FreIndii in FreIndi:
            tij=smp.symbols("t_{"+str(iClas)+"\,"+str(FreIndii)+"}",cls=smp.Symbol,real=True)
            LvAtAtOOSpi[FreIndii][-2]=tij
            rtnFree.append(tij)



    #compute unfree coefficients from free coefficients
    for iClas in range(0,NumClas):
        UnfIndi = HopRelAltClas[iClas][1]
        UnfVali = HopRelAltClas[iClas][2]
        NumUnf = len(UnfIndi)
        for i in range(0,NumUnf):
            unfIndTmp=UnfIndi[i]
            tUnfTmp=0
            if len(UnfVali[i])>0:
                for elem in UnfVali[i]:
                    freeIndTmp,freeCoefTmp=elem
                    tUnfTmp+=LvAtAtOOSp[iClas][freeIndTmp][-2]*freeCoefTmp
            LvAtAtOOSp[iClas][unfIndTmp][-2]=tUnfTmp

    AtOrbInd = ParaHmtR["AtOrbInd"]
    NumStt = len(AtOrbInd)
    atOrb2Num=dict()
    #map entries in AtOrbInd to ordinal numbers
    for i in range(0,len(AtOrbInd)):
        atm,_=AtOrbInd[i]
        atOrb2Num[atm]=dict()

    for i in range(0,len(AtOrbInd)):
        atm,orb=AtOrbInd[i]
        atOrb2Num[atm][orb]=i

    Hk=smp.zeros(NumStt,NumStt)
    k1,k2,k3=smp.symbols("k1,k2,k3",cls=smp.Symbol,real=True)

    for LvAtAtOOSpi in LvAtAtOOSp:
        for rowMat in LvAtAtOOSpi:
            l1,l2,l3,iAtom,jAtom,aOrb,bOrb,tRe,_=rowMat
            xTmp,yTmp,zTmp=l1*Lv[0,:]+l2*Lv[1,:]+l3*Lv[2,:]
            # print(iAtom,aOrb)
            rowNum=atOrb2Num[iAtom][aOrb]
            colNum=atOrb2Num[jAtom][bOrb]
            Hk[rowNum,colNum]+=smp.exp(smp.I*(k1*xTmp+k2*yTmp+k3*zTmp))*tRe



    #return Hk and symbols
    return Hk,k1,k2,k3,rtnFree






    
    