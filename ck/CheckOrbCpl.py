import numpy as np


def CheckOrbitalCompleteness(ParaIn,ParaSym):
    
    # Parameters
    AtOrb        = ParaIn["AtomOrbital"]
    OrbIdv       = ParaIn["OrbIdv"]
    AtName       = ParaIn["AtomName"]
    AtNum        = ParaIn["AtomNumber"]
    SymOrb       = ParaSym["SymOrb"]
    NumAtType, NumOrb = AtOrb.shape
    NumSym = len(SymOrb[0])
    error = 1e-6
    
    # Write Symmetries acting on SPDF together
    IndSPDF = np.array([1,1,3,1,3,5,1,3,5,7,1,3,5,7,1,3,5,7,1,3,5,7,1,3,5,7])
    SymSPDF = np.zeros((NumSym,94,94))
    count = 0
    for Indi in IndSPDF:
        SymSPDF[:,count:count+Indi,count:count+Indi] = SymOrb[(Indi-1)//2]
        count += Indi
    
    # Classify symmetry-related orbitals by off-diagonal
    IndNonZero = np.sum(abs(SymSPDF),axis=0) > error # Nonzero off-diagnoal element in SymSPDF
    AtOrbNew = np.copy(AtOrb)
    for iAt in range(NumAtType):
        IndOne = np.where(AtOrb[iAt])[0] # Input orbital
        for IndOnei in IndOne:
            AtOrbNew[iAt,np.where(IndNonZero[IndOnei])[0]] = 1
    
    # Write Symmetries acting on input orbitals
    Ind  = [np.ix_(np.where(AtOrbNew[iAt])[0],np.where(AtOrbNew[iAt])[0]) for iAt in range(NumAtType)]
    SymAtOrb = []
    for iAt in range(NumAtType):
        SymAtOrbi = []
        for iSym in range(NumSym):
            SymAtOrbi.append(SymSPDF[iSym][Ind[iAt]])
        for i in range(AtNum[iAt]):
            SymAtOrb.append(np.array(SymAtOrbi))
    
    # Inform the user whether input orbital is complete
    dAtOrb = AtOrbNew - AtOrb
    if np.sum(dAtOrb) < error:
        print("Input orbital is complete under Symmetry.")
    else:
        print("Input orbital is incomplete under Symmetry.")
        print("The following orbitals are automatically added:")
        for iAt in range(NumAtType):
            IndOne = np.where(dAtOrb[iAt])[0]
            String = "Atom " + AtName[iAt] + ": "
            if len(IndOne):
                for IndOnei in IndOne[:-1]:
                    String += OrbIdv[IndOnei] + ", "
                String += OrbIdv[IndOne[-1]]
                print(String)
    
    # OutPut
    ParaIn["AtomOrbital"] = AtOrbNew
    ParaSym["SymSPDF"] = SymSPDF
    ParaSym["SymAtOrb"] = SymAtOrb
    
    return ParaIn,ParaSym