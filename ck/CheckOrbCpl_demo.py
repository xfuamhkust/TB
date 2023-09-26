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


def CheckOrbitalCompleteness(ParaIn,ParaSym):
    
    # Parameters
    AtOrb        = ParaIn["AtomOrbital"]
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