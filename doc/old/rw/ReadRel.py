import numpy as np


def ReadRelation(FileName):
    
    f = [line.strip() for line in open(FileName,"r").readlines()]
    for i in range(len(f)):
        if f[i] == "###HoppingTerms###":
            IndHop = i
        elif f[i] == "###HoppingRelations###":
            IndHopRel = i
    
    # Read hopping terms
    LvAtAtOO = []
    LvAtAtOOi = []
    for i in range(IndHop+2,IndHopRel):
        if len(f[i]):
            fsi = f[i].split()
            LvAtAtOOi.append([int(fsi[j]) for j in range(len(fsi))])
            
        else:
            LvAtAtOO.append(np.array(LvAtAtOOi))
            LvAtAtOOi = []
    
    # Read hopping relations
    HopRelAltClas = []
    IndClas = []
    for i in range(IndHopRel+2,len(f)):
        if len(f[i]) and f[i][0] == "[":
            IndClas.append(i)
    IndClas = IndClas[::2]
    for IndClasi in IndClas:
        # Indices of free terms
        fs = f[IndClasi  ][1:-1].split()
        FreIndi = [int(fs[j]) for j in range(len(fs))]
        # Indices of unfree terms
        fs = f[IndClasi+1][1:-1].split()
        UnfIndi = [int(fs[j]) for j in range(len(fs))]
        # Values of unfree terms.
        UnfVali = []
        for i in range(len(UnfIndi)):
            fs = f[IndClasi+i+2].split()
            UnfVali.append([])
            for j in range(len(fs)//2):
                UnfVali[-1].append([int(fs[2*j]),float(fs[2*j+1])])
        HopRelAltClas.append([FreIndi,UnfIndi,UnfVali])
    
    # Output
    ParaRel = {"LvAtAtOO":      LvAtAtOO,
               "HopRelAltClas": HopRelAltClas,
               }
    
    return ParaRel