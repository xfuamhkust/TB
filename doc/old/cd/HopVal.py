import numpy as np

def GetHoppingValue(HopValIn, ParaRel):
    
    # Parameters
    LvAtAtOO      = ParaRel["LvAtAtOO"]
    HopRelAltClas = ParaRel["HopRelAltClas"]
    NumClas = len(LvAtAtOO)
    
    # Relate input values to free hopping terms
    HopValInClas = []
    count = 0
    for iClas in range(NumClas):
        FreIndi = HopRelAltClas[iClas][0]
        HopValInClas.append([])
        for i in range(len(FreIndi)):
            if count < len(HopValIn):
                HopValInClas[-1].append(HopValIn[count])
            else:
                HopValInClas[-1].append(0)
            count += 1
    
    # Calculate value of free and unfree hopping terms
    HopValClas = []
    for iClas in range(NumClas):
        FreIndi = HopRelAltClas[iClas][0]
        UnfIndi = HopRelAltClas[iClas][1]
        UnfVali = HopRelAltClas[iClas][2]
        NumHopi = max(FreIndi + UnfIndi) + 1
        HopValInClasi = HopValInClas[iClas]
        HopValClasi = np.zeros(NumHopi)
        for i in range(len(FreIndi)):
            HopValClasi[FreIndi[i]] = HopValInClasi[i]
        for i in range(len(UnfIndi)):
            HopValTemp = 0
            for UnfValiii in UnfVali[i]:
                HopValTemp += UnfValiii[1] * HopValInClasi[np.where(np.array(FreIndi)==UnfValiii[0])[0][0]]
            HopValClasi[UnfIndi[i]] = HopValTemp
        HopValClas.append(HopValClasi)
        
    return HopValClas
    