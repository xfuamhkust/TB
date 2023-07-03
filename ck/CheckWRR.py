import numpy as np

def CheckWriteReadRelation(ParaRelWrite,ParaRelRead, error = 1e-9):
    
    # Parameters
    LvAtAtOOWrite      = ParaRelWrite["LvAtAtOO"]
    LvAtAtOORead       = ParaRelRead["LvAtAtOO"]
    HopRelAltClasWrite = ParaRelWrite["HopRelAltClas"]
    HopRelAltClasRead  = ParaRelRead["HopRelAltClas"]
    
    # Compare LvAtAtOO
    NumClasWrite = len(LvAtAtOOWrite)
    NumClasRead  = len(LvAtAtOORead)
    if NumClasWrite == NumClasRead:
        NumClas = NumClasRead
    else:
        print("Wrong NumClas!")
        return False
    for iClas in range(NumClas):
        if np.sum(abs(LvAtAtOOWrite[iClas]-LvAtAtOORead[iClas])):
            print("Wrong LvAtAtOO! iClas = %d"%iClas)
            return False
    # Compare HopRelAltClas
    NumClasWrite = len(HopRelAltClasWrite)
    NumClasRead  = len(HopRelAltClasRead)
    if NumClasWrite == NumClasRead:
        NumClas = NumClasRead
    else:
        print("Wrong NumClas!")
        return False
    for iClas in range(NumClas):
        FreIndiWrite, UnfIndiWrite, UnfValiWrite = HopRelAltClasWrite[iClas]
        FreIndiRead , UnfIndiRead , UnfValiRead  = HopRelAltClasRead[iClas]
        d1 = np.sum(abs(np.array(FreIndiWrite,int)-np.array(FreIndiRead,int)))
        d2 = np.sum(abs(np.array(UnfIndiWrite,int)-np.array(UnfIndiRead,int)))
        if d1 == 0 and d2 == 0:
            pass
        else:
            print("Wrong HopRelAltClas! iClas = %d"%iClas)
            return False
        for i in range(len(UnfValiWrite)):
            di = np.sum(abs(np.array(UnfValiWrite[i])-np.array(UnfValiRead[i])))
            if di > error: 
                print("Wrong HopRelAltClas! iClas = %d"%iClas)
                return False
    print("True!")
    return True
        