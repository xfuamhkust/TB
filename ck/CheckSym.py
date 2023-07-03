import numpy as np

def CheckSymmetry(ParaSym):
    
    # Parameters
    SymLv  = ParaSym["SymLv"]
    SymXyz = ParaSym["SymXyz"]
    SymOrb = ParaSym["SymOrb"]
    SymSpn = ParaSym["SymSpn"]
    SymS, SymP, SymD, SymF = SymOrb
    
    # MultiForm(SymLv)
    MfLv  = MultiFormLv(SymLv)
    MfXyz = MultiForm(SymXyz)
    MfS   = MultiForm(SymS)
    MfP   = MultiForm(SymP)
    MfD   = MultiForm(SymD)
    MfF   = MultiForm(SymF)
    MfSpn = MultiForm(SymSpn)
    # print(MfLv)

    dMfXyz = CheckMultiForm(SymXyz,MfLv)
    dMfS   = CheckMultiForm(SymS,MfLv)
    dMfP   = CheckMultiForm(SymP,MfLv)
    dMfD   = CheckMultiForm(SymD,MfLv)
    dMfF   = CheckMultiForm(SymF,MfLv)
    dMfSpn = CheckMultiForm(SymSpn,MfLv)
    print(np.linalg.norm(dMfXyz))
    print(np.linalg.norm(dMfS))
    print(np.linalg.norm(dMfP))
    print(np.linalg.norm(dMfD))
    print(np.linalg.norm(dMfF))
    print(np.linalg.norm(dMfSpn))
    # print(SymSpn)
    
def CheckMultiForm(Sym,Mf0):
    NumSym = len(Sym)
    dMf = np.zeros((NumSym,NumSym))
    for i in range(NumSym):
        for j in range(NumSym):
            Symk = Sym[i] @ Sym[j]
            k = Mf0[i,j] - 1
            dMf[i,j] = np.linalg.norm(Symk-Sym[k])
    return dMf
    
def MultiFormLv(Sym):
    NumSym = len(Sym)
    Mf = np.zeros((NumSym,NumSym),int)
    for i in range(NumSym):
        for j in range(NumSym):
            Oi = Sym[i,0:3,0:3]; ti = Sym[i,0:3,3]
            Oj = Sym[j,0:3,0:3]; tj = Sym[j,0:3,3]
            Ok = Oi @ Oj
            tk = Oi @ tj + ti
            Symk = np.zeros((3,4))
            Symk[0:3,0:3] = Ok; Symk[0:3,3] = tk
            Indk = -1
            for k in range(NumSym):
                if np.linalg.norm(Symk-Sym[k]) < 1e-6:
                    Indk = k + 1
                    break
            Mf[i,j] = Indk
    return Mf

def MultiForm(Sym):
    NumSym = len(Sym)
    Mf = np.zeros((NumSym,NumSym),int)
    for i in range(NumSym):
        for j in range(NumSym):
            Symk = Sym[i] @ Sym[j]
            Indk = -1
            for k in range(NumSym):
                if np.linalg.norm(Symk-Sym[k]) < 1e-6:
                    Indk = k + 1
                    break
            Mf[i,j] = Indk
    return Mf
            
