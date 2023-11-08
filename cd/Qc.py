import numpy as np
pi = np.pi

'''
This is a code to calculate the Quantum Conductance of a two-port system by Green's function.

x↓ y→
                     |←---------------L-------------→|
 ________________     _______________________________     __________________
                 |   |                               |   |                   ↑
                 |   |                               |   |                   |
   Left Lead     |   |        Center Material        |   |     Right Lead    W
                 |   |                               |   |                   |
 ________________|   |_______________________________|   |__________________ ↓
  
 □   □   □   □   □   □   □   □   □   □   □   □   □   □   □   □   □   □   □   □  
 
 □   □   □   □   □   □   □   □   □   □   □   □   □   □   □   □   □   □   □   □  
 
 □   □   □   □   □   □   □   □   □   □   □   □   □   □   □   □   □   □   □   □ 

                 □                           □                   □
T0: □        Tx: ↓       Ty: □ → □      Txy:   ↘        Tyx:   ↗
                 □                               □           □

 H10 H00 H01
 --- --- ---
| □ | □ | □ |          | T0   Tx       |          | Ty   Txy      |
|   |   |   |          |               |          |               |
| □ | □ | □ |    H00 = | Tx+  T0   Tx  |    H01 = | Tyx  Ty   Txy |
|   |   |   |          |               |          |               |
| □ | □ | □ |          |      Tx+  T0  |          |      Tyx  Ty  |
 --- --- ---

Only 2D system is available now.
It is assumed that left and right leads are the same material as the center.
The strip finitely extends along x-direction and infinitely extends along y-direction.
The x and y directions do not need to be perpendicular.
The number of cells (□) in x-direction is W, and that of Center Material in x direction is L.
The onsite hopping matrix of each cell is T0,
    and the offsite hopping matrix to the nearest right (down) is Ty (Tx),
    and the offsite hopping matrix to the nearest lower right (upper right) is Txy (Tyx).
The onsite hopping of each column is H00, 
    and the offsite hopping matrix to the nearest right column is H01.
The method to calculate the Surface Green's function can be chosen as Iteration or Eigsolution.

Ref:
    https://zhuanlan.zhihu.com/p/269595149
    https://mattermodeling.stackexchange.com/questions/4114/about-the-surface-greens-function-method-for-calculating-the-surface-state
'''

def GetQc(ParaQcIn,ParaStrip):
    
    # Parameters
    Method     = ParaQcIn["Method"]
    nE         = ParaQcIn["nE"]
    Emin       = ParaQcIn["Emin"]
    Emax       = ParaQcIn["Emax"]
    eta        = ParaQcIn["eta"]
    Len        = ParaQcIn["Len"]
    H00        = ParaStrip["Hmt00"]
    H01        = ParaStrip["Hmt01"]
    
    # Energy
    E0 = np.linspace(Emin,Emax,nE)
    
    # Calculate quantum conductance by Green's function.
    TLR = GetQuanCond(H00,H01,Len,E0,eta,Method)
    
    # Output
    ParaQc = {"Ener":  E0,
              "TLR":   TLR,
              }
    
    return ParaQc

def GetQuanCond(H00,H01,L,E0,eta=1e-6,Method="Iteration"):
    # Preparation
    nE = len(E0)         # The number of energy points
    E  = E0 + 1j*eta     # Energy with imaginary part
    H10 = H01.conj().T   # The Hermitian conjugate of H01
    nstt = len(H00)    # The number of states in one column
    if Method == "Eigsolution" and np.linalg.det(H01) == 0:
        Method="Iteration"
        print("The method Eigsolution is unavailable due to invertible matrix H01.")
    # For loop (E)
    TLR  = np.zeros(nE)  # Transmission
    for i in range(nE):
        n_cut = 50
        EI = E[i]*np.eye(nstt)
        # Surface Green's function
        if Method == "Iteration":
            G00_R = GetG00Iter(EI, H00, H01, H10, n_cut)
            G00_L = GetG00Iter(EI, H00, H10, H01, n_cut)
        elif Method == "Eigsolution":
            G00_R = GetG00Eig( EI, H00, H01, H10)
            G00_L = GetG00Eig( EI, H00, H10, H01)
        # Self energy & Level width function
        SGM_R = H01 @ G00_R @ H10
        SGM_L = H10 @ G00_L @ H01
        GMA_R = (SGM_R-SGM_R.conj().T)*1j
        GMA_L = (SGM_L-SGM_L.conj().T)*1j
        # G1M
        Gi = np.linalg.inv(EI - H00 - SGM_L)
        G1M = np.copy(Gi)
        for j in range(2,L):
            Gi = np.linalg.inv(EI - H00 - H10 @ Gi @ H01)
            G1M = G1M @ H01 @ Gi
        Gi = np.linalg.inv(EI - H00 - H10 @ Gi @ H01 - SGM_R)
        G1M = G1M @ H01 @ Gi
        # Transmission 
        TLR[i] = np.real(np.trace(GMA_L @ G1M @ GMA_R @ G1M.conj().T))
    return TLR

def GetG00Iter(EI, H00, H01, H10, n_cut):
    Ha = np.copy(H00)
    Hb = np.copy(H00)
    A  = np.copy(H01)
    B  = np.copy(H10)
    for j in range(n_cut):
        g   = np.linalg.inv(EI - Hb)
        AgB = A @ g @ B
        BgA = B @ g @ A
        Ha  = Ha + AgB
        Hb  = Hb + AgB + BgA
        A   = A @ g @ A
        B   = B @ g @ B
    G00 = np.linalg.inv(EI - Ha)
    return G00

def GetG00Eig(EI, H00, H01, H10):
    nstt = len(EI)
    H01_inv = np.linalg.inv(H01)
    T = np.zeros((2*nstt,2*nstt),complex)
    T[0*nstt:1*nstt,0*nstt:1*nstt] =  H01_inv @ (EI - H00)
    T[0*nstt:1*nstt,1*nstt:2*nstt] = -H01_inv @ H10 # Mistake of Taylor cat's zhihu-blog 
    T[1*nstt:2*nstt,0*nstt:1*nstt] =  np.eye(nstt)
    eig_val, eig_stt = np.linalg.eig(T)
    Ind1 = np.where(abs(eig_val)<1)[0]
    S1 = eig_stt[0*nstt:1*nstt,Ind1]; S2 = eig_stt[1*nstt:2*nstt,Ind1]
    G00 = np.linalg.inv(EI - H00 - H01 @ S1 @ np.linalg.inv(S2))
    return G00

def FindIndex(x,X,tol = 1e-3):
    dX = np.sum(abs(X - x), axis=1)
    min_dX = min(dX)
    if min_dX < tol:
        ind = np.where(dX == min_dX)[0][0]
        return ind
    else:
        return -1

