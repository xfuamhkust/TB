import numpy as np
pi = np.pi


def GetEigStrip(ParaStrip,nk = 1001):
    # This function is to calculate the energy band of a quasi 1D strip
    
    # Parameters
    H00        = ParaStrip["Hmt00"]
    H01        = ParaStrip["Hmt01"]
    
    # Energy
    Kpt = np.linspace(-1,1,nk) * pi
    
    # Calculate energy band of strip.
    NumStt = len(H00)
    EigStrip = np.zeros((nk,NumStt))
    for ik in range(nk):
        H01k = H01 * np.exp(1j*Kpt[ik])
        Hk = H00 + H01k + H01k.conj().T
        EigStrip[ik] = np.linalg.eigvalsh(Hk)
    
    # Output
    ParaEigStrip = {"Eig": EigStrip,
                    }
    
    return ParaEigStrip



