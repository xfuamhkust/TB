import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['figure.dpi']= 300
pi2 = 2*np.pi

def Ek(E0,t,k):
    dE = t * np.sqrt(3+2*np.cos(pi2*k[0])+2*np.cos(pi2*k[1])+2*np.cos(pi2*k[0]+pi2*k[1]))
    E = np.array([E0-dE,E0+dE])
    return E 

E0 = 0
t  = 1
KptLv  = np.load("KptLv.npy")
EigVal0 = np.load("EigVal.npy")
NumKpt = len(KptLv)
EigVal = np.zeros((NumKpt,2))
for iKpt in range(NumKpt):
    EigVal[iKpt] = Ek(E0,t,KptLv[iKpt])

for iStt in range(2):
    plt.plot(EigVal[:,iStt], "r",lw=2)
    plt.plot(EigVal0[:,iStt],"b",lw=1)
    
print(np.sum(abs(EigVal-EigVal0)))