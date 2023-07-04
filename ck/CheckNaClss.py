import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['figure.dpi']= 300
pi2 = 2*np.pi

def Ek(EA,EB,t,k):
    dE = np.sqrt(((EA-EB)/2)**2 + (4*t**2)*(np.cos(k[0]/2)+np.cos(k[1]/2)+np.cos(k[2]/2))**2)
    E = np.array([(EA+EB)/2-dE,(EA+EB)/2+dE])
    return E 

EA = 1.0
EB = 0.0
t  = 0.3
KptXyz = np.load("KptXyzNaClss.npy")
EigVal0 = np.load("EigValNaClss.npy")
NumKpt = len(KptXyz)
EigVal = np.zeros((NumKpt,2))
for iKpt in range(NumKpt):
    EigVal[iKpt] = Ek(EA,EB,t,KptXyz[iKpt])

for iStt in range(2):
    plt.plot(EigVal[:,iStt], "r",lw=2)
    plt.plot(EigVal0[:,iStt],"b",lw=1)
    
print(np.sum(abs(EigVal-EigVal0)))