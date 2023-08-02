import numpy as np
import sys
import os

def GetSymLvSG(SGN):
    FileName = "SpGrpMat.txt"
    f = [line.strip() for line in open(FileName,"r").readlines()]
    for iff in range(len(f)):
        fispl = f[iff].split(" ")
        if fispl[0] == "_" + "%d"%SGN + "_":
            NumMat = int(fispl[1])
            SGM = np.zeros((NumMat,3,4))
            for iMat in range(NumMat):
                fiispl = f[iff+iMat+1].split(" ")
                MatLine = []
                for iNum in range(12):
                    if fiispl[iNum].find("/") + 1:
                        [stri1,stri2] = fiispl[iNum].split("/")
                        MatLine.append(float(stri1)/float(stri2))
                    else:
                        MatLine.append(int(fiispl[iNum]))
                SGM[iMat] = np.array(MatLine).reshape((3,4))
            break
    return SGM

NumSG0   = 230
NumSGall = 4425
NumSGunq = 786

# Collect all space group symmetries
NumSG = np.zeros(230,int)
SGall = np.zeros((4425,3,4))
count = 0
for i in range(1,230+1):
    SGi = GetSymLvSG(i)
    NumSG[i-1] = len(SGi)
    for j in range(NumSG[i-1]):
        SGall[count] = SGi[j]
        count += 1

# Find unique elements
SGunq = np.zeros((786,3,4))
Indunq = np.zeros(4425,int)
count = 0
for i in range(4425):
    flag = 0
    for j in range(count):
        if np.linalg.norm(SGunq[j] - SGall[i]) < 1e-3:
            Indunq[i] = j
            flag = 1
            break
    if flag == 0:
        SGunq[count] = SGall[i]
        Indunq[i] = count
        count += 1

# Process Indunq
SGind = []
count = 0
for i in range(1,230+1):
    SGind.append(Indunq[count:count+NumSG[i-1]].tolist())
    count += NumSG[i-1]

# Get ID for each space group
SGID = np.zeros((230,786),int)
for i in range(1,230+1):
    SGID[i-1][SGind[i-1]] = 1

# Reduce ID
IfKeep = np.ones(786,int)
for i in range(786):
    i = 83-1-i
    SGIDtemp = np.copy(SGID)
    SGIDtemp[:,i] = 0
    if len(np.unique(SGIDtemp,axis=0)) == 230:
        IfKeep[i] = 0
        SGID = SGIDtemp


SGIDred = SGID[:,np.where(IfKeep)[0]]
SGunqred = SGunq[np.where(IfKeep)[0]]
print(SGIDred.shape)
print(SGunqred.shape)

'''
Write & Read
'''
np.save("SGID",np.array(SGIDred,"int8"))
np.save("SGU",SGunqred)
SGID = np.load("SGID.npy")
SGU  = np.load("SGU.npy")
'''
'''

def CheckSymAll(SGIDi,SGID):
    return np.where(np.sum(abs(SGID-SGIDi),axis=1)==0)[0][0] + 1

def CheckSymOrd(SGIDi,SGID,Order):
    NumSG,NumID = SGID.shape
    count = 0
    Ind = np.ones(NumSG,bool)
    for j in Order:
        if np.var(SGID[:,j]) < 1e-3:
            continue
        count += 1
        IndInvj = np.where(1-(SGID[:,j] == SGIDi[j])*Ind)
        Ind[IndInvj] = False
        if sum(Ind) == 1:
            return count, np.where(Ind)[0][0] + 1
        

# # Check all unique symmetries
# SGN = np.zeros(230,int)
# for i in range(1,230+1):
#     SGN[i-1] = CheckSymAll(SGID[i-1],SGID)
# print(SGN)

    
Order = np.argsort(np.sum(SGID,axis=0))[::-1]
Count = np.zeros(230,int)
SGN = np.zeros(230,int)
for i in range(1,230+1):
    Count[i-1], SGN[i-1] = CheckSymOrd(SGID[i-1],SGID,Order)
print(SGN)
print(np.mean(Count))
print(np.max(Count))