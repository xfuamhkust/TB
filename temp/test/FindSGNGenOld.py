import numpy as np

def GetSymLvSG(SGN):
    FileName = "SpGrpMatGen.txt"
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

# NumSG = []
# for i in range(1,230+1):
#     NumSG.append(len(GetSymLvSG(i)))

# print(NumSG)
# print(sum(NumSG))

# SGall = np.zeros((988,3,4))
# count = 0
# for i in range(1,230+1):
#     SGi = GetSymLvSG(i)
#     for j in range(NumSG[i-1]):
#         SGall[count] = SGi[j]
#         count += 1

# SGunq = np.unique(SGall,axis=0)
# for i in range(988):
#     flag = 0
#     for j in range(83):
#         if np.linalg.norm(SGall[i] - SGunq[j])<1e-3:
#             flag = 1
#     if flag == 0:
#         print(i)
# print(len(np.unique(SGall,axis=0)))
    


NumSG0   = 230
NumSGall = 988
NumSGunq = 83

# Collect all space group symmetries
NumSG = np.zeros(230,int)
SGall = np.zeros((988,3,4))
count = 0
for i in range(1,230+1):
    SGi = GetSymLvSG(i)
    NumSG[i-1] = len(SGi)
    for j in range(NumSG[i-1]):
        SGall[count] = SGi[j]
        count += 1

# Find unique elements
SGunq = np.zeros((83,3,4))
Indunq = np.zeros(988,int)
count = 0
for i in range(988):
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
# NumSGindinv = np.zeros(83,int)
# SGindinv = []
# for i in range(83):
#     SGindinvi = []
#     for j in range(1,230+1):
#         if i in SGind[j-1]:
#             SGindinvi.append(j)
#     SGindinv.append(SGindinvi)
#     NumSGindinv[i] = len(SGindinvi)

# Get ID for each space group
SGID = np.zeros((230,83),int)
for i in range(1,230+1):
    SGID[i-1][SGind[i-1]] = 1

# Reduce ID
IfKeep = np.ones(83,int)
for i in range(83):
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



    
