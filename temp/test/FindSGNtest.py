import numpy as np

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

def CheckSymOrd2(SGi,SGU,SGID,Order):
    NumSG,NumID = SGID.shape
    count = 0
    Ind = np.ones(NumSG,bool)
    for j in Order:
        if np.var(SGID[:,j]) < 1e-3:
            continue
        count += 1
        flag = 0
        for SGii in SGi:
            if np.linalg.norm(SGii-SGU[j]) < 1e-3:
                flag = 1
                break
        IndInvj = np.where(1-(SGID[:,j] == flag)*Ind)
        Ind[IndInvj] = False
        if sum(Ind) == 1:
            return count, np.where(Ind)[0][0] + 1
        # print(sum(Ind))

SGID  = np.load("SGID.npy")
SGU   = np.load("SGU.npy")
Order = np.argsort(np.sum(SGID,axis=0))[::-1]
Count = np.zeros(230,int)
SGN   = np.zeros(230,int)
for i in range(1,230+1):
    SGi = GetSymLvSG(i)
    Count[i-1], SGN[i-1] = CheckSymOrd2(SGi,SGU,SGID,Order)
print(SGN)
print(np.mean(Count))
print(np.max(Count))