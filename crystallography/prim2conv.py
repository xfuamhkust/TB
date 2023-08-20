import numpy as np
import math


# this script transforms a primitive cell to a conventional cell
# TODO: values of errors (epsilon) to be rectified

def inputParametersInDegree2S(a, b, c, alpha, beta, gamma):
    """

    :param a: length
    :param b: length
    :param c: length
    :param alpha: angle between b and c, in degrees
    :param beta: angle between c and a, in degrees
    :param gamma: angle between a and b, in degrees
    :return: matrix S
    """

    # sanity check
    if a <= 0 or b <= 0 or c <= 0 or alpha <= 0 or alpha >= 180 or beta <= 0 or beta >= 180 or gamma <= 0 or gamma >= 180:
        raise ValueError("Invalid input: out of range.")

    # degree to radian
    alphaInRadian = math.radians(alpha)
    betaInRadians = math.radians(beta)
    gammaInRadians = math.radians(gamma)

    # sanity check
    minAlpha = np.arccos(np.cos(betaInRadians - gammaInRadians))
    maxAlpha = np.arccos(np.cos(betaInRadians + gammaInRadians))
    if alphaInRadian < minAlpha or alphaInRadian > maxAlpha:
        raise ValueError("Invalid combination of angles.")

    S = np.array([[a ** 2, b ** 2, c ** 2],
                  [2 * b * c * np.cos(alphaInRadian), 2 * a * c * np.cos(betaInRadians),
                   2 * a * b * np.cos(gammaInRadians)]], dtype=np.float64)

    return S


def inputParametersInRadian2S(a, b, c, alpha, beta, gamma):
    """

       :param a: length
       :param b: length
       :param c: length
       :param alpha: angle between b and c, in radians
       :param beta: angle between c and a, in radians
       :param gamma: angle between a and b, in radians
       :return: matrix S
       """
    # sanity check
    if a <= 0 or b <= 0 or c <= 0 or alpha <= 0 or alpha >= np.pi or beta <= 0 or beta >= np.pi or gamma <= 0 or gamma >= np.pi:
        raise ValueError("Invalid input: out of range.")
    # sanity check
    minAlpha = np.arccos(np.cos(beta - gamma))
    maxAlpha = np.arccos(np.cos(beta + gamma))
    if alpha < minAlpha or alpha > maxAlpha:
        raise ValueError("Invalid combination of angles.")

    S = np.array([[a ** 2, b ** 2, c ** 2],
                  [2 * b * c * np.cos(alpha), 2 * a * c * np.cos(beta),
                   2 * a * b * np.cos(gamma)]], dtype=np.float64)

    return S


def checkAnglesOfS(S):
    """

    :param S: matrix S=[[A, B, C],[xi, eta, zeta]]
    :return: check if angles are in valid combinations
    """
    [[A, B, C], [xi, eta, zeta]] = S
    a = np.sqrt(A)
    b = np.sqrt(B)
    c = np.sqrt(C)

    alpha = np.arccos(xi / (2 * b * c))
    beta = np.arccos(eta / (2 * a * c))
    gamma = np.arccos(zeta / (2 * a * b))
    if np.isnan(alpha) or np.isnan(beta) or np.isnan(gamma):
        raise ValueError("Matrix S with invalid combination of angles.")
    # sanity check
    minAlpha = np.arccos(np.cos(beta - gamma))
    maxAlpha = np.arccos(np.cos(beta + gamma))
    if alpha < minAlpha or alpha > maxAlpha:
        raise ValueError("Matrix S with invalid combination of angles.")


def AlgorithmN(S):
    """

    :param S: matrix S=[[A, B, C],[xi, eta, zeta]]
    :return: matrix S with normalized parameters
    """
    # [[A, B, C], [xi, eta, zeta]] = S
    # a = np.sqrt(A)
    # b = np.sqrt(B)
    # c = np.sqrt(C)
    # cosAlpha = xi / (2 * b * c)
    # cosBeta = eta / (2 * a * c)
    # cosGamma = zeta / (2 * a * b)
    # print("Before N cos(alpha)=" + str(cosAlpha) + ", cos(beta)=" + str(cosBeta) + ", cos(gamma)=" + str(cosGamma))

    epsRel = 1e-8
    epsAbs = 1e-6
    # step 1, sort A B C
    ABCInd = np.argsort(S[0, :])

    S = S[:, ABCInd]

    # step 2, sort xi and eta

    if np.isclose(S[0, 0], S[0, 1], rtol=epsRel, atol=epsAbs, equal_nan=False):
        xiEtaInd = np.argsort(np.abs(S[1, [0, 1]]))
        S[1, [0, 1]] = S[1, xiEtaInd]

    # step 3, sort eta and zeta
    if np.isclose(S[0, 1], S[0, 2], rtol=epsRel, atol=epsAbs, equal_nan=False):
        etaZetaInd = np.argsort(
            np.abs(S[1, [1, 2]])) + 1  # +1 because the actual indices are 1 and 2, but argsort returns 0 and 1
        S[1, [1, 2]] = S[1, etaZetaInd]

    # step 4, giving signs to xi, eta, zeta
    if np.isclose(S[1, 0], 0, rtol=epsRel, atol=epsAbs, equal_nan=False):
        S[1, :] = -np.abs(S[1, :])
        checkAnglesOfS(S)
        return S
    elif np.isclose(S[1, 1], 0, rtol=epsRel, atol=epsAbs, equal_nan=False):
        S[1, :] = -np.abs(S[1, :])
        checkAnglesOfS(S)
        return S
    elif np.isclose(S[1, 2], 0, rtol=epsRel, atol=epsAbs, equal_nan=False):
        S[1, :] = -np.abs(S[1, :])
        checkAnglesOfS(S)
        return S

    if S[1, 0] * S[1, 1] * S[1, 2] > 0:
        S[1, :] = np.abs(S[1, :])
        checkAnglesOfS(S)
        return S
    else:
        S[1, :] = -np.abs(S[1, :])
        checkAnglesOfS(S)
        return S
    # [[A, B, C], [xi, eta, zeta]] = S
    # a = np.sqrt(A)
    # b = np.sqrt(B)
    # c = np.sqrt(C)
    # cosAlpha = xi / (2 * b * c)
    # cosBeta = eta / (2 * a * c)
    # cosGamma = zeta / (2 * a * b)
    # print("After N cos(alpha)="+str(cosAlpha)+", cos(beta)="+str(cosBeta)+", cos(gamma)="+str(cosGamma))
    # print("N: "+str(S[0,0])+", "+str(S[0,1])+", "+str(S[0,2])+", "+str(S[1,0])+", "+str(S[1,1])+", "+str(S[1,2]))
    # return S


def vec2S(aVec, bVec, cVec):
    """

        :param aVec: basis vector of primitive cell
        :param bVec: basis vector of primitive cell
        :param cVec: basis vector of primitive cell
        :return: S matrix
        """
    epsRel = 1e-8
    epsAbs = 1e-6

    A = aVec.dot(aVec)
    B = bVec.dot(bVec)
    C = cVec.dot(cVec)
    if np.isclose(A, 0, rtol=epsRel, atol=epsAbs, equal_nan=False):
        raise ValueError("vector a is 0.")
    if np.isclose(B, 0, rtol=epsRel, atol=epsAbs, equal_nan=False):
        raise ValueError("vector b is 0.")
    if np.isclose(C, 0, rtol=epsRel, atol=epsAbs, equal_nan=False):
        raise ValueError("vector c is 0.")

    xi = 2 * bVec.dot(cVec)
    eta = 2 * aVec.dot(cVec)
    zeta = 2 * aVec.dot(bVec)

    S = np.array([[A, B, C], [xi, eta, zeta]], dtype=np.float64)

    return S


def vec2xi(aVec, bVec, cVec):
    """

        :param aVec: basis vector of primitive cell
        :param bVec: basis vector of primitive cell
        :param cVec: basis vector of primitive cell
        :return: xi=2bVec . cVec
        """
    return 2 * bVec.dot(cVec)


def vec2eta(aVec, bVec, cVec):
    """

        :param aVec: basis vector of primitive cell
        :param bVec: basis vector of primitive cell
        :param cVec: basis vector of primitive cell
        :return: eta=2aVec . cVec
        """
    return 2 * aVec.dot(cVec)


def vec2zeta(aVec, bVec, cVec):
    """

        :param aVec: basis vector of primitive cell
        :param bVec: basis vector of primitive cell
        :param cVec: basis vector of primitive cell
        :return: zeta=2aVec . bVec
        """
    return 2 * aVec.dot(bVec)


def find0(xiEtaZeta):
    """
    :param xiEtaZeta=[xi,eta,zeta]
    :param xi: real number
    :param eta: real number
    :param zeta: real number
    :return: return indices of 0
    """
    epsRel = 1e-8
    epsAbs = 1e-6
    is0 = np.array([np.isclose(tmp, 0, rtol=epsRel, atol=epsAbs, equal_nan=False) for tmp in xiEtaZeta])
    return np.where(is0 == True)[0]


# print(find0([0,0,0]))

# aVec=np.array([1,2,3])
#
# bVec=np.array([0.1,0.2,0.3])
# cVec=np.array([0.5,0.5,0.5])
# print(bVec)
def findNegative(xiEtaZeta):
    """
        :param xiEtaZeta=[xi,eta,zeta]
        :param xi: real number
        :param eta: real number
        :param zeta: real number
        :return: return indices of <0 elems
        """
    isNeg = np.array([tmp < 0 for tmp in xiEtaZeta])
    return np.where(isNeg == True)[0]


# print(findNegative([0,0,-1]))
def findPositive(xiEtaZeta):
    """
        :param xiEtaZeta=[xi,eta,zeta]
        :param xi: real number
        :param eta: real number
        :param zeta: real number
        :return: return indices of >0 elems
        """
    isNeg = np.array([tmp > 0 for tmp in xiEtaZeta])
    return np.where(isNeg == True)[0]


# print(findPositive([0,1,1]))
def changeSign(aVec, bVec, cVec, j):
    """

    :param aVec:
    :param bVec:
    :param cVec:
    :param j: 0,1,2
    :return: change of  the directions of aVec, bVec, cVec
    """
    if j == 0:
        aVec = -aVec
    elif j == 1:
        bVec = -bVec
    elif j == 2:
        cVec = -cVec
    else:
        raise ValueError("Wrong value of j.")
    return aVec, bVec, cVec


# aVec,bVec,cVec=changeSign(aVec,bVec,cVec,1)
# print(bVec)


import copy


def AlgorithmNVectorTransform(aVec, bVec, cVec):
    """

    :param aVec: basis vector of primitive cell
    :param bVec: basis vector of primitive cell
    :param cVec: basis vector of primitive cell
    :return: basis vectors whose matrix S is normalized
    """
    epsRel = 1e-8
    epsAbs = 1e-6
    S = vec2S(aVec, bVec, cVec)
    # step 1, sort A B C
    j0, j1, j2 = np.argsort(S[0, :])
    originalVecs = [copy.deepcopy(aVec), copy.deepcopy(bVec), copy.deepcopy(cVec)]
    aVec = originalVecs[j0]
    bVec = originalVecs[j1]
    cVec = originalVecs[j2]

    # print(aVec)
    # print(bVec)
    # print(cVec)

    # step 2, sort xi and eta

    if np.isclose(S[0, 0], S[0, 1], rtol=epsRel, atol=epsAbs, equal_nan=False):
        xi = vec2xi(aVec, bVec, cVec)
        eta = vec2eta(aVec, bVec, cVec)
        if np.abs(xi) > np.abs(eta):
            aVec, bVec = bVec, aVec
        # step 3, sort eta and zeta
    if np.isclose(S[0, 1], S[0, 2], rtol=epsRel, atol=epsAbs, equal_nan=False):
        eta = vec2eta(aVec, bVec, cVec)
        zeta = vec2zeta(aVec, bVec, cVec)
        if np.abs(eta) > np.abs(zeta):
            bVec, cVec = cVec, bVec

    xi = vec2xi(aVec, bVec, cVec)
    eta = vec2eta(aVec, bVec, cVec)
    zeta = vec2zeta(aVec, bVec, cVec)
    SRow1 = np.array([xi, eta, zeta])
    indsOf0 = find0(SRow1)
    # step 4, giving signs to xi, eta, zeta
    # case (3): xi*eta*zeta=0
    ##(a) xi=eta=zeta=0
    if len(indsOf0) == 3:
        return aVec, bVec, cVec
    ##(b) 2 of xi, eta, zeta are 0
    elif len(indsOf0) == 2:
        j1, j2 = indsOf0
        j0 = list({0, 1, 2} - {j1, j2})[0]
        if SRow1[j0] > 0:
            aVec, bVec, cVec = changeSign(aVec, bVec, cVec, j1)
            return aVec, bVec, cVec
        else:
            return aVec, bVec, cVec
    ##(c) 1 of xi,eta,zeta is 0
    elif len(indsOf0) == 1:
        j2 = indsOf0[0]
        j0, j1 = list({0, 1, 2} - {j2})
        ### (i)
        if SRow1[j0] > 0 and SRow1[j1] > 0:
            aVec, bVec, cVec = changeSign(aVec, bVec, cVec, j2)
            return aVec, bVec, cVec
        ### (ii)
        if SRow1[j0] < 0 and SRow1[j1] < 0:
            return aVec, bVec, cVec
        ### (iii)
        if SRow1[j0] < 0 and SRow1[j1] > 0:
            j3 = j0
            aVec, bVec, cVec = changeSign(aVec, bVec, cVec, j3)
            return aVec, bVec, cVec
        ### (iii)
        if SRow1[j0] > 0 and SRow1[j1] < 0:
            j3 = j1
            aVec, bVec, cVec = changeSign(aVec, bVec, cVec, j3)
            return aVec, bVec, cVec

    # case (2): xi*eta*zeta<0
    if xi * eta * zeta < 0:
        indsNeg = findNegative(SRow1)
        ## (a): all 3 of xi, eta, zeta are <0
        if len(indsNeg) == 3:
            return aVec, bVec, cVec
        ## (b): 1 of xi, eta, zeta is <0
        elif len(indsNeg) == 1:
            j0 = indsNeg[0]
            aVec, bVec, cVec = changeSign(aVec, bVec, cVec, j0)
            return aVec, bVec, cVec
        else:
            raise RuntimeError("Invalid xi, eta , zeta values for xi*eta*zeta<0.")

    # case (1): xi*eta*zeta>0
    if xi * eta * zeta > 0:
        indsPos = findPositive(SRow1)
        ## (a) all 3 of xi, eta, zeta are >0
        if len(indsPos) == 3:
            return aVec, bVec, cVec
        ## (b) only 1 of xi, eta, zeta are>0
        elif len(indsPos) == 1:
            j2 = indsPos[0]
            aVec, bVec, cVec = changeSign(aVec, bVec, cVec, j2)
            return aVec, bVec, cVec
        else:
            raise RuntimeError("Invalid xi, eta , zeta values for xi*eta*zeta>0.")


#
# aVec=np.array([1,2,3])
#
# bVec=np.array([0.1,0.2,0.3])
# cVec=np.array([0.5,0.5,0.5])
#
#
# AlgorithmNVectorTransform(aVec,bVec,cVec)

def AlgorithmB2VectorTransform(aVec, bVec, cVec):
    """

        :param aVec: basis vector of primitive cell
        :param bVec: basis vector of primitive cell
        :param cVec: basis vector of primitive cell
        :return: basis vectors after algorithm B2
        """

    S = vec2S(aVec, bVec, cVec)
    [[_, B, C], [xi, eta, zeta]] = S
    # in python, entier()=math.floor()
    j = math.floor((xi + B) / (2 * B))
    cVec = cVec - j * bVec
    aVec, bVec, cVec = AlgorithmNVectorTransform(aVec, bVec, cVec)
    return aVec, bVec, cVec


def AlgorithmB3VectorTransform(aVec, bVec, cVec):
    """

    :param aVec: basis vector of primitive cell
    :param bVec: basis vector of primitive cell
    :param cVec: basis vector of primitive cell
    :return: basis vectors after algorithm B3
    """
    S = vec2S(aVec, bVec, cVec)
    [[A, _, C], [xi, eta, zeta]] = S
    j = math.floor((eta + A) / (2 * A))
    cVec = cVec - j * aVec
    aVec, bVec, cVec = AlgorithmNVectorTransform(aVec, bVec, cVec)

    return aVec, bVec, cVec


def AlgorithmB4VectorTransform(aVec, bVec, cVec):
    """

    :param aVec: basis vector of primitive cell
    :param bVec: basis vector of primitive cell
    :param cVec: basis vector of primitive cell
    :return: basis vectors after algorithm B4
    """
    S = vec2S(aVec, bVec, cVec)
    [[A, B, _], [xi, eta, zeta]] = S

    j = math.floor((zeta + A) / (2 * A))
    bVec = bVec - j * aVec
    aVec, bVec, cVec = AlgorithmNVectorTransform(aVec, bVec, cVec)

    return aVec, bVec, cVec


def AlgorithmB5VectorTransform(aVec, bVec, cVec):
    """

    :param aVec: basis vector of primitive cell
    :param bVec: basis vector of primitive cell
    :param cVec: basis vector of primitive cell
    :return: basis vectors after algorithm B4
    """
    S = vec2S(aVec, bVec, cVec)
    [[A, B, C], [xi, eta, zeta]] = S
    j = math.floor((xi + eta + zeta + A + B) / (2 * (A + B + zeta)))

    cVec = cVec - j * (aVec + bVec)
    aVec, bVec, cVec = AlgorithmNVectorTransform(aVec, bVec, cVec)

    return aVec, bVec, cVec


# assemble algorithm B's linear transformation form

def AlgorithmBVectorTransform(aVec, bVec, cVec):
    """

    :param aVec: basis vector of primitive cell
    :param bVec: basis vector of primitive cell
    :param cVec: basis vector of primitive cell
    :return: basis of a Buerger cell
    """
    aVec, bVec, cVec = AlgorithmNVectorTransform(aVec, bVec, cVec)
    S = vec2S(aVec, bVec, cVec)
    [[A, B, _], [xi, eta, zeta]] = S

    condB2 = (np.abs(xi) > B)
    condB3 = (np.abs(eta) > A)
    condB4 = (np.abs(zeta) > A)
    condB5 = (xi + eta + zeta + A + B < 0)

    while condB2 or condB3 or condB4 or condB5:
        S = vec2S(aVec, bVec, cVec)
        [[_, B, _], [xi, _, _]] = S
        condB2 = (np.abs(xi) > B)
        if condB2:
            aVec, bVec, cVec = AlgorithmB2VectorTransform(aVec, bVec, cVec)

        S = vec2S(aVec, bVec, cVec)
        [[A, _, _], [_, eta, _]] = S
        condB3 = (np.abs(eta) > A)
        if condB3:
            aVec, bVec, cVec = AlgorithmB3VectorTransform(aVec, bVec, cVec)

        S = vec2S(aVec, bVec, cVec)
        [[A, _, _], [_, _, zeta]] = S
        condB4 = (np.abs(zeta) > A)
        if condB4:
            aVec, bVec, cVec = AlgorithmB4VectorTransform(aVec, bVec, cVec)

        S = vec2S(aVec, bVec, cVec)
        [[A, B, _], [xi, eta, zeta]] = S
        condB5 = (xi + eta + zeta + A + B < 0)
        if condB5:
            aVec, bVec, cVec = AlgorithmB5VectorTransform(aVec, bVec, cVec)

        S = vec2S(aVec, bVec, cVec)

        [[A, B, _], [xi, eta, zeta]] = S
        condB2 = (np.abs(xi) > B)
        condB3 = (np.abs(eta) > A)
        condB4 = (np.abs(zeta) > A)
        condB5 = (xi + eta + zeta + A + B < 0)

    return aVec, bVec, cVec


# components of Algorithm B
def AlgorthmB2(S):
    """

    :param S:
    :return:
    """
    [[_, B, C], [xi, eta, zeta]] = S
    # in python, entier()=math.floor()
    j = math.floor((xi + B) / (2 * B))
    C = C + j ** 2 * B - j * xi
    xi = xi - 2 * j * B
    eta = eta - j * zeta
    S[0, 2] = C
    S[1, 0] = xi
    S[1, 1] = eta
    # print("B2: " + str(S[0, 0]) + ", " + str(S[0, 1]) + ", " + str(S[0, 2]) + ", " + str(S[1, 0]) + ", " + str(
    #     S[1, 1]) + ", " + str(S[1, 2]))
    # normalization
    S = AlgorithmN(S)
    return S


def AlgorithmB3(S):
    """

    :param S:
    :return:
    """
    [[A, _, C], [xi, eta, zeta]] = S
    # a=np.sqrt(A)
    # b=np.sqrt(B)
    # c=np.sqrt(C)
    # cosAlpha=xi/(2*b*c)
    # cosBeta=eta/(2*a*c)
    # cosGamma=zeta/(2*a*b)
    # # print("Before B3: a="+str(a)+", b="+str(b)+", c="+str(c)+", xi="+str(xi)+", eta="+str(eta)+", zeta="+str(zeta) +", cos(alpha)="+str(cosAlpha)+", cos(beta)="+str(cosBeta)+", cos(gamma)="+str(cosGamma))

    j = math.floor((eta + A) / (2 * A))
    # CNew=c**2+j**2*a**2-j*2*a*c*cosBeta

    C = C + j ** 2 * A - j * eta
    # print("In B3: CNew=" + str(CNew) + ", C=" + str(C))
    xi = xi - j * zeta
    eta = eta - 2 * j * A

    S[0, 2] = C
    S[1, 0] = xi
    S[1, 1] = eta

    # [[A, B, C], [xi, eta, zeta]] = S
    # a = np.sqrt(A)
    # b = np.sqrt(B)
    # c = np.sqrt(C)
    # cosAlpha = xi / (2 * b * c)
    # cosBeta = eta / (2 * a * c)
    # cosGamma = zeta / (2 * a * b)
    # print("After B3: a="+str(a)+", b="+str(b)+", c="+str(c)+", xi="+str(xi)+", eta="+str(eta)+", zeta="+str(zeta) +", cos(alpha)="+str(cosAlpha)+", cos(beta)="+str(cosBeta)+", cos(gamma)="+str(cosGamma))

    # print("B3: " + str(S[0, 0]) + ", " + str(S[0, 1]) + ", " + str(S[0, 2]) + ", " + str(S[1, 0]) + ", " + str(
    #     S[1, 1]) + ", " + str(S[1, 2]))
    # normalization
    S = AlgorithmN(S)
    return S


def AlgorithmB4(S):
    """

    :param S:
    :return:
    """
    [[A, B, _], [xi, eta, zeta]] = S

    j = math.floor((zeta + A) / (2 * A))
    B = B + j ** 2 * A - j * zeta
    xi = xi - j * eta
    zeta = zeta - 2 * j * A

    S[0, 1] = B
    S[1, 0] = xi
    S[1, 2] = zeta

    # print("B4: " + str(S[0, 0]) + ", " + str(S[0, 1]) + ", " + str(S[0, 2]) + ", " + str(S[1, 0]) + ", " + str(
    #     S[1, 1]) + ", " + str(S[1, 2]))
    # normalization
    S = AlgorithmN(S)
    return S


def AlgorithmB5(S):
    """

    :param S:
    :return:
    """
    [[A, B, C], [xi, eta, zeta]] = S

    j = math.floor((xi + eta + zeta + A + B) / (2 * (A + B + zeta)))
    C = C + j ** 2 * (A + B + zeta) - j * (xi + eta)
    xi = xi - j * (2 * B + zeta)
    eta = eta - j * (2 * A + zeta)

    S[0, 2] = C
    S[1, 0] = xi
    S[1, 1] = eta

    # print("B5: " + str(S[0, 0]) + ", " + str(S[0, 1]) + ", " + str(S[0, 2]) + ", " + str(S[1, 0]) + ", " + str(
    #     S[1, 1]) + ", " + str(S[1, 2]))
    # normalization
    S = AlgorithmN(S)
    return S


# assemble Algorithm B
def AlgorithmB(S):
    """

    :param S:
    :return: matrix S containing parameters of a Buerger cell
    """
    S = AlgorithmN(S)

    [[A, B, _], [xi, eta, zeta]] = S

    condB2 = (np.abs(xi) > B)
    condB3 = (np.abs(eta) > A)
    condB4 = (np.abs(zeta) > A)
    condB5 = (xi + eta + zeta + A + B < 0)

    # num=0
    while condB2 or condB3 or condB4 or condB5:
        [[_, B, _], [xi, _, _]] = S
        condB2 = (np.abs(xi) > B)
        if condB2:
            S = AlgorthmB2(S)

        [[A, _, _], [_, eta, _]] = S
        condB3 = (np.abs(eta) > A)
        if condB3:
            S = AlgorithmB3(S)

        [[A, _, _], [_, _, zeta]] = S
        condB4 = (np.abs(zeta) > A)
        if condB4:
            S = AlgorithmB4(S)

        [[A, B, _], [xi, eta, zeta]] = S
        condB5 = (xi + eta + zeta + A + B < 0)
        if condB5:
            S = AlgorithmB5(S)

        [[A, B, _], [xi, eta, zeta]] = S
        condB2 = (np.abs(xi) > B)
        condB3 = (np.abs(eta) > A)
        condB4 = (np.abs(zeta) > A)
        condB5 = (xi + eta + zeta + A + B < 0)
        # num+=1
        # if num>10:
        #     break

    return S


def potentialBuergerCharacteristics(S):
    """

    :param S: returned from Algorithm B
    :return: all (unnormalized) Buerger characteristics belonging to the same Niggli cell
    """
    [[A, B, C], [xi, eta, zeta]] = S
    sumVal = A + B + C
    sigma = xi + eta + zeta

    epsRel = 1e-8
    epsAbs = 1e-6
    potentialBgChars = []
    # test condition 1
    A1 = A
    B1 = B
    C1 = C
    if np.isclose(A1 + B1 + C1, sumVal, rtol=epsRel, atol=epsAbs, equal_nan=False):
        xi1 = xi
        eta1 = eta
        zeta1 = zeta
        potentialBgChars.append(np.array([[A1, B1, C1], [xi1, eta1, zeta1]], dtype=np.float64))

    # test condition 2
    A2 = A
    B2 = B
    C2 = A + C + eta
    if np.isclose(A2 + B2 + C2, sumVal, rtol=epsRel, atol=epsAbs, equal_nan=False):
        xi2 = xi + zeta
        eta2 = 2 * A + eta
        zeta2 = zeta
        potentialBgChars.append(np.array([[A2, B2, C2], [xi2, eta2, zeta2]], dtype=np.float64))

    # test condition 3
    A3 = A
    B3 = B
    C3 = A + C - eta
    if np.isclose(A3 + B3 + C3, sumVal, rtol=epsRel, atol=epsAbs, equal_nan=False):
        xi3 = -xi + zeta
        eta3 = 2 * A - eta
        zeta3 = zeta
        potentialBgChars.append(np.array([[A3, B3, C3], [xi3, eta3, zeta3]], dtype=np.float64))

    # test condition 4
    A4 = A
    B4 = B
    C4 = B + C + xi
    if np.isclose(A4 + B4 + C4, sumVal, rtol=epsRel, atol=epsAbs, equal_nan=False):
        xi4 = 2 * B + xi
        eta4 = eta + zeta
        zeta4 = zeta
        potentialBgChars.append(np.array([[A4, B4, C4], [xi4, eta4, zeta4]], dtype=np.float64))

    # test condition 5
    A5 = A
    B5 = B
    C5 = B + C - xi
    if np.isclose(A5 + B5 + C5, sumVal, rtol=epsRel, atol=epsAbs, equal_nan=False):
        xi5 = 2 * B - xi
        eta5 = -eta + zeta
        zeta5 = zeta
        potentialBgChars.append(np.array([[A5, B5, C5], [xi5, eta5, zeta5]], dtype=np.float64))

    # test condition 6
    A6 = A
    B6 = B
    C6 = sumVal + sigma
    if np.isclose(A6 + B6 + C6, sumVal, rtol=epsRel, atol=epsAbs, equal_nan=False):
        xi6 = 2 * B + xi + zeta
        eta6 = 2 * A + eta + zeta
        zeta6 = zeta
        potentialBgChars.append(np.array([[A6, B6, C6], [xi6, eta6, zeta6]], dtype=np.float64))

    # test condition 7
    A7 = A
    B7 = C
    C7 = A + B + zeta
    if np.isclose(A7 + B7 + C7, sumVal, rtol=epsRel, atol=epsAbs, equal_nan=False):
        xi7 = xi + eta
        eta7 = 2 * A + zeta
        zeta7 = eta
        potentialBgChars.append(np.array([[A7, B7, C7], [xi7, eta7, zeta7]], dtype=np.float64))

    # test condition 8
    A8 = A
    B8 = C
    C8 = A + B - zeta
    if np.isclose(A8 + B8 + C8, sumVal, rtol=epsRel, atol=epsAbs, equal_nan=False):
        xi8 = -xi + eta
        eta8 = 2 * A - zeta
        zeta8 = eta
        potentialBgChars.append(np.array([[A8, B8, C8], [xi8, eta8, zeta8]], dtype=np.float64))

    # test condition 9
    A9 = A
    B9 = C
    C9 = B + C + xi
    if np.isclose(A9 + B9 + C9, sumVal, rtol=epsRel, atol=epsAbs, equal_nan=False):
        xi9 = 2 * C + xi
        eta9 = eta + zeta
        zeta9 = eta
        potentialBgChars.append(np.array([[A9, B9, C9], [xi9, eta9, zeta9]], dtype=np.float64))

    # test condition 10
    A10 = A
    B10 = C
    C10 = B + C - xi
    if np.isclose(A10 + B10 + C10, sumVal, rtol=epsRel, atol=epsAbs, equal_nan=False):
        xi10 = -2 * C + xi
        eta10 = -eta + zeta
        zeta10 = eta
        potentialBgChars.append(np.array([[A10, B10, C10], [xi10, eta10, zeta10]], dtype=np.float64))

    # test condition 11
    A11 = A
    B11 = C
    C11 = sumVal + sigma
    if np.isclose(A11 + B11 + C11, sumVal, rtol=epsRel, atol=epsAbs, equal_nan=False):
        xi11 = 2 * C + xi + eta
        eta11 = 2 * A + eta + zeta
        zeta11 = eta
        potentialBgChars.append(np.array([[A11, B11, C11], [xi11, eta11, zeta11]], dtype=np.float64))

    # test condition 12
    A12 = A
    B12 = A + B + zeta
    C12 = B + C + xi
    if np.isclose(A12 + B12 + C12, sumVal, rtol=epsRel, atol=epsAbs, equal_nan=False):
        xi12 = 2 * B + sigma
        eta12 = eta + zeta
        zeta12 = 2 * A + zeta
        potentialBgChars.append(np.array([[A12, B12, C12], [xi12, eta12, zeta12]], dtype=np.float64))

    # test condition 13
    A13 = A
    B13 = A + B + zeta
    C13 = sumVal + sigma
    if np.isclose(A13 + B13 + C13, sumVal, rtol=epsRel, atol=epsAbs, equal_nan=False):
        xi13 = 2 * A + 2 * B + sigma + zeta
        eta13 = 2 * A + eta + zeta
        zeta13 = 2 * A + zeta
        potentialBgChars.append(np.array([[A13, B13, C13], [xi13, eta13, zeta13]], dtype=np.float64))

    # test condition 14
    A14 = A
    B14 = A + B - zeta
    C14 = B + C - xi
    if np.isclose(A14 + B14 + C14, sumVal, rtol=epsRel, atol=epsAbs, equal_nan=False):
        xi14 = -2 * B + xi - eta + zeta
        eta14 = -eta + zeta
        zeta14 = 2 * A - zeta
        potentialBgChars.append(np.array([[A14, B14, C14], [xi14, eta14, zeta14]], dtype=np.float64))

    # test condition 15
    A15 = A
    B15 = A + C + eta
    C15 = sumVal + sigma
    if np.isclose(A15 + B15 + C15, sumVal, rtol=epsRel, atol=epsAbs, equal_nan=False):
        xi15 = 2 * A + 2 * C + sigma + eta
        eta15 = 2 * A + eta + zeta
        zeta15 = 2 * A + eta
        potentialBgChars.append(np.array([[A15, B15, C15], [xi15, eta15, zeta15]], dtype=np.float64))

    # test condition 16
    A16 = A
    B16 = A + C - eta
    C16 = B + C - xi
    if np.isclose(A16 + B16 + C16, sumVal, rtol=epsRel, atol=epsAbs, equal_nan=False):
        xi16 = 2 * C - xi - eta + zeta
        eta16 = -eta + zeta
        zeta16 = 2 * A - eta
        potentialBgChars.append(np.array([[A16, B16, C16], [xi16, eta16, zeta16]], dtype=np.float64))

    # test condition 17
    A17 = B
    B17 = C
    C17 = A + B + zeta
    if np.isclose(A17 + B17 + C17, sumVal, rtol=epsRel, atol=epsAbs, equal_nan=False):
        xi17 = xi + eta
        eta17 = 2 * B + zeta
        zeta17 = xi
        potentialBgChars.append(np.array([[A17, B17, C17], [xi17, eta17, zeta17]], dtype=np.float64))

    # test condition 18
    A18 = B
    B18 = C
    C18 = A + B - zeta
    if np.isclose(A18 + B18 + C18, sumVal, rtol=epsRel, atol=epsAbs, equal_nan=False):
        xi18 = -xi + eta
        eta18 = -2 * B + zeta
        zeta18 = xi
        potentialBgChars.append(np.array([[A18, B18, C18], [xi18, eta18, zeta18]], dtype=np.float64))

    # test condition 19
    A19 = B
    B19 = C
    C19 = A + C - eta
    if np.isclose(A19 + B19 + C19, sumVal, rtol=epsRel, atol=epsAbs, equal_nan=False):
        xi19 = -2 * C + eta
        eta19 = -xi + zeta
        zeta19 = xi
        potentialBgChars.append(np.array([[A19, B19, C19], [xi19, eta19, zeta19]], dtype=np.float64))

    # test condition 20
    A20 = B
    B20 = A + B + zeta
    C20 = sumVal + sigma
    if np.isclose(A20 + B20 + C20, sumVal, rtol=epsRel, atol=epsAbs, equal_nan=False):
        xi20 = 2 * A + 2 * B + sigma + zeta
        eta20 = 2 * B + xi + zeta
        zeta20 = 2 * B + zeta
        potentialBgChars.append(np.array([[A20, B20, C20], [xi20, eta20, zeta20]], dtype=np.float64))

    # test condition 21
    A21 = B
    B21 = A + B - zeta
    C21 = A + C - eta
    if np.isclose(A21 + B21 + C21, sumVal, rtol=epsRel, atol=epsAbs, equal_nan=False):
        xi21 = 2 * A + xi - eta - zeta
        eta21 = -xi + zeta
        zeta21 = -2 * B + zeta
        potentialBgChars.append(np.array([[A21, B21, C21], [xi21, eta21, zeta21]], dtype=np.float64))

    # test condition 22
    A22 = B
    B22 = A + C - eta
    C22 = B + C - xi
    if np.isclose(A22 + B22 + C22, sumVal, rtol=epsRel, atol=epsAbs, equal_nan=False):
        xi22 = 2 * C - xi - eta + zeta
        eta22 = 2 * B - xi
        zeta22 = -xi + zeta
        potentialBgChars.append(np.array([[A22, B22, C22], [xi22, eta22, zeta22]], dtype=np.float64))

    # test condition 23
    A23 = C
    B23 = A + B - zeta
    C23 = A + C - eta
    if np.isclose(A23 + B23 + C23, sumVal, rtol=epsRel, atol=epsAbs, equal_nan=False):
        xi23 = 2 * A + xi - eta - zeta
        eta23 = -2 * C + eta
        zeta23 = -xi + eta
        potentialBgChars.append(np.array([[A23, B23, C23], [xi23, eta23, zeta23]], dtype=np.float64))

    # test condition 24
    A24 = C
    B24 = A + B - zeta
    C24 = B + C - xi
    if np.isclose(A24 + B24 + C24, sumVal, rtol=epsRel, atol=epsAbs, equal_nan=False):
        xi24 = -2 * B + xi - eta + zeta
        eta24 = -2 * C + xi
        zeta24 = -xi + eta
        potentialBgChars.append(np.array([[A24, B24, C24], [xi24, eta24, zeta24]], dtype=np.float64))

    # test condition 25
    A25 = C
    B25 = A + C + eta
    C25 = sumVal + sigma
    if np.isclose(A25 + B25 + C25, sumVal, rtol=epsRel, atol=epsAbs, equal_nan=False):
        xi25 = 2 * A + 2 * C + sigma + eta
        eta25 = 2 * C + xi + eta
        zeta25 = 2 * C + eta
        potentialBgChars.append(np.array([[A25, B25, C25], [xi25, eta25, zeta25]], dtype=np.float64))

    return potentialBgChars


def checkMatMinDist(mat, matList):
    """

        :param mat: a matrix
        :param matList: a list of matrices
        :return: min distance between mat and matrices in matList
    """
    distList = [np.linalg.norm(mat - elem, ord=2) for elem in matList]
    return np.min(distList)


def allNormalizedBuergerCharacteristics(S):
    """

    :param S: returned from Algorithm B
    :return: all normalized Buerger characteristics belonging to the same Niggli cell
    """
    eps = 1e-6
    potentialBgChars = potentialBuergerCharacteristics(S)  # TODO: only unique matrices should remain!
    normalizedBgCharsAll = [AlgorithmN(STmp) for STmp in potentialBgChars]
    normalizedBgCharsReturned = []
    normalizedBgCharsReturned.append(normalizedBgCharsAll[0])
    for elem in normalizedBgCharsAll:
        if checkMatMinDist(elem, normalizedBgCharsReturned) <= eps:
            continue
        else:
            normalizedBgCharsReturned.append(elem)

    return normalizedBgCharsReturned

    # def checkMatMinDist(mat,matList):
    #     """
    #
    #     :param mat: a matrix
    #     :param matList: a list of matrices
    #     :return: min distance between mat and matrices in matList
    #     """
    #     distList=[np.linalg.norm(mat-elem,ord=2) for elem in matList]
    #     return np.min(distList)


def selectLen3(normalizedBgChars):
    """

    :param normalizedBgChars: normalized Buerger characteristics belonging to the same Niggli cell, there are 3 such matrices
    :return: the category number of the Buerger characteristics
    """
    epsRel = 1e-8
    epsAbs = 1e-6
    # if len(normalizedBgChars)!=3:
    #     raise ValueError("Invalid input: out of range.")

    S0, S1, S2 = normalizedBgChars
    [[A, B, C], [_, _, _]] = S0
    # case A=B<=C
    if np.isclose(A, B, rtol=epsRel, atol=epsAbs, equal_nan=False):
        return 3
    # case A<B<C
    if A < B and B < C:
        return 11

    # case A<B=C
    # if xi is always <0
    if S0[1, 0] < 0 and S1[1, 0] < 0 and S2[1, 0] < 0:
        return 19

    # if there exists eta=0
    if np.isclose(S0[1, 1], 0, rtol=epsRel, atol=epsAbs, equal_nan=False) \
            or np.isclose(S1[1, 1], 0, rtol=epsRel, atol=epsAbs, equal_nan=False) \
            or np.isclose(S2[1, 1], 0, rtol=epsRel, atol=epsAbs, equal_nan=False):
        return 11

    # if abs(xi) is always B
    if np.isclose(np.abs(S0[1, 0]), B, rtol=epsRel, atol=epsAbs, equal_nan=False) \
            and np.isclose(np.abs(S1[1, 0]), B, rtol=epsRel, atol=epsAbs, equal_nan=False) \
            and np.isclose(np.abs(S2[1, 0]), B, rtol=epsRel, atol=epsAbs, equal_nan=False):
        return 18

    # otherwise the category is 13
    return 13


def NiggliOfLen3(normalizedBgChars):
    """

    :param normalizedBgChars: normalized Buerger characteristics belonging to the same Niggli cell, there are 3 such matrices
    :return: S matrix of Niggli cell
    """

    category = selectLen3(normalizedBgChars)
    S0, S1, S2 = normalizedBgChars
    if category == 3:
        # find the matrix with the largest value of xi
        indOfXi = np.argsort([S0[1, 0], S1[1, 0], S2[1, 0]])
        return normalizedBgChars[indOfXi[-1]]
    if category == 11:
        if S0[1, 0] > 0:
            return S0
        elif S1[1, 0] > 0:
            return S1
        else:
            return S2

    if category == 13:
        if S0[1, 1] > 0:
            return S0
        elif S1[1, 1] > 0:
            return S1
        else:
            return S2
    if category == 18:
        indOfEta = np.argsort([S0[1, 1], S1[1, 1], S2[1, 1]])
        return normalizedBgChars[indOfEta[-1]]

    # remaining case: category=19
    if category == 19:
        indOfEta = np.argsort([S0[1, 1], S1[1, 1], S2[1, 1]])
        return normalizedBgChars[indOfEta[0]]

    # if not returned until now, then something is wrong
    raise RuntimeError("Something is wrong for 3-Buerger-cell case.")


def selectLen2(normalizedBgChars):
    """

    :param normalizedBgChars: normalized Buerger characteristics belonging to the same Niggli cell, there are 2 such matrices
    :return: the representative category number of the Buerger characteristics
    """
    S0, S1 = normalizedBgChars

    # all categories except for category 21
    if S0[1, 0] > S1[1, 0] or S1[1, 0] > S0[1, 0]:
        return 1
    else:
        return 21


def NiggliOfLen2(normalizedBgChars):
    """

    :param normalizedBgChars: normalized Buerger characteristics belonging to the same Niggli cell, there are 2 such matrices
    :return: S matrix of Niggli cell
    """
    category = selectLen2(normalizedBgChars)
    S0, S1 = normalizedBgChars
    if category == 1:
        if S0[1, 0] > 0:
            return S0
        else:
            return S1
    if category == 21:
        if S0[1, 1] > S1[1, 1]:
            return S0
        else:
            return S1


def Buerger2Niggli(normalizedBgChars):
    """

    :param normalizedBgChars: normalized Buerger characteristics belonging to the same Niggli cell
    :return: S matrix of Niggli cell
    """

    # if there is only 1 Buerger cell, then this is the Niggli cell

    if len(normalizedBgChars) == 1:
        return normalizedBgChars[0]

    # if there are 5 Buerger cells, they belong to category 20
    if len(normalizedBgChars) == 5:
        candidates = [STmp for STmp in normalizedBgChars if STmp[1, 0] > 0]
        S0 = candidates[0]
        S1 = candidates[1]
        if S0[1, 1] > S1[1, 1]:
            return S0
        else:
            return S1

    # if there are 4 Buerger cells, they belong to category 23
    if len(normalizedBgChars) == 4:
        candidates = [STmp for STmp in normalizedBgChars if STmp[1, 0] > 0]
        S0 = candidates[0]
        S1 = candidates[1]
        if S0[1, 1] > S1[1, 1]:
            return S0
        else:
            return S1

    # if there are 3 Buerger cells
    if len(normalizedBgChars) == 3:
        return NiggliOfLen3(normalizedBgChars)
    # if there are 2 Buerger cells
    if len(normalizedBgChars) == 2:
        return NiggliOfLen2(normalizedBgChars)

    raise RuntimeError("Wrong number of matrices.")


def checkNiggliCell(S):
    """

    :param S: output of Buerger2Niggli(), supposed to be a Niggli cell
    :return: check if the execution of Buerger2Niggli() is correct.
    """

    [[A, B, C], [xi, eta, zeta]] = S
    epsRel = 1e-8
    epsAbs = 1e-6
    # type-I cell
    if xi > 0 and eta > 0 and zeta > 0:
        # main conditions
        if not (
                (A < B or np.isclose(A, B, rtol=epsRel, atol=epsAbs, equal_nan=False))
                and (B < C or np.isclose(B, C, rtol=epsRel, atol=epsAbs, equal_nan=False))
        ):
            raise ValueError("Main condition A<=B<=C not satisfied in type-I.")
        if not (
                (np.abs(xi) < B or np.isclose(np.abs(xi), B, rtol=epsRel, atol=epsAbs, equal_nan=False))
                and (np.abs(eta) < A or np.isclose(np.abs(eta), A, rtol=epsRel, atol=epsAbs, equal_nan=False))
                and (np.abs(zeta) < A or np.isclose(np.abs(zeta), A, rtol=epsRel, atol=epsAbs, equal_nan=False))
        ):
            raise ValueError("Main condition |xi|<=B, |eta|<=A, |zeta|<=A not satisfied in type-I.")
        # special conditions
        # 1
        if np.isclose(A, B, rtol=epsRel, atol=epsAbs, equal_nan=False):
            if not (xi < eta or np.isclose(xi, eta, rtol=epsRel, atol=epsAbs, equal_nan=False)):
                raise ValueError("Special condition xi<=eta in type-I cell not satisfied.")
        # 2
        if np.isclose(B, C, rtol=epsRel, atol=epsAbs, equal_nan=False):
            if not (eta < zeta or np.isclose(eta, zeta, rtol=epsRel, atol=epsAbs, equal_nan=False)):
                raise ValueError("Special condition eta<=zeta in type-I cell not satisfied.")
        # 3
        if np.isclose(xi, B, rtol=epsRel, atol=epsAbs, equal_nan=False):
            if not (1 / 2 * zeta < eta or np.isclose(1 / 2 * zeta, eta, rtol=epsRel, atol=epsAbs, equal_nan=False)):
                raise ValueError("Special condition 1/2 zeta<=eta in type-I cell not satisfied.")

        # 4
        if np.isclose(eta, A, rtol=epsRel, atol=epsAbs, equal_nan=False):
            if not (1 / 2 * zeta < xi or np.isclose(1 / 2 * zeta, xi, rtol=epsRel, atol=epsAbs, equal_nan=False)):
                raise ValueError("Special condition 1/2 zeta<=xi in type-I cell not satisfied.")

        # 5
        if np.isclose(zeta, A, rtol=epsRel, atol=epsAbs, equal_nan=False):
            if not (1 / 2 * eta < xi or np.isclose(1 / 2 * eta, xi, rtol=epsRel, atol=epsAbs, equal_nan=False)):
                raise ValueError("Special condition 1/2 eta<=xi in type-I cell not satisfied.")

    # type-II cell
    if ((xi < 0 or np.isclose(xi, 0, rtol=epsRel, atol=epsAbs, equal_nan=False))
            and (eta < 0 or np.isclose(eta, 0, rtol=epsRel, atol=epsAbs, equal_nan=False))
            and (zeta < 0 or np.isclose(zeta, 0, rtol=epsRel, atol=epsAbs, equal_nan=False))):
        # main conditions
        if not (
                (A < B or np.isclose(A, B, rtol=epsRel, atol=epsAbs, equal_nan=False))
                and (B < C or np.isclose(B, C, rtol=epsRel, atol=epsAbs, equal_nan=False))
        ):
            raise ValueError("Main condition A<=B<=C not satisfied in type II.")
        if not (
                (np.abs(xi) < B or np.isclose(np.abs(xi), B, rtol=epsRel, atol=epsAbs, equal_nan=False))
                and (np.abs(eta) < A or np.isclose(np.abs(eta), A, rtol=epsRel, atol=epsAbs, equal_nan=False))
                and (np.abs(zeta) < A or np.isclose(np.abs(zeta), A, rtol=epsRel, atol=epsAbs, equal_nan=False))
        ):
            raise ValueError("Main condition |xi|<=B, |eta|<=A, |zeta|<=A not satisfied in type II.")
        if not (
                np.abs(xi) + np.abs(eta) + np.abs(zeta) < A + B
                or np.isclose(np.abs(xi) + np.abs(eta) + np.abs(zeta), A + B, rtol=epsRel, atol=epsAbs, equal_nan=False)
        ):
            raise ValueError("Main condition |xi|+|eta|+|zeta|<=A+B not satisfied in type II.")
        # special conditions
        # 1
        if np.isclose(A, B, rtol=epsRel, atol=epsAbs, equal_nan=False):
            if not (
                    np.abs(xi) < np.abs(eta)
                    or np.isclose(np.abs(xi), np.abs(eta), rtol=epsRel, atol=epsAbs, equal_nan=False)
            ):
                raise ValueError("Special condition |xi|<=|eta| in type-II cell not satisfied.")

        # 2
        if np.isclose(B, C, rtol=epsRel, atol=epsAbs, equal_nan=False):
            if not (np.abs(eta) < np.abs(zeta)
                    or np.isclose(np.abs(eta), np.abs(zeta), rtol=epsRel, atol=epsAbs, equal_nan=False)
            ):
                raise ValueError("Special condition |eta|<=|zeta| in type-II cell not satisfied.")
        # 3
        if np.isclose(np.abs(xi), B, rtol=epsRel, atol=epsAbs, equal_nan=False):
            if not np.isclose(zeta, 0, rtol=epsRel, atol=epsAbs, equal_nan=False):
                raise ValueError("Special condition zeta=0 in type-II cell not satisfied.")
        # 4
        if np.isclose(np.abs(eta), A, rtol=epsRel, atol=epsAbs, equal_nan=False):
            if not np.isclose(zeta, 0, rtol=epsRel, atol=epsAbs, equal_nan=False):
                raise ValueError("Special condition zeta=0 in type-II cell not satisfied.")
        # 5
        if np.isclose(np.abs(zeta), A, rtol=epsRel, atol=epsAbs, equal_nan=False):
            if not np.isclose(eta, 0, rtol=epsRel, atol=epsAbs, equal_nan=False):
                raise ValueError("Special condition eta=0 in type-II cell not satisfied.")

        # 6
        if np.isclose(np.abs(xi) + np.abs(eta) + np.abs(zeta), A + B, rtol=epsRel, atol=epsAbs, equal_nan=False):
            if not (
                    A < np.abs(eta) + 1 / 2 * np.abs(zeta)
                    or np.isclose(A, np.abs(eta) + 1 / 2 * np.abs(zeta), rtol=epsRel, atol=epsAbs, equal_nan=False)
            ):
                raise ValueError("Special condition A<=|eta|+1/2|zeta| in type-II cell not satisfied.")


def NiggliType(S):
    """

    :param S: S matrix of a Niggli cell
    :return: type 1 or 2
    """
    epsRel = 1e-8
    epsAbs = 1e-6
    [[_, _, _], [xi, eta, zeta]] = S
    if xi > 0 and eta > 0 and zeta > 0:
        return 1
    elif (xi < 0 or np.isclose(xi, 0, rtol=epsRel, atol=epsAbs, equal_nan=False)) \
            and (eta < 0 or np.isclose(eta, 0, rtol=epsRel, atol=epsAbs, equal_nan=False)) \
            and (zeta < 0 or np.isclose(zeta, 0, rtol=epsRel, atol=epsAbs, equal_nan=False)):
        return 2
    else:
        raise ValueError("Invalid Niggli matrix.")


def Niggli2ConventionalRow(S):
    """

    :param S: S matrix of Niggli cell returned from Buerger2Niggli()
    :return: Row in the table 9.2.5.1 in  International Tables for Crystallography Volume A_ Space-group symmetry-Springer Netherlands (2002)
    """
    [[A, B, C], [xi, eta, zeta]] = S
    D = 1 / 2 * xi
    E = 1 / 2 * eta
    F = 1 / 2 * zeta
    epsRel = 1e-8
    epsAbs = 1e-6

    # A=B=C
    if np.isclose(A, B, rtol=epsRel, atol=epsAbs, equal_nan=False) \
            and np.isclose(B, C, rtol=epsRel, atol=epsAbs, equal_nan=False):
        # 1
        if NiggliType(S) == 1 \
                and np.isclose(D, A / 2, rtol=epsRel, atol=epsAbs, equal_nan=False) \
                and np.isclose(E, A / 2, rtol=epsRel, atol=epsAbs, equal_nan=False) \
                and np.isclose(F, A / 2, rtol=epsRel, atol=epsAbs, equal_nan=False):
            return {"No": 1, "Lattice symmetry": "Cubic", "Bravais type": "cF",
                    "Transformation matrix": np.array([[1, -1, 1], [1, 1, -1], [-1, 1, 1]], dtype=np.float64)}
        # 2
        if NiggliType(S) == 1 \
                and np.isclose(D, D, rtol=epsRel, atol=epsAbs, equal_nan=False) \
                and np.isclose(E, D, rtol=epsRel, atol=epsAbs, equal_nan=False) \
                and np.isclose(F, D, rtol=epsRel, atol=epsAbs, equal_nan=False):
            return {"No": 2, "Lattice symmetry": "Rhombohedral", "Bravais type": "hR",
                    "Transformation matrix": np.array([[1, -1, 0], [-1, 0, 1], [-1, -1, -1]], dtype=np.float64)}

        # 3
        if NiggliType(S) == 2 \
                and np.isclose(D, 0, rtol=epsRel, atol=epsAbs, equal_nan=False) \
                and np.isclose(E, 0, rtol=epsRel, atol=epsAbs, equal_nan=False) \
                and np.isclose(F, 0, rtol=epsRel, atol=epsAbs, equal_nan=False):
            return {"No": 3, "Lattice symmetry": "Cubic", "Bravais type": "cP",
                    "Transformation matrix": np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=np.float64)}

        # 5
        if NiggliType(S) == 2 \
                and np.isclose(D, -A / 3, rtol=epsRel, atol=epsAbs, equal_nan=False) \
                and np.isclose(E, -A / 3, rtol=epsRel, atol=epsAbs, equal_nan=False) \
                and np.isclose(F, -A / 3, rtol=epsRel, atol=epsAbs, equal_nan=False):
            return {"No": 5, "Lattice symmetry": "Cubic", "Bravais type": "cI",
                    "Transformation matrix": np.array([[1, 0, 1], [1, 1, 0], [0, 1, 1]], dtype=np.float64)}

        # 4
        if NiggliType(S) == 2 \
                and np.isclose(D, D, rtol=epsRel, atol=epsAbs, equal_nan=False) \
                and np.isclose(E, D, rtol=epsRel, atol=epsAbs, equal_nan=False) \
                and np.isclose(F, D, rtol=epsRel, atol=epsAbs, equal_nan=False):
            return {"No": 4, "Lattice symmetry": "Rhombohedral", "Bravais type": "hR",
                    "Transformation matrix": np.array([[1, -1, 0], [-1, 0, 1], [-1, -1, -1]], dtype=np.float64)}

        # 6
        if NiggliType(S) == 2 \
                and np.isclose(D, D, rtol=epsRel, atol=epsAbs, equal_nan=False) \
                and np.isclose(E, D, rtol=epsRel, atol=epsAbs, equal_nan=False) \
                and np.isclose(F, F, rtol=epsRel, atol=epsAbs, equal_nan=False) \
                and np.isclose(2 * np.abs(D + E + F), A + B, rtol=epsRel, atol=epsAbs, equal_nan=False):
            return {"No": 6, "Lattice symmetry": "Tetragonal", "Bravais type": "tI",
                    "Transformation matrix": np.array([[0, 1, 1], [1, 0, 1], [1, 1, 0]], dtype=np.float64)}

        # 7
        if NiggliType(S) == 2 \
                and np.isclose(D, D, rtol=epsRel, atol=epsAbs, equal_nan=False) \
                and np.isclose(E, E, rtol=epsRel, atol=epsAbs, equal_nan=False) \
                and np.isclose(F, E, rtol=epsRel, atol=epsAbs, equal_nan=False) \
                and np.isclose(2 * np.abs(D + E + F), A + B, rtol=epsRel, atol=epsAbs, equal_nan=False):
            return {"No": 7, "Lattice symmetry": "Tetragonal", "Bravais type": "tI",
                    "Transformation matrix": np.array([[1, 0, 1], [1, 1, 0], [0, 1, 1]], dtype=np.float64)}

        # 8
        if NiggliType(S) == 2 \
                and np.isclose(D, D, rtol=epsRel, atol=epsAbs, equal_nan=False) \
                and np.isclose(E, E, rtol=epsRel, atol=epsAbs, equal_nan=False) \
                and np.isclose(F, F, rtol=epsRel, atol=epsAbs, equal_nan=False) \
                and np.isclose(2 * np.abs(D + E + F), A + B, rtol=epsRel, atol=epsAbs, equal_nan=False):
            return {"No": 8, "Lattice symmetry": "Orthorhombic", "Bravais type": "oI",
                    "Transformation matrix": np.array([[-1, -1, 0], [-1, 0, -1], [0, -1, -1]], dtype=np.float64)}

    # A=B, no conditions on C
    if np.isclose(A, B, rtol=epsRel, atol=epsAbs, equal_nan=False):
        # 9
        if NiggliType(S) == 1 \
                and np.isclose(D, A / 2, rtol=epsRel, atol=epsAbs, equal_nan=False) \
                and np.isclose(E, A / 2, rtol=epsRel, atol=epsAbs, equal_nan=False) \
                and np.isclose(F, A / 2, rtol=epsRel, atol=epsAbs, equal_nan=False):
            return {"No": 9, "Lattice symmetry": "Rhombohedral", "Bravais type": "hR",
                    "Transformation matrix": np.array([[1, 0, 0], [-1, 1, 0], [-1, -1, 3]], dtype=np.float64)}

        # 10
        if NiggliType(S) == 1 \
                and np.isclose(D, D, rtol=epsRel, atol=epsAbs, equal_nan=False) \
                and np.isclose(E, D, rtol=epsRel, atol=epsAbs, equal_nan=False) \
                and np.isclose(F, F, rtol=epsRel, atol=epsAbs, equal_nan=False):
            return {"No": 10, "Lattice symmetry": "Monoclinic", "Bravais type": "mC",
                    "Transformation matrix": np.array([[1, 1, 0], [1, -1, 0], [0, 0, -1]], dtype=np.float64)}

        # 11
        if NiggliType(S) == 2 \
                and np.isclose(D, 0, rtol=epsRel, atol=epsAbs, equal_nan=False) \
                and np.isclose(E, 0, rtol=epsRel, atol=epsAbs, equal_nan=False) \
                and np.isclose(F, 0, rtol=epsRel, atol=epsAbs, equal_nan=False):
            return {"No": 11, "Lattice symmetry": "Tetragonal", "Bravais type": "tP",
                    "Transformation matrix": np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=np.float64)}

        # 12
        if NiggliType(S) == 2 \
                and np.isclose(D, 0, rtol=epsRel, atol=epsAbs, equal_nan=False) \
                and np.isclose(E, 0, rtol=epsRel, atol=epsAbs, equal_nan=False) \
                and np.isclose(F, -A / 2, rtol=epsRel, atol=epsAbs, equal_nan=False):
            return {"No": 12, "Lattice symmetry": "Hexagonal", "Bravais type": "hP",
                    "Transformation matrix": np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=np.float64)}

        # 13
        if NiggliType(S) == 2 \
                and np.isclose(D, 0, rtol=epsRel, atol=epsAbs, equal_nan=False) \
                and np.isclose(E, 0, rtol=epsRel, atol=epsAbs, equal_nan=False) \
                and np.isclose(F, F, rtol=epsRel, atol=epsAbs, equal_nan=False):
            return {"No": 13, "Lattice symmetry": "Orthorhombic", "Bravais type": "oC",
                    "Transformation matrix": np.array([[1, 1, 0], [-1, 1, 0], [0, 0, 1]], dtype=np.float64)}

        # 15
        if NiggliType(S) == 2 \
                and np.isclose(D, -A / 2, rtol=epsRel, atol=epsAbs, equal_nan=False) \
                and np.isclose(E, -A / 2, rtol=epsRel, atol=epsAbs, equal_nan=False) \
                and np.isclose(F, 0, rtol=epsRel, atol=epsAbs, equal_nan=False):
            return {"No": 15, "Lattice symmetry": "Tetragonal", "Bravais type": "tI",
                    "Transformation matrix": np.array([[1, 0, 0], [0, 1, 0], [1, 1, 2]], dtype=np.float64)}

        # 16
        if NiggliType(S) == 2 \
                and np.isclose(D, D, rtol=epsRel, atol=epsAbs, equal_nan=False) \
                and np.isclose(E, D, rtol=epsRel, atol=epsAbs, equal_nan=False) \
                and np.isclose(F, F, rtol=epsRel, atol=epsAbs, equal_nan=False) \
                and np.isclose(2 * np.abs(D + E + F), A + B, rtol=epsRel, atol=epsAbs, equal_nan=False):
            return {"No": 16, "Lattice symmetry": "Orthorhombic", "Bravais type": "oF",
                    "Transformation matrix": np.array([[-1, -1, 0], [1, -1, 0], [1, 1, 2]], dtype=np.float64)}

        # 14
        if NiggliType(S) == 2 \
                and np.isclose(D, D, rtol=epsRel, atol=epsAbs, equal_nan=False) \
                and np.isclose(E, D, rtol=epsRel, atol=epsAbs, equal_nan=False) \
                and np.isclose(F, F, rtol=epsRel, atol=epsAbs, equal_nan=False):
            return {"No": 14, "Lattice symmetry": "Monoclinic", "Bravais type": "mC",
                    "Transformation matrix": np.array([[1, 1, 0], [-1, 1, 0], [0, 0, 1]], dtype=np.float64)}

        # 17
        if NiggliType(S) == 2 \
                and np.isclose(D, D, rtol=epsRel, atol=epsAbs, equal_nan=False) \
                and np.isclose(E, E, rtol=epsRel, atol=epsAbs, equal_nan=False) \
                and np.isclose(F, F, rtol=epsRel, atol=epsAbs, equal_nan=False) \
                and np.isclose(2 * np.abs(D + E + F), A + B, rtol=epsRel, atol=epsAbs, equal_nan=False):
            return {"No": 17, "Lattice symmetry": "Monoclinic", "Bravais type": "mC",
                    "Transformation matrix": np.array([[1, -1, 0], [1, 1, 0], [-1, 0, -1]], dtype=np.float64)}

    # B=C, no conditions on A
    if np.isclose(B, C, rtol=epsRel, atol=epsAbs, equal_nan=False):
        # 18
        if NiggliType(S) == 1 \
                and np.isclose(D, A / 4, rtol=epsRel, atol=epsAbs, equal_nan=False) \
                and np.isclose(E, A / 2, rtol=epsRel, atol=epsAbs, equal_nan=False) \
                and np.isclose(F, A / 2, rtol=epsRel, atol=epsAbs, equal_nan=False):
            return {"No": 18, "Lattice symmetry": "Tetragonal", "Bravais type": "tI",
                    "Transformation matrix": np.array([[0, -1, 1], [1, -1, -1], [1, 0, 0]], dtype=np.float64)}
        # 19
        if NiggliType(S) == 1 \
                and np.isclose(D, D, rtol=epsRel, atol=epsAbs, equal_nan=False) \
                and np.isclose(E, A / 2, rtol=epsRel, atol=epsAbs, equal_nan=False) \
                and np.isclose(F, A / 2, rtol=epsRel, atol=epsAbs, equal_nan=False):
            return {"No": 19, "Lattice symmetry": "Orthorhombic", "Bravais type": "oI",
                    "Transformation matrix": np.array([[-1, 0, 0], [0, -1, 1], [-1, 1, 1]], dtype=np.float64)}

        # 20
        if NiggliType(S) == 1 \
                and np.isclose(D, D, rtol=epsRel, atol=epsAbs, equal_nan=False) \
                and np.isclose(E, E, rtol=epsRel, atol=epsAbs, equal_nan=False) \
                and np.isclose(F, E, rtol=epsRel, atol=epsAbs, equal_nan=False):
            return {"No": 20, "Lattice symmetry": "Monoclinic", "Bravais type": "mC",
                    "Transformation matrix": np.array([[0, 1, 1], [0, 1, -1], [-1, 0, 0]], dtype=np.float64)}
        # 21
        if NiggliType(S) == 2 \
                and np.isclose(D, 0, rtol=epsRel, atol=epsAbs, equal_nan=False) \
                and np.isclose(E, 0, rtol=epsRel, atol=epsAbs, equal_nan=False) \
                and np.isclose(F, 0, rtol=epsRel, atol=epsAbs, equal_nan=False):
            return {"No": 21, "Lattice symmetry": "Tetragonal", "Bravais type": "tP",
                    "Transformation matrix": np.array([[0, 1, 0], [0, 0, 1], [1, 0, 0]], dtype=np.float64)}

        # 22
        if NiggliType(S) == 2 \
                and np.isclose(D, -B / 2, rtol=epsRel, atol=epsAbs, equal_nan=False) \
                and np.isclose(E, 0, rtol=epsRel, atol=epsAbs, equal_nan=False) \
                and np.isclose(F, 0, rtol=epsRel, atol=epsAbs, equal_nan=False):
            return {"No": 22, "Lattice symmetry": "Hexagonal", "Bravais type": "hP",
                    "Transformation matrix": np.array([[0, 1, 0], [0, 0, 1], [1, 0, 0]], dtype=np.float64)}
        # 23
        if NiggliType(S) == 2 \
                and np.isclose(D, D, rtol=epsRel, atol=epsAbs, equal_nan=False) \
                and np.isclose(E, 0, rtol=epsRel, atol=epsAbs, equal_nan=False) \
                and np.isclose(F, 0, rtol=epsRel, atol=epsAbs, equal_nan=False):
            return {"No": 23, "Lattice symmetry": "Orthorhombic", "Bravais type": "oC",
                    "Transformation matrix": np.array([[0, 1, 1], [0, -1, 1], [1, 0, 0]], dtype=np.float64)}

        # 24
        if NiggliType(S) == 2 \
                and np.isclose(D, D, rtol=epsRel, atol=epsAbs, equal_nan=False) \
                and np.isclose(E, -A / 3, rtol=epsRel, atol=epsAbs, equal_nan=False) \
                and np.isclose(F, -A / 3, rtol=epsRel, atol=epsAbs, equal_nan=False) \
                and np.isclose(2 * np.abs(D + E + F), A + B, rtol=epsRel, atol=epsAbs, equal_nan=False):
            return {"No": 24, "Lattice symmetry": "Rhombohedral", "Bravais type": "hR",
                    "Transformation matrix": np.array([[1, 2, 1], [0, -1, 1], [1, 0, 0]], dtype=np.float64)}

        # 25
        if NiggliType(S) == 2 \
                and np.isclose(D, D, rtol=epsRel, atol=epsAbs, equal_nan=False) \
                and np.isclose(E, E, rtol=epsRel, atol=epsAbs, equal_nan=False) \
                and np.isclose(F, E, rtol=epsRel, atol=epsAbs, equal_nan=False):
            return {"No": 25, "Lattice symmetry": "Monoclinic", "Bravais type": "mC",
                    "Transformation matrix": np.array([[0, 1, 1], [0, -1, 1], [1, 0, 0]], dtype=np.float64)}

    # no conditions on A, B, C

    # 26
    if NiggliType(S) == 1 \
            and np.isclose(D, A / 4, rtol=epsRel, atol=epsAbs, equal_nan=False) \
            and np.isclose(E, A / 2, rtol=epsRel, atol=epsAbs, equal_nan=False) \
            and np.isclose(F, A / 2, rtol=epsRel, atol=epsAbs, equal_nan=False):
        return {"No": 26, "Lattice symmetry": "Orthorhombic", "Bravais type": "oF",
                "Transformation matrix": np.array([[1, 0, 0], [-1, 2, 0], [-1, 0, 2]], dtype=np.float64)}

    # 27
    if NiggliType(S) == 1 \
            and np.isclose(D, D, rtol=epsRel, atol=epsAbs, equal_nan=False) \
            and np.isclose(E, A / 2, rtol=epsRel, atol=epsAbs, equal_nan=False) \
            and np.isclose(F, A / 2, rtol=epsRel, atol=epsAbs, equal_nan=False):
        return {"No": 27, "Lattice symmetry": "Monoclinic", "Bravais type": "mC",
                "Transformation matrix": np.array([[-1, 2, 0], [-1, 0, 0], [0, -1, 1]], dtype=np.float64)}

    # 28
    if NiggliType(S) == 1 \
            and np.isclose(D, D, rtol=epsRel, atol=epsAbs, equal_nan=False) \
            and np.isclose(E, A / 2, rtol=epsRel, atol=epsAbs, equal_nan=False) \
            and np.isclose(F, 2 * D, rtol=epsRel, atol=epsAbs, equal_nan=False):
        return {"No": 28, "Lattice symmetry": "Monoclinic", "Bravais type": "mC",
                "Transformation matrix": np.array([[-1, 0, 0], [-1, 0, 2], [0, 1, 0]], dtype=np.float64)}

    # 29
    if NiggliType(S) == 1 \
            and np.isclose(D, D, rtol=epsRel, atol=epsAbs, equal_nan=False) \
            and np.isclose(E, 2 * D, rtol=epsRel, atol=epsAbs, equal_nan=False) \
            and np.isclose(F, A / 2, rtol=epsRel, atol=epsAbs, equal_nan=False):
        return {"No": 29, "Lattice symmetry": "Monoclinic", "Bravais type": "mC",
                "Transformation matrix": np.array([[1, 0, 0], [1, -2, 0], [0, 0, -1]], dtype=np.float64)}

    # 30
    if NiggliType(S) == 1 \
            and np.isclose(D, B / 2, rtol=epsRel, atol=epsAbs, equal_nan=False) \
            and np.isclose(E, E, rtol=epsRel, atol=epsAbs, equal_nan=False) \
            and np.isclose(F, 2 * E, rtol=epsRel, atol=epsAbs, equal_nan=False):
        return {"No": 30, "Lattice symmetry": "Monoclinic", "Bravais type": "mC",
                "Transformation matrix": np.array([[0, 1, 0], [0, 1, -2], [-1, 0, 0]], dtype=np.float64)}

    # 31
    if NiggliType(S) == 1 \
            and np.isclose(D, D, rtol=epsRel, atol=epsAbs, equal_nan=False) \
            and np.isclose(E, E, rtol=epsRel, atol=epsAbs, equal_nan=False) \
            and np.isclose(F, F, rtol=epsRel, atol=epsAbs, equal_nan=False):
        return {"No": 31, "Lattice symmetry": "Triclinic", "Bravais type": "aP",
                "Transformation matrix": np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=np.float64)}

    # 32
    if NiggliType(S) == 2 \
            and np.isclose(D, 0, rtol=epsRel, atol=epsAbs, equal_nan=False) \
            and np.isclose(E, 0, rtol=epsRel, atol=epsAbs, equal_nan=False) \
            and np.isclose(F, 0, rtol=epsRel, atol=epsAbs, equal_nan=False):
        return {"No": 32, "Lattice symmetry": "Orthorhombic", "Bravais type": "oP",
                "Transformation matrix": np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=np.float64)}

    # 40
    if NiggliType(S) == 2 \
            and np.isclose(D, -B / 2, rtol=epsRel, atol=epsAbs, equal_nan=False) \
            and np.isclose(E, 0, rtol=epsRel, atol=epsAbs, equal_nan=False) \
            and np.isclose(F, 0, rtol=epsRel, atol=epsAbs, equal_nan=False):
        return {"No": 40, "Lattice symmetry": "Orthorhombic", "Bravais type": "oC",
                "Transformation matrix": np.array([[0, -1, 0], [0, 1, 2], [-1, 0, 0]], dtype=np.float64)}

    # 35
    if NiggliType(S) == 2 \
            and np.isclose(D, D, rtol=epsRel, atol=epsAbs, equal_nan=False) \
            and np.isclose(E, 0, rtol=epsRel, atol=epsAbs, equal_nan=False) \
            and np.isclose(F, 0, rtol=epsRel, atol=epsAbs, equal_nan=False):
        return {"No": 35, "Lattice symmetry": "Monoclinic", "Bravais type": "mP",
                "Transformation matrix": np.array([[0, -1, 0], [-1, 0, 0], [0, 0, -1]], dtype=np.float64)}

    # 36
    if NiggliType(S) == 2 \
            and np.isclose(D, 0, rtol=epsRel, atol=epsAbs, equal_nan=False) \
            and np.isclose(E, -A / 2, rtol=epsRel, atol=epsAbs, equal_nan=False) \
            and np.isclose(F, 0, rtol=epsRel, atol=epsAbs, equal_nan=False):
        return {"No": 36, "Lattice symmetry": "Orthorhombic", "Bravais type": "oC",
                "Transformation matrix": np.array([[1, 0, 0], [-1, 0, -2], [0, 1, 0]], dtype=np.float64)}

    # 33
    if NiggliType(S) == 2 \
            and np.isclose(D, 0, rtol=epsRel, atol=epsAbs, equal_nan=False) \
            and np.isclose(E, E, rtol=epsRel, atol=epsAbs, equal_nan=False) \
            and np.isclose(F, 0, rtol=epsRel, atol=epsAbs, equal_nan=False):
        return {"No": 33, "Lattice symmetry": "Monoclinic", "Bravais type": "mP",
                "Transformation matrix": np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=np.float64)}

    # 38
    if NiggliType(S) == 2 \
            and np.isclose(D, 0, rtol=epsRel, atol=epsAbs, equal_nan=False) \
            and np.isclose(E, 0, rtol=epsRel, atol=epsAbs, equal_nan=False) \
            and np.isclose(F, -A / 2, rtol=epsRel, atol=epsAbs, equal_nan=False):
        return {"No": 38, "Lattice symmetry": "Orthorhombic", "Bravais type": "oC",
                "Transformation matrix": np.array([[-1, 0, 0], [1, 2, 0], [0, 0, -1]], dtype=np.float64)}

    # 34
    if NiggliType(S) == 2 \
            and np.isclose(D, 0, rtol=epsRel, atol=epsAbs, equal_nan=False) \
            and np.isclose(E, 0, rtol=epsRel, atol=epsAbs, equal_nan=False) \
            and np.isclose(F, F, rtol=epsRel, atol=epsAbs, equal_nan=False):
        return {"No": 34, "Lattice symmetry": "Monoclinic", "Bravais type": "mP",
                "Transformation matrix": np.array([[-1, 0, 0], [0, 0, -1], [0, -1, 0]], dtype=np.float64)}

    # 42
    if NiggliType(S) == 2 \
            and np.isclose(D, -B / 2, rtol=epsRel, atol=epsAbs, equal_nan=False) \
            and np.isclose(E, -A / 2, rtol=epsRel, atol=epsAbs, equal_nan=False) \
            and np.isclose(F, 0, rtol=epsRel, atol=epsAbs, equal_nan=False):
        return {"No": 42, "Lattice symmetry": "Orthorhombic", "Bravais type": "oI",
                "Transformation matrix": np.array([[-1, 0, 0], [0, -1, 0], [1, 1, 2]], dtype=np.float64)}

    # 41
    if NiggliType(S) == 2 \
            and np.isclose(D, -B / 2, rtol=epsRel, atol=epsAbs, equal_nan=False) \
            and np.isclose(E, E, rtol=epsRel, atol=epsAbs, equal_nan=False) \
            and np.isclose(F, 0, rtol=epsRel, atol=epsAbs, equal_nan=False):
        return {"No": 41, "Lattice symmetry": "Monoclinic", "Bravais type": "mC",
                "Transformation matrix": np.array([[0, -1, -2], [0, -1, 0], [-1, 0, 0]], dtype=np.float64)}

    # 37
    if NiggliType(S) == 2 \
            and np.isclose(D, D, rtol=epsRel, atol=epsAbs, equal_nan=False) \
            and np.isclose(E, -A / 2, rtol=epsRel, atol=epsAbs, equal_nan=False) \
            and np.isclose(F, 0, rtol=epsRel, atol=epsAbs, equal_nan=False):
        return {"No": 37, "Lattice symmetry": "Monoclinic", "Bravais type": "mC",
                "Transformation matrix": np.array([[1, 0, 2], [1, 0, 0], [0, 1, 0]], dtype=np.float64)}

    # 39
    if NiggliType(S) == 2 \
            and np.isclose(D, D, rtol=epsRel, atol=epsAbs, equal_nan=False) \
            and np.isclose(E, 0, rtol=epsRel, atol=epsAbs, equal_nan=False) \
            and np.isclose(F, -A / 2, rtol=epsRel, atol=epsAbs, equal_nan=False):
        return {"No": 39, "Lattice symmetry": "Monoclinic", "Bravais type": "mC",
                "Transformation matrix": np.array([[-1, -2, 0], [-1, 0, 0], [0, 0, -1]], dtype=np.float64)}

    # 43
    if NiggliType(S) == 2 \
            and np.isclose(D, D, rtol=epsRel, atol=epsAbs, equal_nan=False) \
            and np.isclose(E, E, rtol=epsRel, atol=epsAbs, equal_nan=False) \
            and np.isclose(F, F, rtol=epsRel, atol=epsAbs, equal_nan=False) \
            and np.isclose(2 * np.abs(D + E + F), A + B, rtol=epsRel, atol=epsAbs, equal_nan=False) \
            and np.isclose(np.abs(2 * D + F), B, rtol=epsRel, atol=epsAbs, equal_nan=False):
        return {"No": 43, "Lattice symmetry": "Monoclinic", "Bravais type": "mI",
                "Transformation matrix": np.array([[-1, 0, 0], [-1, -1, -2], [0, -1, 0]], dtype=np.float64)}

    # 44
    if NiggliType(S) == 2 \
            and np.isclose(D, D, rtol=epsRel, atol=epsAbs, equal_nan=False) \
            and np.isclose(E, E, rtol=epsRel, atol=epsAbs, equal_nan=False) \
            and np.isclose(F, F, rtol=epsRel, atol=epsAbs, equal_nan=False):
        return {"No": 44, "Lattice symmetry": "Triclinic", "Bravais type": "aP",
                "Transformation matrix": np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=np.float64)}

    # if none of above is satisfied
    raise ValueError("Conventional cell is not found.")


def getConventional(S):
    """

    :param S: S matrix of Niggli cell returned from Buerger2Niggli()
    :return: information concerning the conventional cell
    """
    convInfoList = Niggli2ConventionalRow(S)
    # No=convInfoList["No"]
    # crystalSystem=convInfoList["Lattice symmetry"]
    # BrvTyp=convInfoList["Bravais type"]
    M = convInfoList["Transformation matrix"]

    [[A, B, C], [xi, eta, zeta]] = S
    D = 1 / 2 * xi
    E = 1 / 2 * eta
    F = 1 / 2 * zeta
    m0 = M[0, :]
    m1 = M[1, :]
    m2 = M[2, :]
    # for definition of R2, see notes
    R2 = np.array([[A, F, E],
                   [F, B, D],
                   [E, D, C]])

    AConv = m0 @ R2 @ m0

    BConv = m1 @ R2 @ m1

    CConv = m2 @ R2 @ m2

    DConv = m1 @ R2 @ m2

    EConv = m2 @ R2 @ m0

    FConv = m0 @ R2 @ m1

    aConv = np.sqrt(AConv)
    bConv = np.sqrt(BConv)
    cConv = np.sqrt(CConv)

    alphaConv = np.arccos(DConv / (bConv * cConv))

    betaConv = np.arccos(EConv / (aConv * cConv))

    gammaConv = np.arccos(FConv / (aConv * bConv))

    convParamsAndInfo = {"a": aConv, "b": bConv, "c": cConv, "alpha": alphaConv, "beta": betaConv, "gamma": gammaConv,
                         "Crystal System": convInfoList["Lattice symmetry"],
                         "Bravais type": convInfoList["Bravais type"]}
    # the angles are in radians
    return convParamsAndInfo


def getConventionalInDegree(S):
    """

    :param S: S matrix of Niggli cell returned from Buerger2Niggli()
    :return: information concerning the conventional cell, with angles in degrees
    """
    convParamsAndInfo = getConventional(S)
    alphaInRadian = convParamsAndInfo["alpha"]
    betaInRadian = convParamsAndInfo["beta"]
    gammaInRadian = convParamsAndInfo["gamma"]

    alphaInDegree = math.degrees(alphaInRadian)
    betaInDegree = math.degrees(betaInRadian)
    gammaInDegree = math.degrees(gammaInRadian)

    convParamsAndInfo["alpha"] = alphaInDegree
    convParamsAndInfo["beta"] = betaInDegree
    convParamsAndInfo["gamma"] = gammaInDegree

    return convParamsAndInfo


def prim2convWithS(a, b, c, alpha, beta, gamma):
    """

    :param a: length
    :param b: length
    :param c: length
    :param alpha: angle between b and c, in degrees
    :param beta: angle between c and a, in degrees
    :param gamma: angle between a and b, in degrees
    :return: information about the conventional cell, with angles in degrees
    """
    # input parameters (angle in degrees) to matrix S
    S = inputParametersInDegree2S(a, b, c, alpha, beta, gamma)

    # find Buerger cell
    SofB = AlgorithmB(S)
    # [[A,B,C],[_,_,_]]=SofB
    # print("Lengths after algorithm B: "+str(np.sqrt(A))+", "+str(np.sqrt(B))+", "+str(np.sqrt(C)))

    # find all Buerger cells
    normalizedBgChars = allNormalizedBuergerCharacteristics(SofB)

    # find Niggli cell
    SofN = Buerger2Niggli(normalizedBgChars)
    checkAnglesOfS(SofN)

    checkNiggliCell(SofN)

    # Niggli cell to conventional cell (angle in degrees)
    convParamsAndInfo = getConventionalInDegree(SofN)

    return convParamsAndInfo


def chooseLen5Matrix(SofVec, normalizedBgChars):
    """
    there are 5 Buerger cells
    :param SofVec: S matrix of aVec, bVec, cVec
    :param normalizedBgChars: all possible Buerger cells corresponding to basis aVec, bVec, cVec
    :return: category number, row number, transformation matrix
    """
    epsRel = 1e-8
    epsAbs = 1e-6
    category = 20
    [[A, _, _], [_, _, _]] = SofVec
    # to determine each matrix of normalizedBgChars is which row in  category 20

    # first partition the matrices in normalizedBgChars by the value of zeta
    positiveA = []
    negativeA = []
    notA = []
    for mat in normalizedBgChars:
        if np.isclose(mat[1, 2], A, rtol=epsRel, atol=epsAbs, equal_nan=False):
            positiveA.append(mat)
        elif np.isclose(mat[1, 2], -A, rtol=epsRel, atol=epsAbs, equal_nan=False):
            negativeA.append(mat)
        else:
            notA.append(mat)

    if not (len(positiveA) == 2 and len(negativeA) == 2 and len(notA) == 1):
        raise RuntimeError("Not correct input with 5 Buerger matrices.")
    # determine the row number of  matrices in positiveA in category 20, by using A-q/2>q/2
    S0, S1 = positiveA
    if S0[1, 1] > S1[1, 1]:
        row1Mat = S0
        row2Mat = S1
    else:
        row1Mat = S1
        row2Mat = S0
    # determine the row number of matrices in negativeA in category 20, by using A-q/2>q/2
    S2, S3 = negativeA
    if S2[1, 0] > S3[1, 0]:
        row5Mat = S2
        row4Mat = S3
    else:
        row5Mat = S3
        row4Mat = S2
    row3Mat = notA[0]

    eps = 1e-5
    if np.linalg.norm(SofVec - row1Mat, ord=2) <= eps:
        M = np.identity(3, dtype=np.float64)
        return category, 1, M
    if np.linalg.norm(SofVec - row2Mat, ord=2) <= eps:
        M = np.array([[1, 0, 0], [0, 1, 0], [0, 1, -1]], dtype=np.float64)
        return category, 2, M
    if np.linalg.norm(SofVec - row3Mat, ord=2) <= eps:
        M = np.array([[1, 0, 0], [0, -1, -1], [0, -1, 0]], dtype=np.float64)
        return category, 3, M
    if np.linalg.norm(SofVec - row4Mat, ord=2) <= eps:
        M = np.array([[1, 0, 0], [1, 1, 0], [1, 1, 1]])
        return category, 4, M
    if np.linalg.norm(SofVec - row5Mat, ord=2) <= eps:
        M = np.array([[1, 0, 0], [1, 1, 0], [0, 0, -1]])
        return category, 5, M

    raise RuntimeError("SofVec not in normalizedBgChars with length 5")


def chooseLen4Matrix(SofVec, normalizedBgChars):
    """
        there are 4 Buerger cells
        :param SofVec: S matrix of aVec, bVec, cVec
        :param normalizedBgChars: all possible Buerger cells corresponding to basis aVec, bVec, cVec
        :return: category number, row number, transformation matrix
        """
    epsRel = 1e-8
    epsAbs = 1e-6
    category = 23
    [[A, _, _], [_, _, _]] = SofVec

    # first partition the matrices in normalizedBgChars by the value of zeta
    positiveA = []
    negativeA = []
    for mat in normalizedBgChars:
        if np.isclose(mat[1, 2], A, rtol=epsRel, atol=epsAbs, equal_nan=False):
            positiveA.append(mat)
        elif np.isclose(mat[1, 2], -A, rtol=epsRel, atol=epsAbs, equal_nan=False):
            negativeA.append(mat)
    if not (len(positiveA) == 2 and len(negativeA) == 2):
        raise RuntimeError("Not correct input with 4 Buerger matrices.")
    #determine the row of matrices in category 23, by using A-q/2>q/2
    S0,S1=positiveA
    if S0[1,1]>S1[1,1]:
        row1Mat=S0
        row2Mat=S1
    else:
        row1Mat = S1
        row2Mat = S0
    S2,S3=negativeA
    if S2[1,1]>S3[1,1]:
        row3Mat=S2
        row4Mat=S3
    else:
        row3Mat = S3
        row4Mat = S2
    eps = 1e-5
    if np.linalg.norm(SofVec-row1Mat,ord=2)<=eps:
        M = np.identity(3, dtype=np.float64)
        return category, 1, M
    if np.linalg.norm(SofVec - row2Mat, ord=2) <= eps:
        M = np.array([[1,0,0],[0,1,0],[0,1,-1]], dtype=np.float64)
        return category, 2, M
    if np.linalg.norm(SofVec - row3Mat, ord=2) <= eps:
        M = np.array([[1,0,0],[1,1,0],[1,1,1]], dtype=np.float64)
        return category, 3, M
    if np.linalg.norm(SofVec - row4Mat, ord=2) <= eps:
        M = np.array([[1,0,0],[1,1,0],[0,0,-1]])
        return category, 4, M
    raise RuntimeError("SofVec not in normalizedBgChars with length 4")

def chooseLen3Category3Matrix(SofVec, normalizedBgChars):
    """

    :param SofVec: S matrix of aVec, bVec, cVec
    :param normalizedBgChars: all possible Buerger cells corresponding to basis aVec, bVec, cVec
    :return: row number, transformation matrix for category 3
    """
    epsRel = 1e-8
    epsAbs = 1e-6
    [[A, _, _], [_, _, _]] = SofVec
    # first partition the matrices in normalizedBgChars by the value of zeta
    positiveA = []
    negativeA = []
    for mat in normalizedBgChars:
        if np.isclose(mat[1, 2], A, rtol=epsRel, atol=epsAbs, equal_nan=False):
            positiveA.append(mat)
        elif np.isclose(mat[1, 2], -A, rtol=epsRel, atol=epsAbs, equal_nan=False):
            negativeA.append(mat)
    if not (len(positiveA) == 2 and len(negativeA) == 1):
        raise RuntimeError("Not correct input with 3 Buerger matrices, category 3.")
    #in positiveA, using q-p/2>p/2
    S0,S1=positiveA
    if S0[1,0]>S1[1,0]:
        row1Mat=S0
        row2Mat=S1
    else:
        row1Mat = S1
        row2Mat = S0
    row3Mat=negativeA[0]
    eps = 1e-5
    if np.linalg.norm(SofVec - row1Mat, ord=2) <= eps:
        M = np.identity(3, dtype=np.float64)
        return 1, M
    if np.linalg.norm(SofVec - row2Mat, ord=2) <= eps:
        M = np.array([[1,0,0],[1,-1,0],[0,0,1]],dtype=np.float64)
        return 2, M
    if np.linalg.norm(SofVec - row3Mat, ord=2) <= eps:
        M = np.array([[1,1,0],[1,0,0],[0,0,-1]],dtype=np.float64)
        return 3, M

    raise RuntimeError("SofVec not in normalizedBgChars for category 3")




def partitionLen2(normalizedBgChars):
    """

    :param normalizedBgChars: all possible Buerger cells with length 2
    :return: possible categories for 2 Buerger cells
    """
    [[A,B,C],[_,_,_]]=normalizedBgChars[0]
    epsRel = 1e-8
    epsAbs = 1e-6
    S0,S1=normalizedBgChars
    #case A=B
    if np.isclose(A,B, rtol=epsRel, atol=epsAbs, equal_nan=False):
        # 1
        xi0=S0[1,0]
        xi1=S1[1,0]

        if (xi0>0 and (xi1<0 or np.isclose(xi1,0, rtol=epsRel, atol=epsAbs, equal_nan=False)))\
            or (xi1>0 and (xi0<0 or np.isclose(xi0,0, rtol=epsRel, atol=epsAbs, equal_nan=False))):
            return [1,2,4,6,7,8]
        #2
        eta0=S0[1,1]
        eta1=S1[1,1]
        if eta0>0 and eta1>0:
            return [5]
    #case A<B=C
    if A<B and np.isclose(B,C,rtol=epsRel, atol=epsAbs, equal_nan=False):
        # 1
        xi0 = S0[1, 0]
        xi1 = S1[1, 0]
        if (xi0 > 0 and (xi1 < 0 or np.isclose(xi1, 0, rtol=epsRel, atol=epsAbs, equal_nan=False))) \
                or (xi1 > 0 and (xi0 < 0 or np.isclose(xi0, 0, rtol=epsRel, atol=epsAbs, equal_nan=False))):
            return [1,10,12,14,15]
        # 2
        if xi0>0 and xi1>0:
            return [9]
        #3
        if xi0<0 and xi1<0:
            return [16,17]
    #case A<B<C
    if A<B and B<C:
        # 1
        xi0 = S0[1, 0]
        xi1 = S1[1, 0]
        if (xi0 > 0 and (xi1 < 0 or np.isclose(xi1, 0, rtol=epsRel, atol=epsAbs, equal_nan=False))) \
                or (xi1 > 0 and (xi0 < 0 or np.isclose(xi0, 0, rtol=epsRel, atol=epsAbs, equal_nan=False))):
            return [4, 10, 24, 25, 26, 27, 28]
        #2
        if xi0>0 and xi1>0:
            return [5,9,21]
        #
        return [22]

    raise RuntimeError("Wrong normalized matrices for length 2.")

def chooseLen2Partition1Matrix(SofVec,normalizedBgChars,ptn):
    """
    len(normalizedBgChars)==2
    :param SofVec
    :param normalizedBgChars: all 2 possible Buerger cells
    :param ptn: partition1 list=[1,2,4,6,7,8]
    :return: category, row, M
    """
    epsRel = 1e-8
    epsAbs = 1e-6
    if ptn!=[1,2,4,6,7,8]:
        raise ValueError("Wrong partition type for [1,2,4,6,7,8].")
    #category 1, 4
    S0,S1=normalizedBgChars
    [[A,B,C],[_,_,_]]=S0
    if np.isclose(S0[1,0],0,rtol=epsRel, atol=epsAbs, equal_nan=False)\
            or np.isclose(S1[1,0],0,rtol=epsRel, atol=epsAbs, equal_nan=False):
        #TODO: for case q=A, category 1 and 4 is not distinguishable!
        if np.isclose(B,C,rtol=epsRel, atol=epsAbs, equal_nan=False):
            if np.isclose(SofVec[1,0],0,rtol=epsRel, atol=epsAbs, equal_nan=False):
                M=np.array([[1,0,0],[1,1,0],[0,0,-1]])
                row=2
                return 1, row, M
            if SofVec[1,0]>0:
                M=np.identity(3,dtype=np.float64)
                row=1
                return 1,row,M

        if np.abs(S0[1,2])>np.abs(S1[1,1]):
            if np.isclose(SofVec[1,0],0,rtol=epsRel, atol=epsAbs, equal_nan=False):
                M=np.array([[1,0,0],[1,1,0],[0,0,-1]])
                row=2
                return 1, row, M
            if SofVec[1,0]>0:
                M=np.identity(3,dtype=np.float64)
                row=1
                return 1,row,M
        if np.abs(S1[1,1])>np.abs(S0[1,2]):
            if np.isclose(SofVec[1, 0], 0, rtol=epsRel, atol=epsAbs, equal_nan=False):
                M=np.array([[1,0,0],[0,-1,0],[1,0,1]],dtype=np.float64)
                row=2
                return 4,row,M
            if SofVec[1,0]>0:
                M=np.identity(3,dtype=np.float64)
                row=1
                return 4, row,M


        raise RuntimeError("Uncomparable case for category 1 and 4")
    #category 2
    if not np.isclose(np.abs(S0[1,1]),np.abs(S1[1,1]),rtol=epsRel, atol=epsAbs, equal_nan=False):
        if SofVec[1,0]>0:
            M=np.identity(3,dtype=np.float64)
            row=1
            return 2,row,M
        if SofVec[1,0]<0:
            M=np.array([[1,1,0],[0,1,0],[0,0,-1]],dtype=np.float64)
            row=2
            return 2, row, M

    #category 6,7,8
    #category 6
    minAbs,maxAbs=sorted(np.abs([S0[1,0],S1[1,0]]))
    if np.isclose(maxAbs/minAbs,2,rtol=epsRel, atol=epsAbs, equal_nan=False):
        category=6
        if SofVec[1,0]>0:
            row=1
            M = np.identity(3, dtype=np.float64)
            return category,row,M
        elif SofVec[1,0]<0:
            row=2
            M=np.array([[1,0,0],[0,-1,0],[1,0,1]],dtype=np.float64)
            return category, row, M
        else:
            raise RuntimeError("Results not found for category 6.")
    elif maxAbs/minAbs<2:
        category =7
        if SofVec[1,0]>0:
            row=1
            M = np.identity(3, dtype=np.float64)
            return category, row, M
        elif SofVec[1,0]<0:
            row=2
            M=np.array([[1,0,0],[0,-1,0],[1,0,1]],dtype=np.float64)
            return category, row, M
        else:
            raise RuntimeError("Results not found for category 7.")
    elif maxAbs/minAbs>2:
        category = 8
        if SofVec[1,0]>0:
            row=1
            M = np.identity(3, dtype=np.float64)
            return category, row, M
        elif SofVec[1, 0] < 0:
            row=2
            M=np.array([[1,0,0],[0,-1,0],[1,0,1]],dtype=np.float64)
            return category, row, M
        else:
            raise RuntimeError("Results not found for category 8.")

    raise ValueError("Buerger matrices for partition [1,2,4,6,7,8] is wrong.")

def chooseLen2Partition2Matrix(SofVec,normalizedBgChars,ptn):
    """
    len(normalizedBgChars)==2
    :param SofVec:
    :param normalizedBgChars:   all 2 possible Buerger cells
    :param ptn: partition2 list=[5]
    :return: category, row, M
    """
    epsRel = 1e-8
    epsAbs = 1e-6
    if ptn != [5]:
        raise ValueError("Wrong partition type for [5].")
    category=5
    SMin,SMax=sorted(normalizedBgChars,key=lambda mat: mat[1,0])
    if np.isclose(SofVec[1,0],SMin[1,0],rtol=epsRel, atol=epsAbs, equal_nan=False):
        row=2
        M=np.array([[1,0,0],[0,1,0],[1,0,-1]],dtype=np.float64)
        return category,row,M
    elif np.isclose(SofVec[1,0],SMax[1,0],rtol=epsRel, atol=epsAbs, equal_nan=False):
        row=1
        M = np.identity(3, dtype=np.float64)
        return category, row, M
    else:
        raise RuntimeError("Category 5 results not found.")


def chooseLen2Partition3Matrix(SofVec,normalizedBgChars,ptn):
    """
    len(normalizedBgChars)==2
    :param SofVec:
    :param normalizedBgChars: all 2 possible Buerger cells
    :param ptn: partition3 list=[1,10,12,14,15]
    :return: category, row, M
    """
    epsRel = 1e-8
    epsAbs = 1e-6
    if ptn != [1,10,12,14,15]:
        raise ValueError("Wrong partition type for [1,10,12,14,15].")
    S0, S1 = normalizedBgChars
    #category 15
    if not np.isclose(np.abs(S0[1,2]),np.abs(S1[1,2]),rtol=epsRel, atol=epsAbs, equal_nan=False):
        category=15
        if SofVec[1,0]>0:
            row=1
            M= np.identity(3, dtype=np.float64)
            return category,row,M
        elif SofVec[1,0]<0:
            row=2
            M=np.array([[1,0,0],[0,-1,-1],[0,0,-1]],dtype=np.float64)
            return category, row, M
        else:
            raise RuntimeError("Results not found for category 15.")
    # category 1
    if np.isclose(S0[1,0],0,rtol=epsRel, atol=epsAbs, equal_nan=False)\
        or np.isclose(S1[1,0],0,rtol=epsRel, atol=epsAbs, equal_nan=False):
        category = 1
        if SofVec[1,2]>0:
            row=1
            M = np.identity(3, dtype=np.float64)
            return category, row, M
        elif SofVec[1,2]<0:
            row=2
            M=np.array([[1,0,0],[1,1,0],[0,0,-1]],dtype=np.float64)
            return category, row, M
        else:
            raise RuntimeError("Results not found for category 1 when A<B=C.")
    # category 14
    if np.isclose(S0[1,1],0,rtol=epsRel, atol=epsAbs, equal_nan=False)\
        or np.isclose(S1[1,1],0,rtol=epsRel, atol=epsAbs, equal_nan=False):
        category = 14
        if SofVec[1,0]>0:
            row=1
            M = np.identity(3, dtype=np.float64)
            return category, row, M
        elif SofVec[1,0]<0:
            row=2
            M=np.array([[1,0,0],[0,-1,0],[0,-1,-1]],dtype=np.float64)
            return category, row, M
        else:
            raise RuntimeError("Results not found for category 14 when A<B=C.")
    #category 10,12, they have common transformation matrices, no further distinction is needed
    category=12
    if SofVec[1,2]>0:
        row=1
        M = np.identity(3, dtype=np.float64)
        return category, row, M
    elif SofVec[1,2]<0:
        row=2
        M=np.array([[1,0,0],[1,1,0],[0,0,-1]],dtype=np.float64)
        return category, row, M

    raise RuntimeError("Results not found for partition [1,10,12,14,15].")



def chooseLen2Partition4Matrix(SofVec,normalizedBgChars,ptn):
    """
    len(normalizedBgChars)==2
    :param SofVec:
    :param normalizedBgChars: all 2 possible Buerger cells
    :param ptn: partition4 list=[9]
    :return: category, row, M
    """
    epsRel = 1e-8
    epsAbs = 1e-6
    category=9
    if ptn != [9]:
        raise ValueError("Wrong partition type for [9].")
    SMin,SMax=sorted(normalizedBgChars,key=lambda mat:mat[1,0])
    if np.isclose(SofVec[1,0],SMin[1,0],rtol=epsRel, atol=epsAbs, equal_nan=False):
        row=2
        M=np.array([[1,0,0],[1,-1,0],[0,0,1]],dtype=np.float64)
        return category, row, M
    if np.isclose(SofVec[1,0],SMax[1,0],rtol=epsRel, atol=epsAbs, equal_nan=False):
        row=1
        M = np.identity(3, dtype=np.float64)
        return category, row, M
    raise RuntimeError("Results not found for partition [9].")


def chooseLen2Partition5Matrix(SofVec,normalizedBgChars,ptn):
    """
    len(normalizedBgChars)==2
    :param SofVec:
    :param normalizedBgChars: all 2 possible Buerger cells
    :param ptn: partition5 list=[16,17]
    :return: category, row, M
    """
    epsRel = 1e-8
    epsAbs = 1e-6
    if ptn != [16,17]:
        raise ValueError("Wrong partition type for [16,17].")
    #sort by the value of eta
    SEtaMin,SEtaMax=sorted(normalizedBgChars,key=lambda mat:mat[1,1])
    if np.isclose(SEtaMin[1,2],SEtaMax[1,2],rtol=epsRel, atol=epsAbs, equal_nan=False):
        category=16
        if np.isclose(SofVec[1,1],SEtaMin[1,1],rtol=epsRel, atol=epsAbs, equal_nan=False):
            row=1
            M = np.identity(3, dtype=np.float64)
            return category, row, M
        elif np.isclose(SofVec[1,1],SEtaMax[1,1],rtol=epsRel, atol=epsAbs, equal_nan=False):
            row=2
            M=np.array([[1,0,0],[0,1,0],[-1,-1,-1]],dtype=np.float64)
            return category, row, M
        else:
            raise RuntimeError("Results not found for category 16 when A<B=C.")
    #sort by the value of zeta
    SZetaMin,SZetaMax=sorted(normalizedBgChars,key=lambda mat: mat[1,2])
    if np.isclose(SZetaMin[1,1],SZetaMax[1,1],rtol=epsRel, atol=epsAbs, equal_nan=False):
        category=17
        if np.isclose(SofVec[1,2],SZetaMin[1,2],rtol=epsRel, atol=epsAbs, equal_nan=False):
            row=1
            M = np.identity(3, dtype=np.float64)
            return category, row, M
        elif np.isclose(SofVec[1,2],SZetaMax[1,2],rtol=epsRel, atol=epsAbs, equal_nan=False):
            row=2
            M=np.array([[1,0,0],[-1,-1,-1],[0,0,1]],dtype=np.float64)
            return category, row, M
        else:
            raise RuntimeError("Results not found for category 17 when A<B=C.")

    raise RuntimeError("Results not found for partition [16,17].")





def chooseLen2Partition6Matrix(SofVec,normalizedBgChars,ptn):
    """
    len(normalizedBgChars)==2
    :param SofVec:
    :param normalizedBgChars: all 2 possible Buerger cells
    :param ptn: partition6 list=[4, 10, 24, 25, 26, 27, 28]
    :return: category, row, M
    """
    epsRel = 1e-8
    epsAbs = 1e-6
    if ptn !=[4, 10, 24, 25, 26, 27, 28]:
        raise ValueError("Wrong partition type for [4, 10, 24, 25, 26, 27, 28].")
    S0,S1=normalizedBgChars
    #category 4,28, when one xi=0

    if np.isclose(S0[1,0],0,rtol=epsRel, atol=epsAbs, equal_nan=False)\
        or np.isclose(S1[1,0],0,rtol=epsRel, atol=epsAbs, equal_nan=False):

        SXiMin,SXiMax=sorted(normalizedBgChars,key=lambda mat:mat[1,0])
        # when q=A
        if np.isclose(SXiMax[1,0],SXiMax[1,1],rtol=epsRel, atol=epsAbs, equal_nan=False):
            if np.isclose(SXiMax[1,1],SXiMax[1,2],rtol=epsRel, atol=epsAbs, equal_nan=False):
                category=4
                if np.isclose(SofVec[1,0],SXiMax[1,0],rtol=epsRel, atol=epsAbs, equal_nan=False):
                    row=1
                    M = np.identity(3, dtype=np.float64)
                    return category, row, M
                elif np.isclose(SofVec[1,0],SXiMin[1,0],rtol=epsRel, atol=epsAbs, equal_nan=False):
                    row=2
                    M=np.array([[1,0,0],[0,-1,0],[1,0,1]],dtype=np.float64)
                    return category, row, M
                else:
                    raise RuntimeError("Results not found for category 4 when q=A<B<C.")
            if SXiMax[1,1]<SXiMax[1,2]:
                category=28
                if np.isclose(SofVec[1,0],SXiMax[1,0],rtol=epsRel, atol=epsAbs, equal_nan=False):
                    row=1
                    M = np.identity(3, dtype=np.float64)
                    return category, row, M
                elif np.isclose(SofVec[1,0],SXiMin[1,0],rtol=epsRel, atol=epsAbs, equal_nan=False):
                    row=2
                    M=np.array([[1,0,0],[1,1,0],[0,0,-1]],dtype=np.float64)
                    return category, row, M
                else:
                    raise RuntimeError("Results not found for category 28 when q=A<B<C.")
        # when q<A
        elif SXiMax[1,0]<SXiMax[1,1]:
            category=4
            if np.isclose(SofVec[1,0],SXiMax[1,0],rtol=epsRel, atol=epsAbs, equal_nan=False):
                row=1
                M = np.identity(3, dtype=np.float64)
                return category, row, M
            elif np.isclose(SofVec[1,0],SXiMin[1,0],rtol=epsRel, atol=epsAbs, equal_nan=False):
                row=2
                M=np.array([[1,0,0],[0,-1,0],[1,0,1]],dtype=np.float64)
                return category, row, M
            else:
                raise RuntimeError("Results not found for category 4 when q<A<B<C.")
        else:
            raise RuntimeError("Results not found for category 4,28, one  xi =0.")

    #category 10,24,25,26,27,28, no xi=0

    SZetaMin,SZetaMax=sorted(normalizedBgChars,key=lambda mat:mat[1,2])
    # category 10,24,28
    if SZetaMax[1,1]<SZetaMax[1,2]:
        #transformation matrices are the same for category 10,24,28
        category=10
        if np.isclose(SofVec[1,2],SZetaMax[1,2],rtol=epsRel, atol=epsAbs, equal_nan=False):
            row=1
            M = np.identity(3, dtype=np.float64)
            return category, row, M
        elif np.isclose(SofVec[1,2],SZetaMin[1,2],rtol=epsRel, atol=epsAbs, equal_nan=False):
            row=2
            M=np.array([[1,0,0],[1,1,0],[0,0,-1]],dtype=np.float64)
            return category, row, M
        else:
            raise RuntimeError("Results not found for category 10,24,28.")
    #category 25,26,27
    else:
        sumOfXi=SZetaMax[1,0]+SZetaMin[1,0]
        #category 25
        if np.isclose(sumOfXi,0,rtol=epsRel, atol=epsAbs, equal_nan=False):
            category=25
            if np.isclose(SofVec[1,0],SZetaMax[1,0],rtol=epsRel, atol=epsAbs, equal_nan=False):
                row=1
                M = np.identity(3, dtype=np.float64)
                return category, row, M
            elif np.isclose(SofVec[1,0],SZetaMin[1,0],rtol=epsRel, atol=epsAbs, equal_nan=False):
                row=2
                M=np.array([[1,0,0],[0,-1,0],[0,-1,-1]],dtype=np.float64)
                return category, row, M
            else:
                raise RuntimeError("Results not found for category 25.")
        #category 26,27, no further distinction is needed
        elif sumOfXi>0:
            category=26
            if np.isclose(SofVec[1,1],SZetaMax[1,1],rtol=epsRel, atol=epsAbs, equal_nan=False):
                row = 1
                M = np.identity(3, dtype=np.float64)
                return category, row, M
            elif np.isclose(SofVec[1,1],SZetaMin[1,1],rtol=epsRel, atol=epsAbs, equal_nan=False):
                row=2
                M=np.array([[1,0,0],[0,-1,0],[1,0,1]],dtype=np.float64)
                return category, row, M
        else:
            raise RuntimeError("Results not found for category 25,26,27.")




def chooseLen2Partition7Matrix(SofVec,normalizedBgChars,ptn):
    """
    len(normalizedBgChars)==2
    :param SofVec:
    :param normalizedBgChars: all 2 possible Buerger cells
    :param ptn: partition7 list=[5,9,21]
    :return: category, row, M
    """
    epsRel = 1e-8
    epsAbs = 1e-6
    if ptn !=[5,9,21]:
        raise ValueError("Wrong partition type for [5,9,21].")
    SXiMin,SXiMax=sorted(normalizedBgChars,key=lambda mat:mat[1,0])
    #category 21
    if np.isclose(SXiMin[1,0],SXiMax[1,0],rtol=epsRel, atol=epsAbs, equal_nan=False):
        category=21
        SEtaMin,SEtaMax=sorted(normalizedBgChars,key=lambda mat:mat[1,1])
        if np.isclose(SofVec[1,1],SEtaMax[1,1],rtol=epsRel, atol=epsAbs, equal_nan=False):
            row = 1
            M = np.identity(3, dtype=np.float64)
            return category, row, M
        elif np.isclose(SofVec[1,1],SEtaMin[1,1],rtol=epsRel, atol=epsAbs, equal_nan=False):
            row=2
            M=np.array([[1,0,0],[0,1,0],[0,1,-1]],dtype=np.float64)
            return category, row, M
        else:
            raise RuntimeError("Results not found for category 21.")

    #category 9
    if SXiMax[1,1]<SXiMax[1,2]:
        category=9
        if np.isclose(SofVec[1,0],SXiMax[1,0],rtol=epsRel, atol=epsAbs, equal_nan=False):
            row = 1
            M = np.identity(3, dtype=np.float64)
            return category, row, M
        elif np.isclose(SofVec[1,0],SXiMin[1,0],rtol=epsRel, atol=epsAbs, equal_nan=False):
            row=2
            M=np.array([[1,0,0],[1,-1,0],[0,0,1]],dtype=np.float64)
            return category, row, M
        else:
            raise RuntimeError("Results not found for category 9.")

    #category 5
    category=5
    if np.isclose(SofVec[1,0],SXiMax[1,0],rtol=epsRel, atol=epsAbs, equal_nan=False):
        row = 1
        M = np.identity(3, dtype=np.float64)
        return category, row, M
    elif np.isclose(SofVec[1,0],SXiMin[1,0],rtol=epsRel, atol=epsAbs, equal_nan=False):
        row=2
        M=np.array([[1,0,0],[0,1,0],[1,0,-1]],dtype=np.float64)
        return category, row, M
    else:
        raise RuntimeError("Results not found for category 5.")



def chooseLen2Partition8Matrix(SofVec,normalizedBgChars,ptn):
    """
    len(normalizedBgChars)==2
    :param SofVec:
    :param normalizedBgChars: all 2 possible Buerger cells
    :param ptn: partition8 list=[22]
    :return: category, row, M
    """
    epsRel = 1e-8
    epsAbs = 1e-6
    if ptn !=[22]:
        raise ValueError("Wrong partition type for [22].")
    category=22
    SXiMin,SXiMax=sorted(normalizedBgChars,key=lambda mat:mat[1,0])

    if np.isclose(SofVec[1,0],SXiMax[1,0],rtol=epsRel, atol=epsAbs, equal_nan=False):
        row = 1
        M = np.identity(3, dtype=np.float64)
        return category, row, M
    elif np.isclose(SofVec[1,0],SXiMin[1,0],rtol=epsRel, atol=epsAbs, equal_nan=False):
        row=2
        M=np.array([[1,0,0],[0,1,0],[-1,-1,-1]],dtype=np.float64)
        return category, row, M
    else:
        raise RuntimeError("Results not found for category 22.")















def Buerger2NiggliBasis(SofVec, normalizedBgChars):
    """

    :param SofVec: S matrix of aVec, bVec, cVec
    :param normalizedBgChars: all possible Buerger cells corresponding to basis aVec, bVec, cVec
    :return: basis of Niggli cell
    """
    # if there is only 1 Buerger cell, then this is the Niggli cell
    if len(normalizedBgChars) == 1:
        return [np.identity(3, dtype=np.float64), []]
    # if there are 5 Buerger cells, they belong to category 20
    if len(normalizedBgChars) == 5:
        category, row, M = chooseLen5Matrix(SofVec, normalizedBgChars)
        return [M,[category,row]]

    # if there are 4 Buerger cells, they belong to category 23
    if len(normalizedBgChars)==4:
        category, row, M = chooseLen4Matrix(SofVec, normalizedBgChars)
        return [M, [category, row]]
    #if there are 3 Buerger cells
    if len(normalizedBgChars)==3:
        category=selectLen3(normalizedBgChars)
        row,M=chooseLen3Category3Matrix(SofVec, normalizedBgChars)#TODO: not finished
        return [M,[category,row]]

    #if there are 2 Buerger cells
    if len(normalizedBgChars)==2:
        partitionList=partitionLen2(normalizedBgChars)
        if partitionList==[1,2,4,6,7,8]:
            category,row,M=chooseLen2Partition1Matrix(SofVec,normalizedBgChars,partitionList)
            return [M,[category,row]]
        if partitionList==[5]:
            category,row,M=chooseLen2Partition2Matrix(SofVec,normalizedBgChars,partitionList)
            return [M, [category, row]]
        if partitionList==[1,10,12,14,15]:
            category, row, M=chooseLen2Partition3Matrix(SofVec,normalizedBgChars,partitionList)
            return [M, [category, row]]
        if partitionList==[9]:
            category, row, M = chooseLen2Partition4Matrix(SofVec,normalizedBgChars,partitionList)
            return [M, [category, row]]
        if partitionList==[16,17]:
            category, row, M = chooseLen2Partition5Matrix(SofVec,normalizedBgChars,partitionList)
            return [M, [category, row]]
        if partitionList==[4, 10, 24, 25, 26, 27, 28]:
            category, row, M = chooseLen2Partition6Matrix(SofVec,normalizedBgChars,partitionList)
            return [M, [category, row]]
        if partitionList==[5,9,21]:
            category, row, M = chooseLen2Partition7Matrix(SofVec,normalizedBgChars,partitionList)
            return [M, [category, row]]
        if partitionList==[22]:
            category, row, M = chooseLen2Partition8Matrix(SofVec,normalizedBgChars,partitionList)
            return [M, [category, row]]










def prim2convWithBasis(aVec, bVec, cVec):
    """

    :param aVec: basis vector of primitive cell
    :param bVec: basis vector of primitive cell
    :param cVec: basis vector of primitive cell
    :return: information about the conventional cell, with angles in degrees
    """

    SIn = vec2S(aVec, bVec, cVec)  # for sanity check

    aVecOfB, bVecOfB, cVecOfB = AlgorithmBVectorTransform(aVec, bVec, cVec)

    SofVec = vec2S(aVecOfB, bVecOfB, cVecOfB)  # basis for one Buerger cell
    # find all Buerger cells
    normalizedBgChars = allNormalizedBuergerCharacteristics(SofVec)

# test data for aP #passed
# a=np.array([ 0.23901248, -0.24064142,  0.64521519])
# b=np.array( [ 0.88164565, -0.46578039,  0.22768941])
# c=np.array([ 1.0949293,  -0.39030657 , 0.76405586])
# test data for mP #passed
# a=np.array([1.80000000e+00,0.00000000e+00,0.00000000e+00])
# b=np.array([0.00000000e+00,1.00000000e+00,0.00000000e+00])
# b=a+b
# c=np.array([2.26182247e-17,3.69383642e-01,3.68151541e+00])

# #test data for mS #passed
# a=np.array([-0.0562122,  -0.13247576,  0.64146155])
# b=np.array([-0.47766339 , 0.43155659 , 0.13333848])
# c=np.array( [-0.01755021,  0.61170953,  0.69357314])

# test data for oP
# direct primitive #passed
# a=np.array([ 0.24182656,0.07482303,0.29679825])
# b=np.array([ 0.3111148,-0.51866721,-0.12273512])
# c=np.array( [ 0.24497009,0.20649208,-0.25165457])

# #diagonal primitive#passed
# a=np.array([ 0.24182656,0.07482303,0.29679825])
# b=np.array([ 0.3111148,-0.51866721,-0.12273512])
#
# c=np.array( [ 0.24497009,0.20649208,-0.25165457])
# c=b+c


# test data for oI
# direct primitive #passed
# a=np.array([0.11821031,-0.36746163,0.35876362])
# b=np.array([0.28558464,0.30117927,-0.32472376])
# c=np.array([-0.15474425,0.44772062,0.2308972])
# diagonal primitive#passed
# a=np.array([0.11821031,-0.36746163,0.35876362])
# b=np.array([0.28558464,0.30117927,-0.32472376])
# b=a+b
# c=np.array([-0.15474425,0.44772062,0.2308972])

# test data for oS
# # direct primitive # passed
# test data 1
# a=np.array([-0.03336087,0.32522261,  0.28289925])
# b=np.array([ 0.20039596, -0.23046743, -0.30600801])
# c=np.array([ 0.36131684, -0.48934562,  0.60516279])
# test data 2
# vec1=np.array([1,0,0])
# vec2=np.array([0,1,0])*1.2
# vec3=np.array([0,0,1])*1.8
# a=1/2*vec1+1/2*vec2
# b=1/2*vec1-1/2*vec2
# c=vec3


# test data for oF
# direct primitive #passed
# a=np.array([ 0.00074449 , 0.34546068 , 0.16319794])
# b=np.array([ 0.16904866, -0.00225161,  0.1483042 ])
# c=np.array( [ 0.1603072 ,  0.35175292 , 0.00484137])

# test data for tP
# direct primitive #passed
# a=np.array([-0.01177838,  0.02464893,  0.03547314])
# b=np.array([ 0.04023329, -0.00712276,  0.01830825])
# c=np.array([ 0.22313218,  0.52073779, -0.28775274])
# diagonal primitive#passed
# a=np.array([-0.01177838,  0.02464893,  0.03547314])
# b=np.array([ 0.04023329, -0.00712276,  0.01830825])
# b=a+b
# c=np.array([ 0.22313218,  0.52073779, -0.28775274])


# test data for tI
# direct primitive #passed
# a=np.array([ 0.25168644,  0.13322825, -0.0714029 ])
# b=np.array([-0.28808599, -0.05343177, -0.01859346])
# c=np.array([-0.03187261,  0.16103162,  0.24340739])
# diagonal primitive##passed
# a=np.array([ 0.25168644,  0.13322825, -0.0714029 ])
# b=np.array([-0.28808599, -0.05343177, -0.01859346])
# b=a+b
# c=np.array([-0.03187261,  0.16103162,  0.24340739])

# test data for hR
# direct primitive#passed
# a=np.array([-0.98242063,  0.09226284,  0.62566037])
# b=np.array([ 0.3316406,   1.09929108,  0.21607231])
# c=np.array([-0.52470639,0.27509317,-1.00703553])
# diagonal primitive#passed
# a=np.array([-0.98242063,  0.09226284,  0.62566037])
# b=np.array([ 0.3316406,   1.09929108,  0.21607231])
# b=a+b
# c=np.array([-0.52470639,0.27509317,-1.00703553])


# test data for hP
# direct primitive# passed
# a=np.array([ 0.30446612,  0.56274179,  0.10568754])
# b=np.array([ 0.30516758, -0.56044481,  0.11542165])
# c=np.array( [ 0.087114 ,  -0.002027 ,  -0.24016629])

# test data for cP
# direct primitive#passed
# a=np.array([ 0.35027557, -0.08592507,  0.02108803])
# b=np.array([ 0.08777518,  0.32668059, -0.12687042])
# c=np.array([ 0.01110589,  0.12813067,  0.33760922])
# diagonal primitive#passed
# a=np.array([ 0.35027557, -0.08592507,  0.02108803])
# b=np.array([ 0.08777518,  0.32668059, -0.12687042])
# b=a+b
# c=np.array([ 0.01110589,  0.12813067,  0.33760922])

# test data for cI
# direct primitive#passed
# a=np.array([-0.15960944,  0.17685658,  0.13952733])
# b=np.array([ 0.14663659, -0.15176054,  0.17801053])
# c=np.array( [ 0.17166328,  0.14075557, -0.16413678])
# diagonal primitive#passed
# a=np.array([-0.15960944,  0.17685658,  0.13952733])
# b=np.array([ 0.14663659, -0.15176054,  0.17801053])
# b=a+b
# c=np.array( [ 0.17166328,  0.14075557, -0.16413678])

# test data for cF
# direct primitive #passed
# a=np.array([-0.02491603,  0.55194334, -0.10047505])
# b=np.array([ 0.13000506,  0.36546768,  0.40606594])
# c=np.array( [ 0.47296096,  0.30085982, -0.0338938 ])
# diagonal primitive #passed

# a=np.array([-0.02491603,  0.55194334, -0.10047505])
# b=np.array([ 0.13000506,  0.36546768,  0.40606594])
# b=a+b
# c=np.array( [ 0.47296096,  0.30085982, -0.0338938 ])

### use vector transformation part
# prim2convWithBasis(a,b,c)

####to input S matrix
# aLen=np.linalg.norm(a,ord=2)
# bLen=np.linalg.norm(b,ord=2)
# cLen=np.linalg.norm(c,ord=2)
# # print("original lenths are: "+str(aLen)+", "+str(bLen)+", "+str(cLen))
# alpha=np.arccos(c.dot(b)/(cLen*bLen))*180/np.pi
# beta=np.arccos(a.dot(c)/(aLen*cLen))*180/np.pi
# gamma=np.arccos(a.dot(b)/(aLen*bLen))*180/np.pi
#
# rst=prim2conv(aLen,bLen,cLen,alpha,beta,gamma)
#
# print(rst)
