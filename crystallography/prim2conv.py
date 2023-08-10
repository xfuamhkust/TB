import numpy as np
import math


# this script transforms a primitive cell to a conventional cell
#TODO: values of errors (epsilon) to be rectified

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

    #sanity check
    minAlpha=np.arccos(np.cos(betaInRadians-gammaInRadians))
    maxAlpha=np.arccos(np.cos(betaInRadians+gammaInRadians))
    if alphaInRadian<minAlpha or alphaInRadian>maxAlpha:
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
    if alpha< minAlpha or alpha > maxAlpha:
        raise ValueError("Invalid combination of angles.")

    S = np.array([[a ** 2, b ** 2, c ** 2],
                  [2 * b * c * np.cos(alpha), 2 * a * c * np.cos(beta),
                   2 * a * b * np.cos(gamma)]], dtype=np.float64)

    return S


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
        etaZetaInd = np.argsort(np.abs(S[1, [1, 2]]))+1#+1 because the actual indices are 1 and 2, but argsort returns 0 and 1
        S[1, [1, 2]] = S[1, etaZetaInd]

    # step 4, giving signs to xi, eta, zeta
    if S[1, 0] * S[1, 1] * S[1, 2] > 0:
        S[1, :] = np.abs(S[1, :])
    else:
        S[1, :] = -np.abs(S[1, :])
    # [[A, B, C], [xi, eta, zeta]] = S
    # a = np.sqrt(A)
    # b = np.sqrt(B)
    # c = np.sqrt(C)
    # cosAlpha = xi / (2 * b * c)
    # cosBeta = eta / (2 * a * c)
    # cosGamma = zeta / (2 * a * b)
    # print("After N cos(alpha)="+str(cosAlpha)+", cos(beta)="+str(cosBeta)+", cos(gamma)="+str(cosGamma))
    # print("N: "+str(S[0,0])+", "+str(S[0,1])+", "+str(S[0,2])+", "+str(S[1,0])+", "+str(S[1,1])+", "+str(S[1,2]))
    return S


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
    [[A,B,C],[xi,eta,zeta]]=S
    sumVal=A+B+C
    sigma=xi+eta+zeta

    epsRel = 1e-8
    epsAbs = 1e-6
    potentialBgChars=[]
    #test condition 1
    A1=A
    B1=B
    C1=C
    if np.isclose(A1+B1+C1, sumVal, rtol=epsRel, atol=epsAbs, equal_nan=False):
        xi1=xi
        eta1=eta
        zeta1=zeta
        potentialBgChars.append(np.array([[A1,B1,C1],[xi1,eta1,zeta1]], dtype=np.float64))

    #test condition 2
    A2=A
    B2=B
    C2=A+C+eta
    if np.isclose(A2+B2+C2,sumVal, rtol=epsRel, atol=epsAbs, equal_nan=False):
        xi2=xi+zeta
        eta2=2*A+eta
        zeta2=zeta
        potentialBgChars.append(np.array([[A2,B2,C2],[xi2,eta2,zeta2]], dtype=np.float64))

    #test condition 3
    A3=A
    B3=B
    C3=A+C-eta
    if np.isclose(A3+B3+C3,sumVal,rtol=epsRel, atol=epsAbs, equal_nan=False):
        xi3=-xi+zeta
        eta3=2*A-eta
        zeta3=zeta
        potentialBgChars.append(np.array([[A3,B3,C3],[xi3,eta3,zeta3]], dtype=np.float64))

    #test condition 4
    A4=A
    B4=B
    C4=B+C+xi
    if np.isclose(A4+B4+C4,sumVal,rtol=epsRel, atol=epsAbs, equal_nan=False):
        xi4=2*B+xi
        eta4=eta+zeta
        zeta4=zeta
        potentialBgChars.append(np.array([[A4,B4,C4],[xi4,eta4,zeta4]], dtype=np.float64))

    #test condition 5
    A5=A
    B5=B
    C5=B+C-xi
    if np.isclose(A5+B5+C5,sumVal,rtol=epsRel, atol=epsAbs, equal_nan=False):
        xi5=2*B-xi
        eta5=-eta+zeta
        zeta5=zeta
        potentialBgChars.append(np.array([[A5,B5,C5],[xi5,eta5,zeta5]], dtype=np.float64))

    #test condition 6
    A6=A
    B6=B
    C6=sumVal+sigma
    if np.isclose(A6+B6+C6,sumVal,rtol=epsRel, atol=epsAbs, equal_nan=False):
        xi6=2*B+xi+zeta
        eta6=2*A+eta+zeta
        zeta6=zeta
        potentialBgChars.append(np.array([[A6,B6,C6],[xi6,eta6,zeta6]], dtype=np.float64))

    #test condition 7
    A7=A
    B7=C
    C7=A+B+zeta
    if np.isclose(A7+B7+C7,sumVal,rtol=epsRel, atol=epsAbs, equal_nan=False):
        xi7=xi+eta
        eta7=2*A+zeta
        zeta7=eta
        potentialBgChars.append(np.array([[A7,B7,C7],[xi7,eta7,zeta7]], dtype=np.float64))

    #test condition 8
    A8=A
    B8=C
    C8=A+B-zeta
    if np.isclose(A8+B8+C8,sumVal,rtol=epsRel, atol=epsAbs, equal_nan=False):
        xi8=-xi+eta
        eta8=2*A-zeta
        zeta8=eta
        potentialBgChars.append(np.array([[A8,B8,C8],[xi8,eta8,zeta8]], dtype=np.float64))

    #test condition 9
    A9=A
    B9=C
    C9=B+C+xi
    if np.isclose(A9+B9+C9,sumVal,rtol=epsRel, atol=epsAbs, equal_nan=False):
        xi9=2*C+xi
        eta9=eta+zeta
        zeta9=eta
        potentialBgChars.append(np.array([[A9,B9,C9],[xi9,eta9,zeta9]], dtype=np.float64))

    #test condition 10
    A10=A
    B10=C
    C10=B+C-xi
    if np.isclose(A10+B10+C10,sumVal,rtol=epsRel, atol=epsAbs, equal_nan=False):
        xi10=-2*C+xi
        eta10=-eta+zeta
        zeta10=eta
        potentialBgChars.append(np.array([[A10,B10,C10],[xi10,eta10,zeta10]], dtype=np.float64))

    #test condition 11
    A11=A
    B11=C
    C11=sumVal+sigma
    if np.isclose(A11+B11+C11,sumVal,rtol=epsRel, atol=epsAbs, equal_nan=False):
        xi11=2*C+xi+eta
        eta11=2*A+eta+zeta
        zeta11=eta
        potentialBgChars.append(np.array([[A11,B11,C11],[xi11,eta11,zeta11]], dtype=np.float64))

    #test condition 12
    A12=A
    B12=A+B+zeta
    C12=B+C+xi
    if np.isclose(A12+B12+C12,sumVal,rtol=epsRel, atol=epsAbs, equal_nan=False):
        xi12=2*B+sigma
        eta12=eta+zeta
        zeta12=2*A+zeta
        potentialBgChars.append(np.array([[A12,B12,C12],[xi12,eta12,zeta12]], dtype=np.float64))

    #test condition 13
    A13=A
    B13=A+B+zeta
    C13=sumVal+sigma
    if np.isclose(A13+B13+C13,sumVal,rtol=epsRel, atol=epsAbs, equal_nan=False):
        xi13=2*A+2*B+sigma+zeta
        eta13=2*A+eta+zeta
        zeta13=2*A+zeta
        potentialBgChars.append(np.array([[A13,B13,C13],[xi13,eta13,zeta13]], dtype=np.float64))

    #test condition 14
    A14=A
    B14=A+B-zeta
    C14=B+C-xi
    if np.isclose(A14+B14+C14,sumVal,rtol=epsRel, atol=epsAbs, equal_nan=False):
        xi14=-2*B+xi-eta+zeta
        eta14=-eta+zeta
        zeta14=2*A-zeta
        potentialBgChars.append(np.array([[A14,B14,C14],[xi14,eta14,zeta14]], dtype=np.float64))

    #test condition 15
    A15=A
    B15=A+C+eta
    C15=sumVal+sigma
    if np.isclose(A15+B15+C15,sumVal,rtol=epsRel, atol=epsAbs, equal_nan=False):
        xi15=2*A+2*C+sigma+eta
        eta15=2*A+eta+zeta
        zeta15=2*A+eta
        potentialBgChars.append(np.array([[A15,B15,C15],[xi15,eta15,zeta15]], dtype=np.float64))

    #test condition 16
    A16=A
    B16=A+C-eta
    C16=B+C-xi
    if np.isclose(A16+B16+C16,sumVal,rtol=epsRel, atol=epsAbs, equal_nan=False):
        xi16=2*C-xi-eta+zeta
        eta16=-eta+zeta
        zeta16=2*A-eta
        potentialBgChars.append(np.array([[A16,B16,C16],[xi16,eta16,zeta16]], dtype=np.float64))

    #test condition 17
    A17=B
    B17=C
    C17=A+B+zeta
    if np.isclose(A17+B17+C17,sumVal,rtol=epsRel, atol=epsAbs, equal_nan=False):
        xi17=xi+eta
        eta17=2*B+zeta
        zeta17=xi
        potentialBgChars.append(np.array([[A17,B17,C17],[xi17,eta17,zeta17]], dtype=np.float64))

    #test condition 18
    A18=B
    B18=C
    C18=A+B-zeta
    if np.isclose(A18+B18+C18,sumVal,rtol=epsRel, atol=epsAbs, equal_nan=False):
        xi18=-xi+eta
        eta18=-2*B+zeta
        zeta18=xi
        potentialBgChars.append(np.array([[A18,B18,C18],[xi18,eta18,zeta18]], dtype=np.float64))

    #test condition 19
    A19=B
    B19=C
    C19=A+C-eta
    if np.isclose(A19+B19+C19,sumVal,rtol=epsRel, atol=epsAbs, equal_nan=False):
        xi19=-2*C+eta
        eta19=-xi+zeta
        zeta19=xi
        potentialBgChars.append(np.array([[A19,B19,C19],[xi19,eta19,zeta19]], dtype=np.float64))

    #test condition 20
    A20=B
    B20=A+B+zeta
    C20=sumVal+sigma
    if np.isclose(A20+B20+C20,sumVal,rtol=epsRel, atol=epsAbs, equal_nan=False):
        xi20=2*A+2*B+sigma+zeta
        eta20=2*B+xi+zeta
        zeta20=2*B+zeta
        potentialBgChars.append(np.array([[A20,B20,C20],[xi20,eta20,zeta20]], dtype=np.float64))

    #test condition 21
    A21=B
    B21=A+B-zeta
    C21=A+C-eta
    if np.isclose(A21+B21+C21,sumVal,rtol=epsRel, atol=epsAbs, equal_nan=False):
        xi21=2*A+xi-eta-zeta
        eta21=-xi+zeta
        zeta21=-2*B+zeta
        potentialBgChars.append(np.array([[A21,B21,C21],[xi21,eta21,zeta21]], dtype=np.float64))

    #test condition 22
    A22=B
    B22=A+C-eta
    C22=B+C-xi
    if np.isclose(A22+B22+C22,sumVal,rtol=epsRel, atol=epsAbs, equal_nan=False):
        xi22=2*C-xi-eta+zeta
        eta22=2*B-xi
        zeta22=-xi+zeta
        potentialBgChars.append(np.array([[A22,B22,C22],[xi22,eta22,zeta22]], dtype=np.float64))


    #test condition 23
    A23=C
    B23=A+B-zeta
    C23=A+C-eta
    if np.isclose(A23+B23+C23,sumVal,rtol=epsRel, atol=epsAbs, equal_nan=False):
        xi23=2*A+xi-eta-zeta
        eta23=-2*C+eta
        zeta23=-xi+eta
        potentialBgChars.append(np.array([[A23,B23,C23],[xi23,eta23,zeta23]], dtype=np.float64))

    #test condition 24
    A24=C
    B24=A+B-zeta
    C24=B+C-xi
    if np.isclose(A24+B24+C24,sumVal,rtol=epsRel, atol=epsAbs, equal_nan=False):
        xi24=-2*B+xi-eta+zeta
        eta24=-2*C+xi
        zeta24=-xi+eta
        potentialBgChars.append(np.array([[A24,B24,C24],[xi24,eta24,zeta24]], dtype=np.float64))


    #test condition 25
    A25=C
    B25=A+C+eta
    C25=sumVal+sigma
    if np.isclose(A25+B25+C25,sumVal,rtol=epsRel, atol=epsAbs, equal_nan=False):
        xi25=2*A+2*C+sumVal+eta
        eta25=2*C+xi+eta
        zeta25=2*C+eta
        potentialBgChars.append(np.array([[A25,B25,C25],[xi25,eta25,zeta25]], dtype=np.float64))

    return potentialBgChars


def allNormalizedBuergerCharacteristics(S):
    """

    :param S: returned from Algorithm B
    :return: all normalized Buerger characteristics belonging to the same Niggli cell
    """

    potentialBgChars=potentialBuergerCharacteristics(S)#TODO: only unique matrices should remain!
    normalizedBgCharsFlattened=[(np.around(AlgorithmN(STmp),8)).flatten() for STmp in potentialBgChars]
    #flatten the matrix to a vector to enable hashing, remove duplicated vectors using set(), then  transform the vector back to array
    uniqueData=list(map(np.array, set(map(tuple, normalizedBgCharsFlattened))))
    normalizedBgChars=[np.reshape(vec,(2,3)) for vec in uniqueData]
    return normalizedBgChars


def selectLen3(normalizedBgChars):
    """

    :param normalizedBgChars: normalized Buerger characteristics belonging to the same Niggli cell, there are 3 such matrices
    :return: the category number of the Buerger characteristics
    """
    epsRel = 1e-8
    epsAbs = 1e-6
    # if len(normalizedBgChars)!=3:
    #     raise ValueError("Invalid input: out of range.")

    S0,S1,S2=normalizedBgChars
    [[A,B,C],[_,_,_]]=S0
    #case A=B<=C
    if np.isclose(A,B, rtol=epsRel, atol=epsAbs, equal_nan=False):
        return 3
    #case A<B<C
    if A<B and B<C:
        return 11

    #case A<B=C
    #if xi is always <0
    if S0[1,0]<0 and S1[1,0]<0 and S2[1,0]<0:
        return 19

    #if there exists eta=0
    if np.isclose(S0[1,1],0, rtol=epsRel, atol=epsAbs, equal_nan=False)\
        or np.isclose(S1[1,1],0, rtol=epsRel, atol=epsAbs, equal_nan=False)\
        or np.isclose(S2[1,1],0, rtol=epsRel, atol=epsAbs, equal_nan=False):
        return 11

    #if abs(xi) is always B
    if np.isclose(np.abs(S0[1,0]),B, rtol=epsRel, atol=epsAbs, equal_nan=False)\
        and np.isclose(np.abs(S1[1,0]),B, rtol=epsRel, atol=epsAbs, equal_nan=False)\
        and np.isclose(np.abs(S2[1,0]),B, rtol=epsRel, atol=epsAbs, equal_nan=False):
        return 18

    #otherwise the category is 13
    return 13

def NiggliOfLen3(normalizedBgChars):
    """

    :param normalizedBgChars: normalized Buerger characteristics belonging to the same Niggli cell, there are 3 such matrices
    :return: S matrix of Niggli cell
    """

    category=selectLen3(normalizedBgChars)
    S0,S1,S2=normalizedBgChars
    if category==3:
        #find the matrix with the largest value of xi
        indOfXi=np.argsort([S0[1,0],S1[1,0],S2[1,0]])
        return normalizedBgChars[indOfXi[-1]]
    if category==11:
        if S0[1,0]>0:
            return S0
        elif S1[1,0]>0:
            return S1
        else:
            return S2

    if category==13:
        if S0[1,1]>0:
            return S0
        elif S1[1,1]>0:
            return S1
        else:
            return S2
    if category==18:
        indOfEta=np.argsort([S0[1,1],S1[1,1],S2[1,1]])
        return normalizedBgChars[indOfEta[-1]]

    #remaining case: category=19
    if category==19:
        indOfEta=np.argsort([S0[1,1],S1[1,1],S2[1,1]])
        return normalizedBgChars[indOfEta[0]]

    #if not returned until now, then something is wrong
    raise RuntimeError("Something is wrong for 3-Buerger-cell case.")

def selectLen2(normalizedBgChars):
    """

    :param normalizedBgChars: normalized Buerger characteristics belonging to the same Niggli cell, there are 2 such matrices
    :return: the representative category number of the Buerger characteristics
    """
    S0, S1= normalizedBgChars

    #all categories except for category 21
    if S0[1,0]>S1[1,0] or S1[1,0]>S0[1,0]:
        return 1
    else:
        return 21

def NiggliOfLen2(normalizedBgChars):
    """

    :param normalizedBgChars: normalized Buerger characteristics belonging to the same Niggli cell, there are 2 such matrices
    :return: S matrix of Niggli cell
    """
    category=selectLen2(normalizedBgChars)
    S0,S1=normalizedBgChars
    if category==1:
        if S0[1,0]>0:
            return S0
        else:
            return S1
    if category==21:
        if S0[1,1]>S1[1,1]:
            return S0
        else:
            return S1






def Buerger2Niggli(normalizedBgChars):
    """

    :param normalizedBgChars: normalized Buerger characteristics belonging to the same Niggli cell
    :return: S matrix of Niggli cell
    """

    #if there is only 1 Buerger cell, then this is the Niggli cell
    if len(normalizedBgChars)==1:
        return normalizedBgChars[0]

    #if there are 5 Buerger cells, they belong to category 20
    if len(normalizedBgChars)==5:
        candidates=[STmp for STmp in normalizedBgChars if STmp[1,0]>0]
        S0=candidates[0]
        S1=candidates[1]
        if S0[1,1]>S1[1,1]:
            return S0
        else:
            return S1

    #if there are 4 Buerger cells, they belong to category 23
    if len(normalizedBgChars)==4:
        candidates = [STmp for STmp in normalizedBgChars if STmp[1, 0] > 0]
        S0 = candidates[0]
        S1 = candidates[1]
        if S0[1,1]>S1[1,1]:
            return S0
        else:
            return S1

    #if there are 3 Buerger cells
    if len(normalizedBgChars)==3:
        return NiggliOfLen3(normalizedBgChars)
    # if there are 2 Buerger cells
    if len(normalizedBgChars)==2:
        return NiggliOfLen2(normalizedBgChars)

    raise RuntimeError("Wrong number of matrices.")



def checkNiggliCell(S):
    """

    :param S: output of Buerger2Niggli(), supposed to be a Niggli cell
    :return: check if the execution of Buerger2Niggli() is correct.
    """

    [[A,B,C],[xi,eta,zeta]]=S
    epsRel = 1e-8
    epsAbs = 1e-6
    #type-I cell
    if xi>0 and eta>0 and zeta>0:
        #main conditions
        if not (A<=B and B<=C):
            raise ValueError("Main condition A<=B<=C not satisfied.")
        if not (np.abs(xi)<=B and np.abs(eta)<=A and np.abs(zeta)<=A):
            raise ValueError("Main condition |xi|<=B, |eta|<=A, |zeta|<=A not satisfied.")
        #special conditions
        #1
        if np.isclose(A,B, rtol=epsRel, atol=epsAbs, equal_nan=False):
            if not xi<=eta:
                raise ValueError("Special condition xi<=eta in type-I cell not satisfied.")
        #2
        if np.isclose(B,C,rtol=epsRel, atol=epsAbs, equal_nan=False):
            if not eta<=zeta:
                raise ValueError("Special condition eta<=zeta in type-I cell not satisfied.")
        #3
        if np.isclose(xi,B,rtol=epsRel, atol=epsAbs, equal_nan=False):
            if not 1/2*zeta<=eta:
                raise ValueError("Special condition 1/2 zeta<=eta in type-I cell not satisfied.")

        #4
        if np.isclose(eta,A,rtol=epsRel, atol=epsAbs, equal_nan=False):
            if not 1/2*zeta<=xi:
                raise ValueError("Special condition 1/2 zeta<=xi in type-I cell not satisfied.")

        #5
        if np.isclose(zeta,A,rtol=epsRel, atol=epsAbs, equal_nan=False):
            if not 1/2*eta<=xi:
                raise ValueError("Special condition 1/2 eta<=xi in type-I cell not satisfied.")


    #type-II cell
    if xi<=0 and eta<=0 and zeta<=0:
        #main conditions
        if not (A<=B and B<=C):
            raise ValueError("Main condition A<=B<=C not satisfied.")
        if not (np.abs(xi)<=B and np.abs(eta)<=A and np.abs(zeta)<=A):
            raise ValueError("Main condition |xi|<=B, |eta|<=A, |zeta|<=A not satisfied.")
        if not (np.abs(xi)+np.abs(eta)+np.abs(zeta)<=A+B):
            raise ValueError("Main condition |xi|+|eta|+|zeta|<=A+B not satisfied.")
        #special conditions
        #1
        if np.isclose(A,B,rtol=epsRel, atol=epsAbs, equal_nan=False):
            if not np.abs(xi)<=np.abs(eta):
                raise ValueError("Special condition |xi|<=|eta| in type-II cell not satisfied.")

        #2
        if np.isclose(B,C,rtol=epsRel, atol=epsAbs, equal_nan=False):
            if not np.abs(eta)<=np.abs(zeta):
                raise ValueError("Special condition |eta|<=|zeta| in type-II cell not satisfied.")
        #3
        if np.isclose(np.abs(xi),B,rtol=epsRel, atol=epsAbs, equal_nan=False):
            if not np.isclose(zeta,0,rtol=epsRel, atol=epsAbs, equal_nan=False):
                raise ValueError("Special condition zeta=0 in type-II cell not satisfied.")
        #4
        if np.isclose(np.abs(eta),A,rtol=epsRel, atol=epsAbs, equal_nan=False):
            if not np.isclose(zeta,0,rtol=epsRel, atol=epsAbs, equal_nan=False):
                raise ValueError("Special condition zeta=0 in type-II cell not satisfied.")
        #5
        if np.isclose(np.abs(zeta),A,rtol=epsRel, atol=epsAbs, equal_nan=False):
            if not np.isclose(eta,0,rtol=epsRel, atol=epsAbs, equal_nan=False):
                raise ValueError("Special condition eta=0 in type-II cell not satisfied.")

        #6
        if np.isclose(np.abs(xi)+np.abs(eta)+np.abs(zeta),A+B,rtol=epsRel, atol=epsAbs, equal_nan=False):
            if not A<=np.abs(eta)+1/2*np.abs(zeta):
                raise ValueError("Special condition A<=|eta|+1/2|zeta| in type-II cell not satisfied.")


    #not the above 2 cells, raise error

    raise ValueError("Wrong S matrix.")


def NiggliType(S):
    """

    :param S: S matrix of a Niggli cell
    :return: type 1 or 2
    """
    epsRel = 1e-8
    epsAbs = 1e-6
    [[_,_,_],[xi,eta,zeta]]=S
    if xi>0 and eta>0 and zeta>0:
        return 1
    elif (xi<0 or np.isclose(xi,0,rtol=epsRel, atol=epsAbs, equal_nan=False))\
        and (eta<0 or np.isclose(eta,0,rtol=epsRel, atol=epsAbs, equal_nan=False))\
        and (zeta<0 or np.isclose(zeta,0,rtol=epsRel, atol=epsAbs, equal_nan=False)):
        return 2
    else:
        raise ValueError("Invalid Niggli matrix.")


def Niggli2ConventionalRow(S):
    """

    :param S: S matrix of Niggli cell returned from Buerger2Niggli()
    :return: Row in the table 9.2.5.1 in  International Tables for Crystallography Volume A_ Space-group symmetry-Springer Netherlands (2002)
    """
    [[A,B,C],[xi,eta,zeta]]=S
    D=1/2*xi
    E=1/2*eta
    F=1/2*zeta
    epsRel = 1e-8
    epsAbs = 1e-6

    #A=B=C
    if np.isclose(A,B,rtol=epsRel, atol=epsAbs, equal_nan=False)\
        and np.isclose(B,C,rtol=epsRel, atol=epsAbs, equal_nan=False):
        #1
        if NiggliType(S)==1\
            and np.isclose(D,A/2,rtol=epsRel, atol=epsAbs, equal_nan=False)\
            and np.isclose(E,A/2,rtol=epsRel, atol=epsAbs, equal_nan=False)\
            and np.isclose(F,A/2,rtol=epsRel, atol=epsAbs, equal_nan=False):
            return {"No":1, "Lattice symmetry":"Cubic", "Bravais type":"cF", "Transformation matrix":np.array([[1,-1,1],[1,1,-1],[-1,1,1]],dtype=np.float64)}
        #2
        if NiggliType(S)==1\
            and np.isclose(D,D,rtol=epsRel, atol=epsAbs, equal_nan=False)\
            and np.isclose(E,D,rtol=epsRel, atol=epsAbs, equal_nan=False)\
            and np.isclose(F,D,rtol=epsRel, atol=epsAbs, equal_nan=False):
            return {"No":2, "Lattice symmetry":"Rhombohedral", "Bravais type":"hR", "Transformation matrix":np.array([[1,-1,0],[-1,0,1],[-1,-1,-1]],dtype=np.float64)}

        #3
        if NiggliType(S)==2\
            and np.isclose(D,0,rtol=epsRel, atol=epsAbs, equal_nan=False)\
            and np.isclose(E,0,rtol=epsRel, atol=epsAbs, equal_nan=False)\
            and np.isclose(F,0,rtol=epsRel, atol=epsAbs, equal_nan=False):
            return {"No":3, "Lattice symmetry":"Cubic", "Bravais type":"cP", "Transformation matrix":np.array([[1,0,0],[0,1,0],[0,0,1]],dtype=np.float64)}

        #5
        if NiggliType(S)==2\
            and np.isclose(D,-A/3,rtol=epsRel, atol=epsAbs, equal_nan=False)\
            and np.isclose(E,-A/3,rtol=epsRel, atol=epsAbs, equal_nan=False)\
            and np.isclose(F,-A/3,rtol=epsRel, atol=epsAbs, equal_nan=False):
            return {"No":5, "Lattice symmetry":"Cubic", "Bravais type":"cI", "Transformation matrix":np.array([[1,0,1],[1,1,0],[0,1,1]],dtype=np.float64)}


        #4
        if NiggliType(S)==2\
            and np.isclose(D,D,rtol=epsRel, atol=epsAbs, equal_nan=False)\
            and np.isclose(E,D,rtol=epsRel, atol=epsAbs, equal_nan=False)\
            and np.isclose(F,D,rtol=epsRel, atol=epsAbs, equal_nan=False):
            return {"No":4, "Lattice symmetry":"Rhombohedral", "Bravais type":"hR", "Transformation matrix":np.array([[1,-1,0],[-1,0,1],[-1,-1,-1]],dtype=np.float64)}

        #6
        if NiggliType(S)==2\
            and np.isclose(D,D,rtol=epsRel, atol=epsAbs, equal_nan=False)\
            and np.isclose(E,D,rtol=epsRel, atol=epsAbs, equal_nan=False)\
            and np.isclose(F,F,rtol=epsRel, atol=epsAbs, equal_nan=False)\
            and np.isclose(2*np.abs(D+E+F),A+B,rtol=epsRel, atol=epsAbs, equal_nan=False):
            return {"No":6, "Lattice symmetry":"Tetragonal", "Bravais type":"tI", "Transformation matrix":np.array([[0,1,1],[1,0,1],[1,1,0]],dtype=np.float64)}


        #7
        if NiggliType(S)==2\
            and np.isclose(D,D,rtol=epsRel, atol=epsAbs, equal_nan=False)\
            and np.isclose(E,E,rtol=epsRel, atol=epsAbs, equal_nan=False)\
            and np.isclose(F,E,rtol=epsRel, atol=epsAbs, equal_nan=False)\
            and np.isclose(2*np.abs(D+E+F),A+B,rtol=epsRel, atol=epsAbs, equal_nan=False):
            return {"No":7, "Lattice symmetry":"Tetragonal", "Bravais type":"tI", "Transformation matrix":np.array([[1,0,1],[1,1,0],[0,1,1]],dtype=np.float64)}

        #8
        if NiggliType(S)==2\
            and np.isclose(D,D,rtol=epsRel, atol=epsAbs, equal_nan=False)\
            and np.isclose(E,E, rtol=epsRel, atol=epsAbs, equal_nan=False)\
            and np.isclose(F,F,rtol=epsRel, atol=epsAbs, equal_nan=False)\
            and np.isclose(2*np.abs(D+E+F),A+B,rtol=epsRel, atol=epsAbs, equal_nan=False):
            return {"No":8, "Lattice symmetry":"Orthorhombic", "Bravais type":"oI", "Transformation matrix":np.array([[-1,-1,0],[-1,0,-1],[0,-1,-1]],dtype=np.float64)}


    #A=B, no conditions on C
    if np.isclose(A,B,rtol=epsRel, atol=epsAbs, equal_nan=False):
        #9
        if NiggliType(S)==1\
            and np.isclose(D,A/2,rtol=epsRel, atol=epsAbs, equal_nan=False)\
            and np.isclose(E,A/2,rtol=epsRel, atol=epsAbs, equal_nan=False)\
            and np.isclose(F,A/2,rtol=epsRel, atol=epsAbs, equal_nan=False):
            return {"No": 9, "Lattice symmetry": "Rhombohedral", "Bravais type": "hR",
                    "Transformation matrix": np.array([[1,0,0], [-1,1,0], [-1-1,3]], dtype=np.float64)}

        #10
        if NiggliType(S)==1\
            and np.isclose(D,D,rtol=epsRel, atol=epsAbs, equal_nan=False)\
            and np.isclose(E,D,rtol=epsRel, atol=epsAbs, equal_nan=False)\
            and np.isclose(F,F,rtol=epsRel, atol=epsAbs, equal_nan=False):
            return {"No": 10, "Lattice symmetry": "Monoclinic", "Bravais type": "mC",
                    "Transformation matrix": np.array([[1,1,0], [1,-1,0], [0,0,-1]], dtype=np.float64)}

        #11
        if NiggliType(S)==2\
            and np.isclose(D,0,rtol=epsRel, atol=epsAbs, equal_nan=False)\
            and np.isclose(E,0,rtol=epsRel, atol=epsAbs, equal_nan=False)\
            and np.isclose(F,0,rtol=epsRel, atol=epsAbs, equal_nan=False):
            return {"No": 11, "Lattice symmetry": "Tetragonal", "Bravais type": "tP",
                    "Transformation matrix": np.array([[1,0,0], [0,1,0], [0, 0, 1]], dtype=np.float64)}

        #12
        if NiggliType(S)==2\
            and np.isclose(D,0,rtol=epsRel, atol=epsAbs, equal_nan=False)\
            and np.isclose(E,0,rtol=epsRel, atol=epsAbs, equal_nan=False)\
            and np.isclose(F,-A/2,rtol=epsRel, atol=epsAbs, equal_nan=False):
            return {"No": 12, "Lattice symmetry": "Hexagonal", "Bravais type": "hP",
                    "Transformation matrix": np.array([[1,0,0], [0,1,0], [0,0,1]], dtype=np.float64)}

        #13
        if NiggliType(S)==2\
            and np.isclose(D,0,rtol=epsRel, atol=epsAbs, equal_nan=False)\
            and np.isclose(E,0,rtol=epsRel, atol=epsAbs, equal_nan=False)\
            and np.isclose(F,F,rtol=epsRel, atol=epsAbs, equal_nan=False):
            return {"No": 13, "Lattice symmetry": "Orthorhombic", "Bravais type": "oC",
                    "Transformation matrix": np.array([[1,1,0], [-1,1,0], [0,0,1]], dtype=np.float64)}

        #15
        if NiggliType(S)==2\
            and np.isclose(D,-A/2,rtol=epsRel, atol=epsAbs, equal_nan=False)\
            and np.isclose(E,-A/2,rtol=epsRel, atol=epsAbs, equal_nan=False)\
            and np.isclose(F,0,rtol=epsRel, atol=epsAbs, equal_nan=False):
            return {"No": 15, "Lattice symmetry": "Tetragonal", "Bravais type": "tI",
                    "Transformation matrix": np.array([[1,0,0], [0,1,0], [1,1,2]], dtype=np.float64)}

        #16
        if NiggliType(S)==2\
            and np.isclose(D,D,rtol=epsRel, atol=epsAbs, equal_nan=False)\
            and np.isclose(E,D,rtol=epsRel, atol=epsAbs, equal_nan=False)\
            and np.isclose(F,F,rtol=epsRel, atol=epsAbs, equal_nan=False) \
            and np.isclose(2 * np.abs(D + E + F), A + B, rtol=epsRel, atol=epsAbs, equal_nan=False):
            return {"No": 16, "Lattice symmetry": "Orthorhombic", "Bravais type": "oF",
                    "Transformation matrix": np.array([[-1,-1,0], [1,-1,0], [1,1,2]], dtype=np.float64)}

        #14
        if NiggliType(S)==2\
            and np.isclose(D,D,rtol=epsRel, atol=epsAbs, equal_nan=False)\
            and np.isclose(E,D,rtol=epsRel, atol=epsAbs, equal_nan=False)\
            and np.isclose(F,F,rtol=epsRel, atol=epsAbs, equal_nan=False):
            return {"No": 14, "Lattice symmetry": "Monoclinic", "Bravais type": "mC",
                    "Transformation matrix": np.array([[1,1,0], [-1,1,0], [0,0,1]], dtype=np.float64)}

        #17
        if NiggliType(S)==2\
            and np.isclose(D,D,rtol=epsRel, atol=epsAbs, equal_nan=False)\
            and np.isclose(E,E,rtol=epsRel, atol=epsAbs, equal_nan=False)\
            and np.isclose(F,F,rtol=epsRel, atol=epsAbs, equal_nan=False) \
            and np.isclose(2 * np.abs(D + E + F), A + B, rtol=epsRel, atol=epsAbs, equal_nan=False):
            return {"No": 17, "Lattice symmetry": "Monoclinic", "Bravais type": "mC",
                    "Transformation matrix": np.array([[1,-1,0], [1,1,0], [-1,0,-1]], dtype=np.float64)}


    #B=C, no conditions on A
    if np.isclose(B,C,rtol=epsRel, atol=epsAbs, equal_nan=False):
        #18
        if NiggliType(S)==1\
            and np.isclose(D,A/4,rtol=epsRel, atol=epsAbs, equal_nan=False)\
            and np.isclose(E,A/2,rtol=epsRel, atol=epsAbs, equal_nan=False)\
            and np.isclose(F,A/2,rtol=epsRel, atol=epsAbs, equal_nan=False):
            return {"No": 18, "Lattice symmetry": "Tetragonal", "Bravais type": "tI",
                    "Transformation matrix": np.array([[0,-1,1], [1,-1,-1], [1,0,0]], dtype=np.float64)}
        #19
        if NiggliType(S)==1\
            and np.isclose(D,D,rtol=epsRel, atol=epsAbs, equal_nan=False)\
            and np.isclose(E,A/2,rtol=epsRel, atol=epsAbs, equal_nan=False)\
            and np.isclose(F,A/2,rtol=epsRel, atol=epsAbs, equal_nan=False):
            return {"No": 19, "Lattice symmetry": "Orthorhombic", "Bravais type": "oI",
                    "Transformation matrix": np.array([[-1,0,0], [0,-1,1], [-1,1,1]], dtype=np.float64)}

        #20
        if NiggliType(S)==1\
            and np.isclose(D,D,rtol=epsRel, atol=epsAbs, equal_nan=False)\
            and np.isclose(E,E,rtol=epsRel, atol=epsAbs, equal_nan=False)\
            and np.isclose(F,E,rtol=epsRel, atol=epsAbs, equal_nan=False):
            return {"No": 20, "Lattice symmetry": "Monoclinic", "Bravais type": "mC",
                    "Transformation matrix": np.array([[0,1,1], [0,1,-1], [-1,0,0]], dtype=np.float64)}
        #21
        if NiggliType(S)==2\
            and np.isclose(D,0,rtol=epsRel, atol=epsAbs, equal_nan=False)\
            and np.isclose(E,0,rtol=epsRel, atol=epsAbs, equal_nan=False)\
            and np.isclose(F,0,rtol=epsRel, atol=epsAbs, equal_nan=False):
            return {"No": 21, "Lattice symmetry": "Tetragonal", "Bravais type": "tP",
                    "Transformation matrix": np.array([[0,1,0], [0,0,1], [1,0,0]], dtype=np.float64)}

        #22
        if NiggliType(S)==2\
            and np.isclose(D,-B/2,rtol=epsRel, atol=epsAbs, equal_nan=False)\
            and np.isclose(E,0,rtol=epsRel, atol=epsAbs, equal_nan=False)\
            and np.isclose(F,0,rtol=epsRel, atol=epsAbs, equal_nan=False):
            return {"No": 22, "Lattice symmetry": "Hexagonal", "Bravais type": "hP",
                    "Transformation matrix": np.array([[0,1,0], [0,0,1], [1,0,0]], dtype=np.float64)}
        #23
        if NiggliType(S)==2\
            and np.isclose(D,D,rtol=epsRel, atol=epsAbs, equal_nan=False)\
            and np.isclose(E,0,rtol=epsRel, atol=epsAbs, equal_nan=False)\
            and np.isclose(F,0,rtol=epsRel, atol=epsAbs, equal_nan=False):
            return {"No": 23, "Lattice symmetry": "Orthorhombic", "Bravais type": "oC",
                    "Transformation matrix": np.array([[0,1,1], [0,-1,1], [1,0,0]], dtype=np.float64)}

        #24
        if NiggliType(S)==2\
            and np.isclose(D,D,rtol=epsRel, atol=epsAbs, equal_nan=False)\
            and np.isclose(E,-A/3,rtol=epsRel, atol=epsAbs, equal_nan=False)\
            and np.isclose(F,-A/3,rtol=epsRel, atol=epsAbs, equal_nan=False) \
            and np.isclose(2 * np.abs(D + E + F), A + B, rtol=epsRel, atol=epsAbs, equal_nan=False):
            return {"No": 24, "Lattice symmetry": "Rhombohedral", "Bravais type": "hR",
                    "Transformation matrix": np.array([[1,2,1], [0,-1,1], [1,0,0]], dtype=np.float64)}

        #25
        if NiggliType(S)==2\
            and np.isclose(D,D,rtol=epsRel, atol=epsAbs, equal_nan=False)\
            and np.isclose(E,E,rtol=epsRel, atol=epsAbs, equal_nan=False)\
            and np.isclose(F,E,rtol=epsRel, atol=epsAbs, equal_nan=False):
            return {"No": 25, "Lattice symmetry": "Monoclinic", "Bravais type": "mC",
                    "Transformation matrix": np.array([[0,1,1], [0,-1,1], [1,0,0]], dtype=np.float64)}


    #no conditions on A, B, C

    #26
    if NiggliType(S)==1\
        and np.isclose(D,A/4,rtol=epsRel, atol=epsAbs, equal_nan=False)\
        and np.isclose(E,A/2,rtol=epsRel, atol=epsAbs, equal_nan=False)\
        and np.isclose(F,A/2, rtol=epsRel, atol=epsAbs, equal_nan=False):
        return {"No": 26, "Lattice symmetry": "Orthorhombic", "Bravais type": "oF",
                "Transformation matrix": np.array([[1,0,0], [-1,2,0], [-1,0,2]], dtype=np.float64)}

    #27
    if NiggliType(S)==1\
        and np.isclose(D,D,rtol=epsRel, atol=epsAbs, equal_nan=False)\
        and np.isclose(E,A/2,rtol=epsRel, atol=epsAbs, equal_nan=False)\
        and np.isclose(F,A/2,rtol=epsRel, atol=epsAbs, equal_nan=False):
        return {"No": 27, "Lattice symmetry": "Monoclinic", "Bravais type": "mC",
                "Transformation matrix": np.array([[-1,2,0], [-1,0,0], [0,-1,1]], dtype=np.float64)}


    #28
    if NiggliType(S)==1\
        and np.isclose(D,D,rtol=epsRel, atol=epsAbs, equal_nan=False)\
        and np.isclose(E,A/2,rtol=epsRel, atol=epsAbs, equal_nan=False)\
        and np.isclose(F,2*D,rtol=epsRel, atol=epsAbs, equal_nan=False):
        return {"No": 28, "Lattice symmetry": "Monoclinic", "Bravais type": "mC",
                "Transformation matrix": np.array([[-1,0,0], [-1,0,2], [0,1,0]], dtype=np.float64)}


    #29
    if NiggliType(S)==1\
        and np.isclose(D,D,rtol=epsRel, atol=epsAbs, equal_nan=False)\
        and np.isclose(E,2*D,rtol=epsRel, atol=epsAbs, equal_nan=False)\
        and np.isclose(F,A/2,rtol=epsRel, atol=epsAbs, equal_nan=False):
        return {"No": 29, "Lattice symmetry": "Monoclinic", "Bravais type": "mC",
                "Transformation matrix": np.array([[1,0,0], [1,-2,0], [0,0,-1]], dtype=np.float64)}

    #30
    if NiggliType(S)==1\
        and np.isclose(D,B/2,rtol=epsRel, atol=epsAbs, equal_nan=False)\
        and np.isclose(E,E,rtol=epsRel, atol=epsAbs, equal_nan=False)\
        and np.isclose(F,2*E,rtol=epsRel, atol=epsAbs, equal_nan=False):
        return {"No": 30, "Lattice symmetry": "Monoclinic", "Bravais type": "mC",
                "Transformation matrix": np.array([[0,1,0], [0,1,-2], [-1,0,0]], dtype=np.float64)}

    #31
    if NiggliType(S)==1\
        and np.isclose(D,D,rtol=epsRel, atol=epsAbs, equal_nan=False)\
        and np.isclose(E,E,rtol=epsRel, atol=epsAbs, equal_nan=False)\
        and np.isclose(F,F,rtol=epsRel, atol=epsAbs, equal_nan=False):
        return {"No": 31, "Lattice symmetry": "Triclinic", "Bravais type": "aP",
                "Transformation matrix": np.array([[1,0,0], [0,1,0], [0,0,1]], dtype=np.float64)}


    #32
    if NiggliType(S)==2\
        and np.isclose(D,0,rtol=epsRel, atol=epsAbs, equal_nan=False)\
        and np.isclose(E,0,rtol=epsRel, atol=epsAbs, equal_nan=False)\
        and np.isclose(F,0,rtol=epsRel, atol=epsAbs, equal_nan=False):
        return {"No": 32, "Lattice symmetry": "Orthorhombic", "Bravais type": "oP",
                "Transformation matrix": np.array([[1,0,0], [0,1,0], [0,0,1]], dtype=np.float64)}


    #40
    if NiggliType(S)==2\
        and np.isclose(D,-B/2,rtol=epsRel, atol=epsAbs, equal_nan=False)\
        and np.isclose(E,0,rtol=epsRel, atol=epsAbs, equal_nan=False)\
        and np.isclose(F,0,rtol=epsRel, atol=epsAbs, equal_nan=False):
        return {"No": 40, "Lattice symmetry": "Orthorhombic", "Bravais type": "oC",
                "Transformation matrix": np.array([[0,-1,0], [0,1,2], [-1,0,0]], dtype=np.float64)}

    #35
    if NiggliType(S)==2\
        and np.isclose(D,D,rtol=epsRel, atol=epsAbs, equal_nan=False)\
        and np.isclose(E,0,rtol=epsRel, atol=epsAbs, equal_nan=False)\
        and np.isclose(F,0,rtol=epsRel, atol=epsAbs, equal_nan=False):
        return {"No": 35, "Lattice symmetry": "Monoclinic", "Bravais type": "mP",
                "Transformation matrix": np.array([[0,-1,0], [-1,0,0], [0,0,-1]], dtype=np.float64)}


    #36
    if NiggliType(S)==2\
        and np.isclose(D,0,rtol=epsRel, atol=epsAbs, equal_nan=False)\
        and np.isclose(E,-A/2,rtol=epsRel, atol=epsAbs, equal_nan=False)\
        and np.isclose(F,0,rtol=epsRel, atol=epsAbs, equal_nan=False):
        return {"No": 36, "Lattice symmetry": "Orthorhombic", "Bravais type": "oC",
                "Transformation matrix": np.array([[1,0,0], [-1,0,-2], [0,1,0]], dtype=np.float64)}

    #33
    if NiggliType(S)==2\
        and np.isclose(D,0,rtol=epsRel, atol=epsAbs, equal_nan=False)\
        and np.isclose(E,E,rtol=epsRel, atol=epsAbs, equal_nan=False)\
        and np.isclose(F,0,rtol=epsRel, atol=epsAbs, equal_nan=False):
        return {"No": 33, "Lattice symmetry": "Monoclinic", "Bravais type": "mP",
                "Transformation matrix": np.array([[1,0,0], [0,1,0], [0,0,1]], dtype=np.float64)}


    #38
    if NiggliType(S)==2\
        and np.isclose(D,0,rtol=epsRel, atol=epsAbs, equal_nan=False)\
        and np.isclose(E,0,rtol=epsRel, atol=epsAbs, equal_nan=False)\
        and np.isclose(F,-A/2,rtol=epsRel, atol=epsAbs, equal_nan=False):
        return {"No": 38, "Lattice symmetry": "Orthorhombic", "Bravais type": "oC",
                "Transformation matrix": np.array([[-1,0,0], [1,2,0], [0,0,-1]], dtype=np.float64)}


    #34
    if NiggliType(S)==2\
        and np.isclose(D,0,rtol=epsRel, atol=epsAbs, equal_nan=False)\
        and np.isclose(E,0,rtol=epsRel, atol=epsAbs, equal_nan=False)\
        and np.isclose(F,F,rtol=epsRel, atol=epsAbs, equal_nan=False):
        return {"No": 34, "Lattice symmetry": "Monoclinic", "Bravais type": "mP",
                "Transformation matrix": np.array([[-1,0,0], [0,0,-1], [0,-1,0]], dtype=np.float64)}

    #42
    if NiggliType(S)==2\
        and np.isclose(D,-B/2,rtol=epsRel, atol=epsAbs, equal_nan=False)\
        and np.isclose(E,-A/2,rtol=epsRel, atol=epsAbs, equal_nan=False)\
        and np.isclose(F,0,rtol=epsRel, atol=epsAbs, equal_nan=False):
        return {"No": 42, "Lattice symmetry": "Orthorhombic", "Bravais type": "oI",
                "Transformation matrix": np.array([[-1,0,0], [0,-1,0], [1,1,2]], dtype=np.float64)}


    #41
    if NiggliType(S)==2\
        and np.isclose(D,-B/2,rtol=epsRel, atol=epsAbs, equal_nan=False)\
        and np.isclose(E,E,rtol=epsRel, atol=epsAbs, equal_nan=False)\
        and np.isclose(F,0,rtol=epsRel, atol=epsAbs, equal_nan=False):
        return {"No": 41, "Lattice symmetry": "Monoclinic", "Bravais type": "mC",
                "Transformation matrix": np.array([[0,-1,-2], [0,-1,0], [-1,0,0]], dtype=np.float64)}


    #37
    if NiggliType(S)==2\
        and np.isclose(D,D,rtol=epsRel, atol=epsAbs, equal_nan=False)\
        and np.isclose(E,-A/2,rtol=epsRel, atol=epsAbs, equal_nan=False)\
        and np.isclose(F,0,rtol=epsRel, atol=epsAbs, equal_nan=False):
        return {"No": 37, "Lattice symmetry": "Monoclinic", "Bravais type": "mC",
                "Transformation matrix": np.array([[1,0,2], [1,0,0], [0,1,0]], dtype=np.float64)}



    #39
    if NiggliType(S)==2\
        and np.isclose(D,D,rtol=epsRel, atol=epsAbs, equal_nan=False)\
        and np.isclose(E,0,rtol=epsRel, atol=epsAbs, equal_nan=False)\
        and np.isclose(F,-A/2,rtol=epsRel, atol=epsAbs, equal_nan=False):
        return {"No": 39, "Lattice symmetry": "Monoclinic", "Bravais type": "mC",
                "Transformation matrix": np.array([[-1,-2,0], [-1,0,0], [0,0,-1]], dtype=np.float64)}


    #43
    if NiggliType(S)==2\
        and np.isclose(D,D,rtol=epsRel, atol=epsAbs, equal_nan=False)\
        and np.isclose(E,E,rtol=epsRel, atol=epsAbs, equal_nan=False)\
        and np.isclose(F,F,rtol=epsRel, atol=epsAbs, equal_nan=False)\
        and np.isclose(2*np.abs(D+E+F),A+B,rtol=epsRel, atol=epsAbs, equal_nan=False)\
        and np.isclose(np.abs(2*D+F),B,rtol=epsRel, atol=epsAbs, equal_nan=False):
        return {"No": 43, "Lattice symmetry": "Monoclinic", "Bravais type": "mI",
                "Transformation matrix": np.array([[-1,0,0], [-1,-1,-2], [0,-1,0]], dtype=np.float64)}


    #44
    if NiggliType(S)==2\
        and np.isclose(D,D,rtol=epsRel, atol=epsAbs, equal_nan=False)\
        and np.isclose(E,E,rtol=epsRel, atol=epsAbs, equal_nan=False)\
        and np.isclose(F,F,rtol=epsRel, atol=epsAbs, equal_nan=False):
        return {"No": 44, "Lattice symmetry": "Triclinic", "Bravais type": "aP",
                "Transformation matrix": np.array([[1,0,0], [0,1,0], [0,0,1]], dtype=np.float64)}

    #if none of above is satisfied
    raise ValueError("Conventional cell is not found.")




def getConventional(S):
    """

    :param S: S matrix of Niggli cell returned from Buerger2Niggli()
    :return: information concerning the conventional cell
    """
    convInfoList=Niggli2ConventionalRow(S)
    # No=convInfoList["No"]
    # crystalSystem=convInfoList["Lattice symmetry"]
    # BrvTyp=convInfoList["Bravais type"]
    M=convInfoList["Transformation matrix"]

    [[A,B,C],[xi,eta,zeta]]=S
    D = 1 / 2 * xi
    E = 1 / 2 * eta
    F = 1 / 2 * zeta
    m0=M[0,:]
    m1=M[1,:]
    m2=M[2,:]
    #for definition of R2, see notes
    R2=np.array([[A,F,E],
                 [F,B,D],
                 [E,D,C]])

    AConv=m0@R2@m0

    BConv=m1@R2@m1

    CConv=m2@R2@m2

    DConv=m1@R2@m2

    EConv=m2@R2@m0

    FConv=m0@R2@m1

    aConv=np.sqrt(AConv)
    bConv=np.sqrt(BConv)
    cConv=np.sqrt(CConv)

    alphaConv=np.arccos(DConv/(bConv*cConv))

    betaConv=np.arccos(EConv/(aConv*cConv))

    gammaConv=np.arccos(FConv/(aConv*bConv))

    convParamsAndInfo={"a":aConv,"b":bConv,"c":cConv,"alpha":alphaConv,"beta":betaConv,"gamma":gammaConv,
                     "Crystal System":convInfoList["Lattice symmetry"],"Bravais type": convInfoList["Bravais type"] }
    #the angles are in radians
    return convParamsAndInfo



def getConventionalInDegree(S):
    """

    :param S: S matrix of Niggli cell returned from Buerger2Niggli()
    :return: information concerning the conventional cell, with angles in degrees
    """
    convParamsAndInfo=getConventional(S)
    alphaInRadian=convParamsAndInfo["alpha"]
    betaInRadian=convParamsAndInfo["beta"]
    gammaInRadian=convParamsAndInfo["gamma"]

    alphaInDegree=math.degrees(alphaInRadian)
    betaInDegree = math.degrees(betaInRadian)
    gammaInDegree = math.degrees(gammaInRadian)

    convParamsAndInfo["alpha"]=alphaInDegree
    convParamsAndInfo["beta"]=betaInDegree
    convParamsAndInfo["gamma"]=gammaInDegree

    return convParamsAndInfo




def prim2conv(a,b,c,alpha,beta,gamma):
    """

    :param a: length
    :param b: length
    :param c: length
    :param alpha: angle between b and c, in degrees
    :param beta: angle between c and a, in degrees
    :param gamma: angle between a and b, in degrees
    :return: information about the conventional cell, with angles in degrees
    """
    #input parameters (angle in degrees) to matrix S
    S=inputParametersInDegree2S(a,b,c,alpha,beta,gamma)

    #find Buerger cell
    SofB=AlgorithmB(S)

    #find all Buerger cells
    normalizedBgChars=allNormalizedBuergerCharacteristics(SofB)

    #find Niggli cell
    SofN=Buerger2Niggli(normalizedBgChars)

    #Niggli cell to conventional cell (angle in degrees)
    convParamsAndInfo=getConventionalInDegree(SofN)

    return convParamsAndInfo





