import numpy as np
import math


# this script transforms a primitive cell to a conventional cell


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
    [[A, B, C], [xi, eta, zeta]] = S
    a = np.sqrt(A)
    b = np.sqrt(B)
    c = np.sqrt(C)
    cosAlpha = xi / (2 * b * c)
    cosBeta = eta / (2 * a * c)
    cosGamma = zeta / (2 * a * b)
    # print("Before N cos(alpha)=" + str(cosAlpha) + ", cos(beta)=" + str(cosBeta) + ", cos(gamma)=" + str(cosGamma))

    epsRel = 1e-8
    epsAbs = 1e-5
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
    [[A, B, C], [xi, eta, zeta]] = S
    a = np.sqrt(A)
    b = np.sqrt(B)
    c = np.sqrt(C)
    cosAlpha = xi / (2 * b * c)
    cosBeta = eta / (2 * a * c)
    cosGamma = zeta / (2 * a * b)
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
    [[A, B, C], [xi, eta, zeta]] = S
    a=np.sqrt(A)
    b=np.sqrt(B)
    c=np.sqrt(C)
    cosAlpha=xi/(2*b*c)
    cosBeta=eta/(2*a*c)
    cosGamma=zeta/(2*a*b)
    # print("Before B3: a="+str(a)+", b="+str(b)+", c="+str(c)+", xi="+str(xi)+", eta="+str(eta)+", zeta="+str(zeta) +", cos(alpha)="+str(cosAlpha)+", cos(beta)="+str(cosBeta)+", cos(gamma)="+str(cosGamma))

    j = math.floor((eta + A) / (2 * A))
    # CNew=c**2+j**2*a**2-j*2*a*c*cosBeta

    C = C + j ** 2 * A - j * eta
    # print("In B3: CNew=" + str(CNew) + ", C=" + str(C))
    xi = xi - j * zeta
    eta = eta - 2 * j * A

    S[0, 2] = C
    S[1, 0] = xi
    S[1, 1] = eta

    [[A, B, C], [xi, eta, zeta]] = S
    a = np.sqrt(A)
    b = np.sqrt(B)
    c = np.sqrt(C)
    cosAlpha = xi / (2 * b * c)
    cosBeta = eta / (2 * a * c)
    cosGamma = zeta / (2 * a * b)
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
        [[A, B, _], [xi, _, _]] = S
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



