import numpy as np
import math
#this script transforms a primitive cell to a conventional cell



def inputParameters2S(a,b,c,alpha,beta,gamma):
    """

    :param a: length
    :param b: length
    :param c: length
    :param alpha: angle between b and c, in degrees
    :param beta: angle between c and a, in degrees
    :param gamma: angle between a and b, in degrees
    :return: matrix S
    """


    #sanity check
    if a<=0 or b<=0 or c<=0 or alpha<=0 or alpha>=180 or beta<=0 or beta>=180 or gamma<=0 or gamma>=180:
        raise ValueError("Invalid input.")


    #degree to radian
    alphaInRadian = math.radians(alpha)
    print(np.cos(alphaInRadian))
    betaInRadians = math.radians(beta)
    gammaInRadians = math.radians(gamma)

    S=np.array([[a**2,b**2,c**2],
                [2*b*c*np.cos(alphaInRadian),2*a*c*np.cos(betaInRadians),2*a*b*np.cos(gammaInRadians)]],dtype=np.float64)

    return S


def AlgorithmN(S):
    """

    :param S: matrix S
    :return:
    """