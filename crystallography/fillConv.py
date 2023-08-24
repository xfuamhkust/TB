import numpy as np
# from prim2conv import prim2convWithBasis
import copy


# This script computes positions of all atoms within a conventional cell


def primitivePoints(origin, a, b, c):
    """

    :param origin: origin of the conventional cell
    :param a: a vector of conventional cell
    :param b: b vector of conventional cell
    :param c: c vector of conventional cell
    :return: points of primitive lattice
    """
    v0 = copy.copy(origin)
    v1 = origin + a
    v2 = origin + a + b
    v3 = origin + b

    v4 = v0 + c
    v5 = v1 + c
    v6 = v2 + c
    v7 = v3 + c

    return [v0, v1, v2, v3, v4, v5, v6, v7]


def faceCenteredPoints(origin, a, b, c):
    """

    :param origin: origin of the conventional cell
    :param a: a vector of conventional cell
    :param b: b vector of conventional cell
    :param c: c vector of conventional cell
    :return: points at each center face
    """

    vec1 = origin + 1 / 2 * a + 1 / 2 * b
    vec2 = origin + 1 / 2 * a + 1 / 2 * c
    vec3 = origin + 1 / 2 * b + 1 / 2 * c

    vec4 = vec1 + c
    vec5 = vec2 + b
    vec6 = vec3 + a

    return [vec1, vec2, vec3, vec4, vec5, vec6]


def baseCenteredPoints(origin, a, b, c):
    """

        :param origin: origin of the conventional cell
        :param a: a vector of conventional cell
        :param b: b vector of conventional cell
        :param c: c vector of conventional cell
        :return: points at each  base center face
        """
    vec1 = origin + 1 / 2 * a + 1 / 2 * b
    vec2 = vec1 + c
    return [vec1, vec2]


def bodyCenteredPoint(origin, a, b, c):
    """

            :param origin: origin of the conventional cell
            :param a: a vector of conventional cell
            :param b: b vector of conventional cell
            :param c: c vector of conventional cell
            :return: point at body center
            """
    vec = origin + 1 / 2 * a + 1 / 2 * b + 1 / 2 * c
    return [vec]


def pointscF(origin, a, b, c):
    """

    :param origin: origin of the conventional cell cF
    :param a: a vector of conventional cell
    :param b: b vector of conventional cell
    :param c: c vector of conventional cell
    :return: lattice  points of cF
    """

    # TODO: check condition on a, b, c to be an cF

    # 8 vertices
    # v0=copy.copy(origin)
    # v1=origin+a
    # v2=origin+a+b
    # v3=origin+b
    #
    # v4=v0+c
    # v5=v1+c
    # v6=v2+c
    # v7=v3+c

    # 6 points on faces
    # v8=origin+1/2*a+1/2*b
    # v9=origin+1/2*a+1/2*c
    # v10=origin+1/2*b+1/2*c
    #
    # v11=v8+c
    # v12=v9+b
    # v13=v10+a

    # return [v0,v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12,v13]
    return primitivePoints(origin, a, b, c) + faceCenteredPoints(origin, a, b, c)


def pointshP(origin, a, b, c):
    """

    :param origin: origin of the conventional cell hP
    :param a: a vector of conventional cell
    :param b: b vector of conventional cell
    :param c: c vector of conventional cell
    :return: lattice  points of hP
    """
    return primitivePoints(origin, a, b, c)


def pointscP(origin, a, b, c):
    """

    :param origin: origin of the conventional cell cP
    :param a: a vector of conventional cell
    :param b: b vector of conventional cell
    :param c: c vector of conventional cell
    :return: lattice  points of cP
    """

    return primitivePoints(origin, a, b, c)


def pointsaP(origin, a, b, c):
    """

        :param origin: origin of the conventional cell aP
        :param a: a vector of conventional cell
        :param b: b vector of conventional cell
        :param c: c vector of conventional cell
        :return: lattice  points of aP
        """
    return primitivePoints(origin, a, b, c)


def pointsmP(origin, a, b, c):
    """

        :param origin: origin of the conventional cell mP
        :param a: a vector of conventional cell
        :param b: b vector of conventional cell
        :param c: c vector of conventional cell
        :return: lattice  points of mP
        """
    return primitivePoints(origin, a, b, c)


def pointsmS(origin, a, b, c):
    """

        :param origin: origin of the conventional cell mS
        :param a: a vector of conventional cell
        :param b: b vector of conventional cell
        :param c: c vector of conventional cell
        :return: lattice  points of mS
        """
    return primitivePoints(origin, a, b, c) + baseCenteredPoints(origin, a, b, c)


def pointsoP(origin, a, b, c):
    """

        :param origin: origin of the conventional cell oP
        :param a: a vector of conventional cell
        :param b: b vector of conventional cell
        :param c: c vector of conventional cell
        :return: lattice  points of oP
        """
    return primitivePoints(origin, a, b, c)


def pointsoS(origin, a, b, c):
    """

        :param origin: origin of the conventional cell oS
        :param a: a vector of conventional cell
        :param b: b vector of conventional cell
        :param c: c vector of conventional cell
        :return: lattice  points of oS
        """
    return primitivePoints(origin, a, b, c) + baseCenteredPoints(origin, a, b, c)


def pointsoI(origin, a, b, c):
    """

        :param origin: origin of the conventional cell oI
        :param a: a vector of conventional cell
        :param b: b vector of conventional cell
        :param c: c vector of conventional cell
        :return: lattice  points of oI
        """
    return primitivePoints(origin, a, b, c) + bodyCenteredPoint(origin, a, b, c)


def pointsoF(origin, a, b, c):
    """

           :param origin: origin of the conventional cell oF
           :param a: a vector of conventional cell
           :param b: b vector of conventional cell
           :param c: c vector of conventional cell
           :return: lattice  points of oF
           """
    return primitivePoints(origin, a, b, c) + faceCenteredPoints(origin, a, b, c)


def pointstP(origin, a, b, c):
    """

           :param origin: origin of the conventional cell tP
           :param a: a vector of conventional cell
           :param b: b vector of conventional cell
           :param c: c vector of conventional cell
           :return: lattice  points of tP
           """
    return primitivePoints(origin,a,b,c)

def pointstI(origin, a, b, c):
    """

           :param origin: origin of the conventional cell tI
           :param a: a vector of conventional cell
           :param b: b vector of conventional cell
           :param c: c vector of conventional cell
           :return: lattice  points of tI
           """
    return primitivePoints(origin,a,b,c)+bodyCenteredPoint(origin,a,b,c)

def pointshR(origin, a, b, c):
    """

           :param origin: origin of the conventional cell hR
           :param a: a vector of conventional cell
           :param b: b vector of conventional cell
           :param c: c vector of conventional cell
           :return: lattice  points of hR
           """
    inPoints=[origin+2/3*a+1/3*b+1/3*c,origin+1/3*a+2/3*b+2/3*c]#TODO: double check
    return primitivePoints(origin,a,b,c)+inPoints

def pointscI(origin, a, b, c):
    """

           :param origin: origin of the conventional cell cI
           :param a: a vector of conventional cell
           :param b: b vector of conventional cell
           :param c: c vector of conventional cell
           :return: lattice  points of cI
           """
    return primitivePoints(origin,a,b,c)+bodyCenteredPoint(origin,a,b,c)

def transformationOfVector(M, a, b, c):
    """

    :param M: transformation matrix
    :param a: a vector
    :param b: b vector
    :param c: c vector
    :return: vectors trasnformed by M
    """
    [[M00, M01, M02], [M10, M11, M12], [M20, M21, M22]] = M

    aVec = M00 * a + M01 * b + M02 * c
    bVec = M10 * a + M11 * b + M12 * c
    cVec = M20 * a + M21 * b + M22 * c

    return aVec, bVec, cVec


def checkMatMinDist(arr, arrList):
    """
        :param mat: an np.array
        :param matList: a list of np.array
        :return: min distance between arr and arrs in arrList
    """
    distList = [np.linalg.norm(arr - elem, ord=2) for elem in arrList]
    return np.min(distList)


def removeDuplicatedArrs(arrList):
    """

    :param matList: a list of np.array
    :return: unique arrs in the list
    """
    eps = 1e-6
    if len(arrList) == 0:
        return []
    ret = [arrList[0]]
    for elem in arrList:
        if checkMatMinDist(elem, ret) <= eps:
            continue
        else:
            ret.append(elem)
    return ret


def convCellPoints(rfcPosition, a, b, c, BrvType):
    """

    :param rfcPosition: reference point of the conventional cell where an atom is occupied
    #TODO: currently rfcPosition is the coordinates under the basis of a permutation of a,b,c. We assume a,b,c is the correct order but this is not guaranteed.
    :param a: a vector of conventional cell
    :param b: b vector of conventional cell
    :param c: c vector of conventional cell
    :param BrvType: Bravais lattice type
    :return: the lattice points for conventional cell with corner point rfcPosition
    """

    generatingFunc = BrvType2Points(BrvType)
    a = np.array(a)
    b = np.array(b)
    c = np.array(c)

    rfcCartesian = rfcPosition[0] * a + rfcPosition[1] * b + rfcPosition[2] * c

    V0 = generatingFunc(rfcCartesian, a, b, c)
    # V1=[vec-a for vec in V0]
    # V2=[vec-b for vec in V0]
    # V3=[vec-a-b for vec in V0]
    # V4=[vec-c for vec in V0]
    # V5=[vec-a-c for vec in V0]
    # V6=[vec-b-c for vec in V0]
    # V7=[vec-a-b-c for vec in V0]

    # U=V0#+V1+V2+V3+V4+V5+V6+V7
    #
    # uniquePoints=removeDuplicatedArrs(U)

    return V0


def withinRange(a, x, y):
    """

    :param a:
    :param x:
    :return: whether a in [x,y]
    """
    epsRel = 1e-8
    epsAbs = 1e-6

    rst = (np.isclose(a, x, rtol=epsRel, atol=epsAbs, equal_nan=False) or x < a) \
          and (np.isclose(a, y, rtol=epsRel, atol=epsAbs, equal_nan=False) or a < y)

    return rst


def truncatedPointsInConventionalCell(atmPosition, a, b, c, BrvType):
    """

    :param atmPosition: atom position within the primitive cell
        #TODO: currently rfcPosition is the coordinates under the basis of a permutation of a,b,c. We assume a,b,c is the correct order but this is not guaranteed.

    :param a: a vector of conventional cell
    :param b: b vector of conventional cell
    :param c: c vector of conventional cell
    :param BrvType: BrvType: Bravais lattice type
    :return:  positions (modulo a, b, c) of atmPosition within the conventional cell with reference point rfcPosition
    """

    convAtmPoints = convCellPoints(atmPosition, a, b, c, BrvType)  # under Cartesian coordinates
    # print("conv points:"+str(convAtmPoints))

    aLen = np.linalg.norm(a, ord=2)
    bLen = np.linalg.norm(b, ord=2)
    cLen = np.linalg.norm(c, ord=2)

    basis = np.array([a, b, c], dtype=np.float64).T

    basisInv = np.linalg.inv(basis)

    # to compute the atoms within the parallelepiped origin+a, origin+b, origin+c
    # coordinates of atoms under a,b,c basis of conventional cell and with origin
    convAtmUnderBasis = [basisInv @ pnt for pnt in convAtmPoints]

    # translation into the parallelepiped
    convAtmTranslated = []
    for pnt in convAtmUnderBasis:
        x = pnt[0] % aLen
        y = pnt[1] % bLen
        z = pnt[2] % cLen
        convAtmTranslated.append(np.array([x, y, z], dtype=np.float64))

    uniquePoints = removeDuplicatedArrs(convAtmTranslated)

    # transform to original basis
    ret = [basis @ pnt for pnt in uniquePoints]

    return ret


def BrvType2Points(BrvType):
    """

    :param BrvType: type of Bravais lattice
    :return: function that computes all lattice points for BrvType
    """
    #triclinlic
    if BrvType=="aP":
        return pointsaP
    #monoclinic
    if BrvType=="mP":
        return pointsmP
    if BrvType=="mS":
        return pointsmS
    #orthorhombic
    if BrvType=="oP":
        return pointsoP
    if BrvType=="oS":
        return pointsoS
    if BrvType=="oI":
        return pointsoI
    if BrvType=="oF":
        return pointsoF

    #tetragonal
    if BrvType=="tP":
        return pointstP
    if BrvType=="tI":
        return pointstI
    #rhombohedral
    if BrvType=="hR":
        return pointshR
    #hexagonal
    if BrvType == "hP":
        return pointshP
    #cubic
    if BrvType == "cF":
        return pointscF

    if BrvType == "cP":
        return pointscP
    if BrvType=="cI":
        return pointscI

    raise ValueError("Invalid Bravais type.")

# test data

# NaCl
# a=np.array([ 0.00000000 ,    0.50000000  ,   0.50000000])
# b=np.array([ 0.50000000 ,    0.00000000 ,    0.50000000])
# c=np.array([ 0.50000000,     0.50000000 ,    0.00000000])
# convParamsAndInfo=prim2convWithBasis(a,b,c)
# origin=np.array([0,0,0])
# Na0=np.array([ 0.00000000,     0.00000000,     0.00000000])
# Cl0=np.array([ 0.50000000,     0.50000000,     0.50000000])
# aBasis=convParamsAndInfo["Basis a"]
# bBasis=convParamsAndInfo["Basis b"]
# cBasis=convParamsAndInfo["Basis c"]
# brvType=convParamsAndInfo["Bravais type"]
# ltPointsNa0=truncatedPointsInConventionalCell(origin,Na0,aBasis,bBasis,cBasis,brvType)
# ltPointsCl0=truncatedPointsInConventionalCell(origin,Cl0,aBasis,bBasis,cBasis,brvType)
# print(ltPointsNa0)
# Si
# a=np.array([ 0.00000000,     1/2,     1/2])
# b=np.array([ 1/2,     0.00000000 ,    1/2])
# c=np.array([ 1/2,     1/2,     0.00000000])
# #
# convParamsAndInfo=prim2convWithBasis(a,b,c)
# si0=np.array([0,    0,    0])
# si1=np.array([ 1/4,     1/4,     1/4])
# aBasis=convParamsAndInfo["Basis a"]
# bBasis=convParamsAndInfo["Basis b"]
# cBasis=convParamsAndInfo["Basis c"]
# #
# origin=np.array([0,0,0])
# brvType=convParamsAndInfo["Bravais type"]
# ltPointsSi0=truncatedPointsInConventionalCell(origin,si0,aBasis,bBasis,cBasis,brvType)
# ltPointsSi1=truncatedPointsInConventionalCell(origin,si1,aBasis,bBasis,cBasis,brvType)
# print(ltPointsSi1)

# graphene
# a=np.array([ 0.50000000,    -0.86602540,     0.00000000])
# b=np.array([ 0.50000000,     0.86602540,     0.00000000])
# c=np.array([ 0.00000000,     0.00000000 ,    1.00000000])
# convParamsAndInfo=prim2convWithBasis(a,b,c)
# #
# C0=np.array([ 0.33333333 ,    0.66666666,     0.00000000])
# C1=np.array([ 0.66666666 ,    0.33333333  ,   0.00000000])
# aBasis=convParamsAndInfo["Basis a"]
# bBasis=convParamsAndInfo["Basis b"]
# cBasis=convParamsAndInfo["Basis c"]
# #
# brvType=convParamsAndInfo["Bravais type"]
# origin=np.array([0,0,0])
# ltPointsC0=truncatedPointsInConventionalCell(origin,C0,aBasis,bBasis,cBasis,brvType)
# ltPointsC1=truncatedPointsInConventionalCell(origin,C1,aBasis,bBasis,cBasis,brvType)
# print(ltPointsC1)
