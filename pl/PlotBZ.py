import numpy as np
import matplotlib.pyplot as plt   
pi2 = 2*np.pi              
sr3 = np.sqrt(3)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   
plt.rcParams['figure.dpi']= 300

def PlotBrillouinZone(Dim, BraLat,LvC=[]):
    if Dim == 3:
        return PlotBrillouinZone3D(BraLat,LvC)
    elif Dim == 2:
        return PlotBrillouinZone2D(BraLat,LvC)
    
def PlotBrillouinZone2D(BraLat,LvC=[]):
    
    '''
    This function is to plot the high-symmetry points in 1st Brillouin Zone (2D).
    Input: 
        BraLat: Bravais lattice. tP, hP, oP, oC, mP
        LvC: Conventional lattice vectors. 3*3, else invalid.
    Intermediate:
        See PlotBrillouinZone3D
    Ref: https://wiki.fysik.dtu.dk/ase/ase/dft/kpoints.html#high-symmetry-paths
    '''
    
    # Get a, b, gamma from input LvC or defalt
    if np.array(LvC).shape == (3,3):
        a = np.linalg.norm(LvC[0])
        b = np.linalg.norm(LvC[1])
        gamma = np.arccos(LvC[0]@LvC[1]/a/b)
    else:
        a = 1.
        b = 1.2
        gamma = 76. # only influences mp
    LvC_ab = [a, b, gamma]
    
    # Get LvP, LvK, LatName
    LvP, LvK, LatName = get_LvP_LatName_2d(BraLat,LvC_ab)
    
    # Get HiSymKptName & HiSymKptLv & HiSymLine
    HiSymKptName, HiSymKptLv, HiSymLine = get_high_sym_2d(LatName,LvP,LvC_ab)
    HiSymKpt = HiSymKptLv@LvK
    NumHiSymKpt = len(HiSymKpt)
    
    # Get BZ edges
    Vtx, Edge = get_brillouin_zone_2d(LvK[:2,:2])
    
    # Plot Reciprocal primitive lattice vectors, BZ edges, high symmetry points
    fig = plt.figure()
    plt.axis("equal")
    plt.axis("off")
    rmax = np.max(np.linalg.norm(LvK,axis=0))
    plt.xlim([-2*rmax, 2*rmax])
    plt.ylim([-2*rmax, 2*rmax])
    # Plot reciprocal primitive lattice vectors
    for i in range(2):
        x0, y0, _ = 1/2*LvK[i]
        u , v , _ = 1/4*LvK[i]
        plt.quiver(x0, y0, u, v, color = "k",linewidth=0.5)
        plt.text(x0+u,y0+v,"$b_%d$"%(i+1),size=7)
    # Plot BZ edges
    for Edgei in Edge:
        plt.plot(Edgei[:,0],Edgei[:,1],"k-",lw=0.8)
    # Plot high symmetry points
    for i in range(NumHiSymKpt):
        plt.plot(HiSymKpt[i,0],HiSymKpt[i,1],"ro")
        plt.text(HiSymKpt[i,0],HiSymKpt[i,1],HiSymKptName[i])
    # Plot high symmetry lines
    for HiSymLine0 in HiSymLine:
        x1, y1, _ = HiSymKpt[FindName(HiSymKptName,HiSymLine0[0])]
        for i in range(1,len(HiSymLine0)):
            x2, y2, _ = HiSymKpt[FindName(HiSymKptName,HiSymLine0[i])]
            plt.plot([x1,x2],[y1,y2],"r-")
            x1, y1 = x2, y2
    # Output
    ParaBZPlt = {"FigureBZ": fig,
                  }
    
    return ParaBZPlt

def PlotBrillouinZone3D(BraLat,LvC=[]):
    
    '''
    This function is to plot the high-symmetry points in 1st Brillouin Zone (3D).
    Input: 
        BraLat: Bravais lattice. cP, cF, cI, tP, tI, oP, oF, oI, oS, hP, hR, mP, mS, aP
        LvC: Conventional lattice vectors. 3*3, else invalid.
    Intermediate:
        LvP: Primitive lattice vectors.
        LvK: Reciprocal primitive lattice vectors.
        LatName: Lattice name. CUB, FCC, BCC, TET, BCT, ORC, ORCF(1,2,3), ORCI, ORCC, HEX, RHL(1,2), MCL, MCLC(1,2,3,4,5), TRI(1a,2a,1b,2b), 
        HiSymKptName: The Names of high symmetry k-points
        HiSymKptLv: High symmetry k-points in basis of LvK
        HiSymKpt: High symmetry k-points in basis of [kx,ky,kz]
    Reference:
        https://doi.org/10.1016/j.commatsci.2010.05.010
            a,  b,  c,  alpha,  beta,  gamma:   Conventional lattice vectors
            ak, bk, ck, alphak, betak, gammak:  Reciprocal primitive lattice vectors
        https://eng.libretexts.org/Bookshelves/Materials_Science/TLP_Library_I/15%3A_Crystallography/15.04%3A_Section_4-#:~:text=The%20angles%20between%20the%20crystallographic%20axes%20are%20defined,called%20%E2%80%98unit%20cell%20parameters%E2%80%99%2C%20or%20just%20%E2%80%98cell%20parameters%E2%80%99%29.
            α = the angle between b and c
            β = the angle between a and c
            γ = the angle between a and b
    '''
    
    # Get a, b, c, alpha, beta, gamma from input LvC or defalt
    if np.array(LvC).shape == (3,3):
        a = np.linalg.norm(LvC[0])
        b = np.linalg.norm(LvC[1])
        c = np.linalg.norm(LvC[2])
        alpha = np.arccos(LvC[1]@LvC[2]/b/c)
        beta  = np.arccos(LvC[2]@LvC[0]/c/a)
        gamma = np.arccos(LvC[0]@LvC[1]/a/b)
    else:
        a = 1.
        b = 1.2
        c = 1.5
        alpha = 75.
        beta  = 85.
        gamma = 95.
        # a = 1/np.sqrt(1/b**2 + 1/c**2)-0.1
        # a = 0.5
        # b = 0.9
        # c = 0.8
    LvC_abc = [a, b, c, alpha, beta, gamma]
    
    # Get LvP, LvK, LatName
    LvP, LvK, LatName = get_LvP_LatName_3d(BraLat,LvC_abc)
    
    # Get HiSymKptName & HiSymKptLv & HiSymLine
    HiSymKptName, HiSymKptLv, HiSymLine = get_high_sym_3d(LatName,LvP,LvC_abc)
    HiSymKpt = HiSymKptLv@LvK
    NumHiSymKpt = len(HiSymKpt)
    
    # Get BZ edges
    Vtx, Edge, Face = get_brillouin_zone_3d(LvK)
    
    # Plot Reciprocal primitive lattice vectors, BZ edges, high symmetry points
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d', proj_type='ortho')
    ax.view_init(elev=6, azim=25)
    ax.set_box_aspect([1, 1, 1])
    ax.axis("off")
    # Plot reciprocal primitive lattice vectors
    for i in range(3):
        x0, y0, z0 = 1/2*LvK[i]
        u , v , w  = 1/5*LvK[i]
        ax.quiver3D(x0, y0, z0, u, v, w, length=1, arrow_length_ratio=.1,color = "k",linewidth=0.5)
        ax.text(x0+u,y0+v,z0+w,"$b_%d$"%(i+1),size=7)
    # Plot BZ edges
    for Edgei in Edge:
        ax.plot(Edgei[:,0],Edgei[:,1],Edgei[:,2],"k-",lw=0.8)
    # Plot high symmetry points
    for i in range(NumHiSymKpt):
        ax.plot(HiSymKpt[i,0],HiSymKpt[i,1],HiSymKpt[i,2],"ro")
        ax.text(HiSymKpt[i,0],HiSymKpt[i,1],HiSymKpt[i,2],HiSymKptName[i])
    # Plot high symmetry lines
    for HiSymLine0 in HiSymLine:
        x1, y1, z1 = HiSymKpt[FindName(HiSymKptName,HiSymLine0[0])]
        for i in range(1,len(HiSymLine0)):
            x2, y2, z2 = HiSymKpt[FindName(HiSymKptName,HiSymLine0[i])]
            ax.plot([x1,x2],[y1,y2],[z1,z2],"r-")
            x1, y1, z1 = x2, y2, z2
    # Output
    ParaBZPlt = {"FigureBZ": fig,
                  }
    
    return ParaBZPlt

def GetRecLv(Lv):
    a1,a2,a3 = Lv
    Vol = np.dot(a1,np.cross(a2,a3))
    b1 = pi2 * np.cross(a2,a3) / Vol
    b2 = pi2 * np.cross(a3,a1) / Vol
    b3 = pi2 * np.cross(a1,a2) / Vol
    RecLv = np.array([b1,b2,b3])
    return RecLv

def FindName(Kname,Kname0):
    for i in range(len(Kname)):
        if Kname0 == Kname[i]:
            return i
    return None

def get_LvP_LatName_2d(BraLat,LvC_ab,error = 1e-3):
    a, b, gamma = LvC_ab
    if   BraLat == "tP":
        LatName = "SQR"
        LvP = np.array([[ a, 0, 0],
                        [ 0, a, 0],
                        [ 0, 0, 1]])
    elif BraLat == "hP":
        LatName = "HEX2D"
        LvP = np.array([[ 1/2*a,-sr3/2*a, 0],
                        [ 1/2*a, sr3/2*a, 0],
                        [     0,       0, 1]])
    elif BraLat == "oP":
        LatName = "RECT"
        LvP = np.array([[ a, 0, 0],
                        [ 0, b, 0],
                        [ 0, 0, 1]])
    elif BraLat == "oC":
        # a < b
        LatName = "CRECT"
        LvP = np.array([[ a/2, b/2, 0],
                        [-a/2, b/2, 0],
                        [   0,   0, 1]])
    elif BraLat == "mP":
        # a < b, gamma < 90
        LatName = "OBL"
        cg  = np.cos(gamma/180*np.pi)
        sg  = np.sin(gamma/180*np.pi)
        LvP = np.array([[    a,    0, 0],
                        [ b*cg, b*sg, 0],
                        [    0,    0, 1]])
    LvK = GetRecLv(LvP)
    return LvP, LvK, LatName

def get_LvP_LatName_3d(BraLat,LvC_abc,error = 1e-3):
    a, b, c, alpha, beta, gamma = LvC_abc
    ca  = np.cos(alpha/180*np.pi)
    sa  = np.sin(alpha/180*np.pi)
    cb  = np.cos(beta /180*np.pi)
    cg  = np.cos(gamma/180*np.pi)
    sg  = np.sin(gamma/180*np.pi)
    sa2 = np.sin(alpha/180*np.pi/2)
    if BraLat[0] == "c":
        if   BraLat == "cP":
            LatName = "CUB"
            LvP = np.array([[ a, 0, 0],
                            [ 0, a, 0],
                            [ 0, 0, a]])
        elif BraLat == "cF":
            LatName = "FCC"
            LvP = np.array([[   0, a/2, a/2],
                            [ a/2,   0, a/2],
                            [ a/2, a/2,   0]])
        elif BraLat == "cI":
            LatName = "BCC"
            LvP = np.array([[-a/2, a/2, a/2],
                            [ a/2,-a/2, a/2],
                            [ a/2, a/2,-a/2]])
        LvK = GetRecLv(LvP)
    elif BraLat[0] == "t":
        if   BraLat == "tP":
            LatName = "TET"
            LvP = np.array([[ a, 0, 0],
                            [ 0, a, 0],
                            [ 0, 0, c]])
        elif BraLat == "tI":
            if   c < a:
                LatName = "BCT1"
            elif c > a:
                LatName = "BCT2"
            LvP = np.array([[-a/2, a/2, c/2],
                            [ a/2,-a/2, c/2],
                            [ a/2, a/2,-c/2]])
        LvK = GetRecLv(LvP)
    elif BraLat[0] == "o":
        if   BraLat == "oP":
            # a < b < c
            LatName = "ORC"
            LvP = np.array([[ a, 0, 0],
                            [ 0, b, 0],
                            [ 0, 0, c]])
        elif BraLat == "oF":
            # a < b < c
            x = 1/a**2 - (1/b**2 + 1/c**2)
            if abs(x)/(1/a**2 + 1/b**2 + 1/c**2) < error:
                LatName = "ORCF3"
            elif x > 0:
                LatName = "ORCF1"
            elif x < 0:
                LatName = "ORCF2"
            LvP = np.array([[   0, b/2, c/2],
                            [ a/2,   0, c/2],
                            [ a/2, b/2,   0]])
        elif BraLat == "oI":
            # a < b < c
            LatName = "ORCI"
            LvP = np.array([[-a/2, b/2, c/2],
                            [ a/2,-b/2, c/2],
                            [ a/2, b/2,-c/2]])
        elif BraLat == "oS":
            # a < b
            LatName = "ORCC"
            LvP = np.array([[ a/2,-b/2, 0],
                            [ a/2, b/2, 0],
                            [   0,   0, c]])
        LvK = GetRecLv(LvP)
    elif BraLat == "hP":
        LatName = "HEX"
        LvP = np.array([[ a/2,-a*sr3/2, 0],
                        [ a/2, a*sr3/2, 0],
                        [   0,       0, c]])
        LvK = GetRecLv(LvP)
    elif BraLat == "hR":
        if   alpha < 90:
            LatName = "RHL1"
        elif alpha > 90:
            LatName = "RHL2"
        # LvP = np.array([[ a*ca2,   -a*sa2,                         0],
        #                 [ a*ca2,    a*sa2,                         0],
        #                 [ a*ca/ca2,     0, a*np.sqrt(1-ca**2/ca2**2)]])
        LvP = np.array([[ a*sa2,-a  /sr3*sa2, a*np.sqrt(1-4/3*sa2**2)],
                        [     0, a*2/sr3*sa2, a*np.sqrt(1-4/3*sa2**2)],
                        [-a*sa2,-a  /sr3*sa2, a*np.sqrt(1-4/3*sa2**2)]])
        LvK = GetRecLv(LvP)
    elif BraLat[0] == "m":
        # a, b <= c, alpha < 90, beta = gamma = 90
        if   BraLat == "mP":
            LatName = "MCL"
            # LvP = np.array([[ a,    0,    0],
            #                 [ 0,    b,    0],
            #                 [ 0, c*ca, c*sa]])
            LvP = np.array([[    0,    0, a],
                            [    b,    0, 0],
                            [ c*ca, c*sa, 0]])
            LvK = GetRecLv(LvP)
        elif BraLat == "mS":
            LvP = np.array([[ a/2,  b/2,    0],
                            [-a/2,  b/2,    0],
                            [   0, c*ca, c*sa]])
            LvK = GetRecLv(LvP)
            ak, bk = LvK[0], LvK[1]
            gammak = np.arccos(ak@bk/np.linalg.norm(ak)/np.linalg.norm(bk))/np.pi*180
            if abs(gammak - 90) < error:
                LatName = "MCLC2"
            elif gammak > 90:
                LatName = "MCLC1"
            elif gammak < 90:
                x = b/c*ca + b**2/a**2*sa**2 - 1
                if abs(x) < error:
                    LatName = "MCLC4"
                elif x < 0:
                    LatName = "MCLC3"
                elif x > 0:
                    LatName = "MCLC5"
    elif BraLat == "aP":
        LvP = np.array([[    a,    0,    0],
                        [ b*cg, b*sg,    0],
                        [ c*cb, c/sg*(ca-cb*cg), c/sg*np.sqrt(sg**2-ca**2-cb**2+2*ca*cb*cg)]])
        LvK = GetRecLv(LvP)
        ak, bk, ck = LvK[0], LvK[1], LvK[2]
        alphak = np.arccos(bk@ck/np.linalg.norm(bk)/np.linalg.norm(ck))/np.pi*180
        betak  = np.arccos(ck@ak/np.linalg.norm(ck)/np.linalg.norm(ak))/np.pi*180
        gammak = np.arccos(ak@bk/np.linalg.norm(ak)/np.linalg.norm(bk))/np.pi*180
        if abs(gammak - 90) < error:
            if   alphak > 90 and betak > 90:
                LatName = "TRI2a"
            elif alphak < 90 and betak < 90:
                LatName = "TRI2b"
        elif alphak > 90 and betak > 90 and gammak > 90 and abs(gammak - np.min([alphak,betak,gammak])) < error:
            LatName = "TRI1a"
        elif alphak < 90 and betak < 90 and gammak < 90 and abs(gammak - np.max([alphak,betak,gammak])) < error:
            LatName = "TRI1b"
        elif abs(gammak - np.min([alphak,betak,gammak])) < error:
            LatName = "TRIa" # only used in this code
        elif abs(gammak - np.min([alphak,betak,gammak])) < error:
            LatName = "TRIb" # only used in this code
            
    return LvP, LvK, LatName

def get_high_sym_2d(LatName,LvP,LvC_ab):
    a, b, gamma = LvC_ab
    cg = np.cos(gamma/180*np.pi)
    if   LatName == "SQR":
        HiSymKptName = ["${\Gamma}$","M","X"]
        HiSymKptLv = np.array([[   0,   0,   0],
                               [ 1/2, 1/2,   0],
                               [   0, 1/2,   0]])
        HiSymLine = [["${\Gamma}$","X","M","${\Gamma}$"]]
    elif LatName == "HEX2D":
        HiSymKptName = ["${\Gamma}$","M","K"]
        HiSymKptLv = np.array([[   0,   0,   0],
                               [ 1/2,   0,   0],
                               [ 1/3, 1/3,   0]])
        HiSymLine = [["${\Gamma}$","M","K","${\Gamma}$"]]
    elif LatName == "RECT":
        HiSymKptName = ["${\Gamma}$","X","Y","S"]
        HiSymKptLv = np.array([[   0,   0,   0],
                               [ 1/2,   0,   0],
                               [   0, 1/2,   0],
                               [ 1/2, 1/2,   0]])
        HiSymLine = [["${\Gamma}$","X","S","Y","${\Gamma}$","S"]]
    elif LatName == "CRECT":
        HiSymKptName = ["${\Gamma}$","X","Y","$A_1$"]
        eta = (1+a**2/b**2)/4
        HiSymKptLv = np.array([[   0,   0,   0],
                               [ eta,-eta,   0],
                               [ 1/2, 1/2,   0],
                               [1-eta,eta,   0]])
        HiSymLine = [["${\Gamma}$","X","$A_1$","Y","${\Gamma}$"]]
    elif LatName == "OBL":
        HiSymKptName = ["${\Gamma}$","X","Y","C","H","$H_1$"]
        eta  = 1/2*(1-a/b*cg)/(1-cg**2)
        zeta = 1/2*(1-b/a*cg)/(1-cg**2)
        HiSymKptLv = np.array([[   0,   0,   0],
                               [ 1/2,   0,   0],
                               [   0, 1/2,   0],
                               [ 1/2, 1/2,   0],
                               [ eta,1-zeta, 0],
                               [1-eta,zeta,  0]])
        HiSymLine = [["${\Gamma}$","Y","H","C","$H_1$","X","${\Gamma}$"]]
    return HiSymKptName, HiSymKptLv, HiSymLine

def get_high_sym_3d(LatName,LvP,LvC_abc):
    a, b, c, alpha, beta, gamma = LvC_abc
    ca  = np.cos(alpha/180*np.pi)
    sa  = np.sin(alpha/180*np.pi)
    ca2 = np.cos(alpha/180*np.pi/2)
    sa2 = np.sin(alpha/180*np.pi/2)
    if   LatName == "CUB":
        HiSymKptName = ["${\Gamma}$","M","R","X"]
        HiSymKptLv = np.array([[   0,   0,   0],
                               [ 1/2, 1/2,   0],
                               [ 1/2, 1/2, 1/2],
                               [   0, 1/2,   0]])
        HiSymLine = [["${\Gamma}$","X","M","${\Gamma}$","R","X"],
                     ["M","R"]]
    elif LatName == "FCC":
        HiSymKptName = ["${\Gamma}$","K","L","U","W","X"]
        HiSymKptLv = np.array([[   0,   0,   0],
                               [ 3/8, 3/8, 3/4],
                               [ 1/2, 1/2, 1/2],
                               [ 5/8, 1/4, 5/8],
                               [ 1/2, 1/4, 3/4],
                               [ 1/2,   0, 1/2]])
        HiSymLine = [["${\Gamma}$","X","W","K","${\Gamma}$","L","U","W","L","K"],
                     ["U","X"]]
    elif LatName == "BCC":
        HiSymKptName = ["${\Gamma}$","H","P","N"]
        HiSymKptLv = np.array([[   0,   0,   0],
                               [ 1/2,-1/2, 1/2],
                               [ 1/4, 1/4, 1/4],
                               [   0,   0, 1/2]])
        HiSymLine = [["${\Gamma}$","H","N","${\Gamma}$","P","H"],
                     ["P","N"]]
    elif LatName == "TET":
        HiSymKptName = ["${\Gamma}$","A","M","R","X","Z"]
        HiSymKptLv = np.array([[   0,   0,   0],
                               [ 1/2, 1/2, 1/2],
                               [ 1/2, 1/2,   0],
                               [   0, 1/2, 1/2],
                               [   0, 1/2,   0],
                               [   0,   0, 1/2]])
        HiSymLine = [["${\Gamma}$","X","M","${\Gamma}$","Z","R","A","Z"],
                     ["X","R"],["M","A"]]
    elif LatName == "BCT1":
        eta = (1+c**2/a**2)/4
        HiSymKptName = ["${\Gamma}$","M","N","P","X","Z","$Z_1$"]
        HiSymKptLv = np.array([[   0,   0,   0],
                               [-1/2, 1/2, 1/2],
                               [   0, 1/2,   0],
                               [ 1/4, 1/4, 1/4],
                               [   0,   0, 1/2],
                               [ eta, eta,-eta],
                               [-eta,1-eta,eta]])
        HiSymLine = [["${\Gamma}$","X","M","${\Gamma}$","Z","P","N","$Z_1$","M"],
                     ["X","P"]]
    elif LatName == "BCT2":
        eta = (1+a**2/c**2)/4
        zeta = a**2/c**2/2
        HiSymKptName = ["${\Gamma}$","N","P","${\Sigma}$","${\Sigma_1}$","X","Y","$Y_1$","Z"]
        HiSymKptLv = np.array([[   0,   0,   0],
                               [   0, 1/2,   0],
                               [ 1/4, 1/4, 1/4],
                               [-eta, eta, eta],
                               [eta,1-eta,-eta],
                               [   0,   0, 1/2],
                               [-zeta,zeta,1/2],
                               [ 1/2,1/2,-zeta],
                               [ 1/2, 1/2,-1/2]])
        HiSymLine = [["${\Gamma}$","X","Y","${\Sigma}$","${\Gamma}$","Z","${\Sigma_1}$","N","P","$Y_1$","Z"],
                     ["X","P"]]
    elif LatName == "ORC":
        HiSymKptName = ["${\Gamma}$","R","S","T","U","X","Y","Z"]
        HiSymKptLv = np.array([[   0,   0,   0],
                               [ 1/2, 1/2, 1/2],
                               [ 1/2, 1/2,   0],
                               [   0, 1/2, 1/2],
                               [ 1/2,   0, 1/2],
                               [ 1/2,   0,   0],
                               [   0, 1/2,   0],
                               [   0,   0, 1/2]])      
        HiSymLine = [["${\Gamma}$","X","S","Y","${\Gamma}$","Z","U","R","T","Z"],
                     ["Y","T"],["U","X"],["S","R"]]
    elif LatName == "ORCF1":
        zeta = (1+a**2/b**2-a**2/c**2)/4
        eta = (1+a**2/b**2+a**2/c**2)/4
        HiSymKptName = ["${\Gamma}$","A","$A_1$","L","T","X","$X_1$","Y","Z"]
        HiSymKptLv = np.array([[   0,   0,   0],
                               [ 1/2, 1/2+zeta,zeta],
                               [ 1/2, 1/2-zeta,1-zeta],
                               [ 1/2, 1/2, 1/2],
                               [   1, 1/2, 1/2],
                               [   0, eta, eta],
                               [   1,1-eta,1-eta],
                               [ 1/2,   0, 1/2],
                               [ 1/2, 1/2,   0]])
        HiSymLine = [["${\Gamma}$","Y","T","Z","${\Gamma}$","X","$A_1$","Y"],
                     ["T","$X_1$"],["X","A","Z"],["L","${\Gamma}$"]]
    elif LatName == "ORCF2":
        eta   = (1+a**2/b**2-a**2/c**2)/4
        phi   = (1+c**2/b**2-c**2/a**2)/4
        delta = (1+b**2/a**2-b**2/c**2)/4
        HiSymKptName = ["${\Gamma}$","C","$C_1$","D","$D_1$","L","H","$H_1$","X","Y","Z"]
        HiSymKptLv = np.array([[   0,   0,   0],
                               [ 1/2, 1/2-eta,1-eta],
                               [ 1/2, 1/2+eta,eta],
                               [ 1/2-delta, 1/2,1-delta],
                               [ 1/2+delta, 1/2,delta],
                               [ 1/2, 1/2, 1/2],
                               [1-phi,1/2-phi,1/2],
                               [ phi,1/2+phi,1/2],
                               [   0, 1/2, 1/2],
                               [ 1/2,   0, 1/2],
                               [ 1/2, 1/2,   0]])
        HiSymLine = [["${\Gamma}$","Y","C","D","X","${\Gamma}$","Z","$D_1$","H","C"],
                     ["$C_1$","Z"],["X","$H_1$"],["H","Y"],["L","${\Gamma}$"]]
    elif LatName == "ORCF3":
        zeta = (1+a**2/b**2-a**2/c**2)/4
        HiSymKptName = ["${\Gamma}$","A","$A_1$","L","T","X","Y","Z"]
        HiSymKptLv = np.array([[   0,   0,   0],
                               [ 1/2, 1/2+zeta,zeta],
                               [ 1/2, 1/2-zeta,1-zeta],
                               [ 1/2, 1/2, 1/2],
                               [   1, 1/2, 1/2],
                               [   0, 1/2, 1/2],
                               [ 1/2,   0, 1/2],
                               [ 1/2, 1/2,   0]])
        HiSymLine = [["${\Gamma}$","Y","T","Z","${\Gamma}$","X","$A_1$","Y"],
                     ["X","A","Z"],["L","${\Gamma}$"]]
    elif LatName == "ORCI":
        zeta  = (1+a**2/c**2)/4
        eta   = (1+b**2/c**2)/4
        delta = (b**2/c**2-a**2/c**2)/4
        mu    = (b**2/c**2+a**2/c**2)/4
        HiSymKptName = ["${\Gamma}$","L","$L_1$","$L_2$","R","S","T","W","X","$X_1$","Y","$Y_1$","Z"]
        HiSymKptLv = np.array([[   0,   0,   0],
                               [ -mu,  mu,1/2-delta],
                               [  mu, -mu,1/2+delta],
                               [1/2-delta,1/2+delta,-mu],
                               [   0, 1/2,   0],
                               [ 1/2,   0,   0],
                               [   0,   0, 1/2],
                               [ 1/4, 1/4, 1/4],
                               [-zeta,zeta,zeta],
                               [zeta,1-zeta,-zeta],
                               [ eta,-eta, eta],
                               [1-eta,eta,-eta],
                               [ 1/2, 1/2,-1/2]])
        HiSymLine = [["${\Gamma}$","X","L","T","W","R","$X_1$","Z","${\Gamma}$","Y","S","W"],
                     ["$L_1$","Y"],["$Y_1$","Z"]]
    elif LatName == "ORCC":
        zeta = (1+a**2/b**2)/4
        HiSymKptName = ["${\Gamma}$","A","$A_1$","R","S","T","X","$X_1$","Y","Z"]
        HiSymKptLv = np.array([[   0,   0,   0],
                               [zeta,zeta, 1/2],
                               [-zeta,1-zeta,1/2],
                               [   0, 1/2, 1/2],
                               [   0, 1/2,   0],
                               [-1/2, 1/2, 1/2],
                               [zeta,zeta,   0],
                               [-zeta,1-zeta,0],
                               [-1/2, 1/2,   0],
                               [   0,   0, 1/2]])
        HiSymLine = [["${\Gamma}$","X","S","R","A","Z","${\Gamma}$","$X_1$","$A_1$","T","Y"],
                     ["Z","T"]]
    elif LatName == "HEX":
        HiSymKptName = ["${\Gamma}$","A","H","K","L","M"]
        HiSymKptLv = np.array([[   0,   0,   0],
                               [   0,   0, 1/2],
                               [ 1/3, 1/3, 1/2],
                               [ 1/3, 1/3,   0],
                               [ 1/2,   0, 1/2],
                               [ 1/2,   0,   0]])
        HiSymLine = [["${\Gamma}$","M","K","${\Gamma}$","A","L","H","A"],
                     ["L","M"],["K","H"]]
    elif LatName == "RHL1":
        eta = (1+4*ca)/(2+4*ca)
        nu  = 3/4-eta/2
        HiSymKptName = ["${\Gamma}$","B","$B_1$","F","L","$L_1$","P","$P_1$","$P_2$","Q","X","Z"]
        HiSymKptLv = np.array([[   0,   0,   0],
                               [ eta, 1/2,1-eta],
                               [ 1/2,1-eta,eta-1],
                               [ 1/2, 1/2,   0],
                               [ 1/2,   0,   0],
                               [   0,   0,-1/2],
                               [ eta,  nu,  nu],
                               [1-nu,1-nu,1-eta],
                               [  nu,  nu,eta-1],
                               [1-nu,  nu,   0],
                               [  nu,   0, -nu],
                               [ 1/2, 1/2, 1/2]])
        HiSymLine = [["${\Gamma}$","L","$B_1$"],["B","Z","${\Gamma}$","X"],
                     ["Q","F","$P_1$","Z"],["L","P"]]
    elif LatName == "RHL2":
        eta = 1/2*ca2**2/sa2**2
        nu  = 3/4-eta/2
        HiSymKptName = ["${\Gamma}$","F","L","P","$P_1$","Q","$Q_1$","Z"]
        HiSymKptLv = np.array([[   0,   0,   0],
                               [ 1/2,-1/2,   0],
                               [ 1/2,   0,   0],
                               [1-nu, -nu,1-nu],
                               [  nu,nu-1,nu-1],
                               [ eta, eta, eta],
                               [1-eta,-eta,-eta],
                               [ 1/2,-1/2, 1/2]])
        HiSymLine = [["${\Gamma}$","P","Z","Q","${\Gamma}$","F","$P_1$","$Q_1$","L","Z"]]
    elif LatName == "MCL":
        eta = (1-b/c*ca)/(2*sa**2)
        nu  = 1/2-eta*c/b*ca
        HiSymKptName = ["${\Gamma}$","A","C","D","$D_1$","E","H","$H_1$","$H_2$","M","$M_1$","$M_2$","X","Y","$Y_1$","Z"]
        HiSymKptLv = np.array([[   0,   0,   0],
                               [ 1/2, 1/2,   0],
                               [   0, 1/2, 1/2],
                               [ 1/2,   0, 1/2],
                               [ 1/2,   0,-1/2],
                               [ 1/2, 1/2, 1/2],
                               [   0, eta,1-nu],
                               [   0,1-eta, nu],
                               [   0, eta, -nu],
                               [ 1/2, eta,1-nu],
                               [ 1/2,1-eta, nu],
                               [ 1/2, eta, -nu],
                               [   0, 1/2,   0],
                               [   0,   0, 1/2],
                               [   0,   0,-1/2],
                               [ 1/2,   0,   0]])    
        HiSymLine = [["${\Gamma}$","Y","H","C","E","$M_1$","A","X","$H_1$"],
                     ["M","D","Z"],["Y","D"]]         
    elif LatName == "MCLC1":
        zeta = (2-b/c*ca)/(4*sa**2)
        eta  = 1/2+2*zeta*c/b*ca
        psi  = 3/4-a**2/(4*b**2*sa**2)
        phi  = psi+(3/4-psi)*b/c*ca
        HiSymKptName = ["${\Gamma}$","N","$N_1$","F","$F_1$","$F_2$","I","$I_1$","L","M","X","$X_1$","$X_2$","Y","$Y_1$","Z"]
        HiSymKptLv = np.array([[   0,   0,   0],
                               [ 1/2,   0,   0],
                               [   0,-1/2,   0],
                               [1-zeta,1-zeta,1-eta],
                               [zeta,zeta, eta],
                               [-zeta,-zeta,1-eta],
                               [ phi,1-phi,1/2],
                               [1-phi,phi-1,1/2],
                               [ 1/2, 1/2, 1/2],
                               [ 1/2,   0, 1/2],
                               [1-psi,psi-1, 0],
                               [ psi,1-psi,  0],
                               [psi-1,-psi,  0],
                               [ 1/2, 1/2,   0],
                               [-1/2,-1/2,   0],
                               [   0,   0, 1/2]])  
        HiSymLine = [["${\Gamma}$","Y","F","L","I"],["$I_1$","Z","$F_1$"],["Y","$X_1$"],
                     ["X","${\Gamma}$","N"],["M","${\Gamma}$"]]
    elif LatName == "MCLC2":
        zeta = (2-b/c*ca)/(4*sa**2)
        eta  = 1/2+2*zeta*c/b*ca
        psi  = 3/4-a**2/(4*b**2*sa**2)
        phi  = psi+(3/4-psi)*b/c*ca
        HiSymKptName = ["${\Gamma}$","N","$N_1$","F","$F_1$","$F_2$","$F_3$","I","$I_1$","L","M","X","Y","$Y_1$","Z"]
        HiSymKptLv = np.array([[   0,   0,   0],
                               [ 1/2,   0,   0],
                               [   0,-1/2,   0],
                               [1-zeta,1-zeta,1-eta],
                               [zeta,zeta, eta],
                               [-zeta,-zeta,1-eta],
                               [1-zeta,-zeta,1-eta],
                               [ phi,1-phi,1/2],
                               [1-phi,phi-1,1/2],
                               [ 1/2, 1/2, 1/2],
                               [ 1/2,   0, 1/2],
                               [1-psi,psi-1, 0],
                               [ 1/2, 1/2,   0],
                               [-1/2,-1/2,   0],
                               [   0,   0, 1/2]])  
        HiSymLine = [["${\Gamma}$","Y","F","L","I"],["$I_1$","Z","$F_1$"],
                     ["N","${\Gamma}$","M"]]            
    elif LatName == "MCLC3":
        mu    = (1+b**2/a**2)/4
        delta = b*c*ca/(2*a**2)
        zeta  = mu-1/4+(1-b/c*ca)/(4*sa**2)
        eta   = 1/2+2*zeta*c/b*ca
        phi   = 1+zeta-2*mu
        psi   = eta-2*delta
        HiSymKptName = ["${\Gamma}$","F","$F_1$","$F_2$","H","$H_1$","$H_2$","I","M","N","$N_1$","X","Y","$Y_1$","$Y_2$","$Y_3$","Z"]
        HiSymKptLv = np.array([[   0,   0,   0],
                               [1-phi,1-phi,1-psi],
                               [ phi,phi-1,psi],
                               [1-phi,-phi,1-psi],
                               [zeta,zeta, eta],
                               [1-zeta,-zeta,1-eta],
                               [-zeta,-zeta,1-eta],
                               [ 1/2,-1/2, 1/2],
                               [ 1/2,   0, 1/2],
                               [ 1/2,   0,   0],
                               [   0,-1/2,   0],
                               [ 1/2,-1/2,   0],
                               [  mu,  mu,delta],
                               [1-mu, -mu,-delta],
                               [ -mu, -mu,-delta],
                               [  mu,mu-1,delta],
                               [   0,   0, 1/2]])
        HiSymLine = [["${\Gamma}$","Y","F","H","Z","I","$F_1$"],["$H_1$","$Y_1$","X","${\Gamma}$","N"],
                     ["M","${\Gamma}$"]]
    elif LatName == "MCLC4":
        mu    = (1+b**2/a**2)/4
        delta = b*c*ca/(2*a**2)
        zeta  = mu-1/4+(1-b/c*ca)/(4*sa**2)
        eta   = 1/2+2*zeta*c/b*ca
        phi   = 1+zeta-2*mu
        psi   = eta-2*delta
        HiSymKptName = ["${\Gamma}$","F","H","$H_1$","$H_2$","I","M","N","$N_1$","X","Y","$Y_1$","$Y_2$","$Y_3$","Z"]
        HiSymKptLv = np.array([[   0,   0,   0],
                               [1-phi,1-phi,1-psi],
                               [zeta,zeta, eta],
                               [1-zeta,-zeta,1-eta],
                               [-zeta,-zeta,1-eta],
                               [ 1/2,-1/2, 1/2],
                               [ 1/2,   0, 1/2],
                               [ 1/2,   0,   0],
                               [   0,-1/2,   0],
                               [ 1/2,-1/2,   0],
                               [  mu,  mu,delta],
                               [1-mu, -mu,-delta],
                               [ -mu, -mu,-delta],
                               [  mu,mu-1,delta],
                               [   0,   0, 1/2]])
        HiSymLine = [["${\Gamma}$","Y","F","H","Z","I"],["$H_1$","$Y_1$","X","${\Gamma}$","N"],
                     ["M","${\Gamma}$"]]
    elif LatName == "MCLC5":
        zeta  = (b**2/a**2+(1-b/c*ca)/sa**2)/4
        eta   = 1/2+2*zeta*c/b*ca
        mu    = eta/2+b**2/(4*a**2)-b*c*ca/(2*a**2)
        nu    = 2*mu-zeta
        omega = (4*mu-1-b**2*sa**2/a**2)*c/(2*b*ca)
        delta = zeta*c/b*ca+omega/2-1/4
        rho   = 1-zeta*a**2/b**2
        HiSymKptName = ["${\Gamma}$","F","$F_1$","$F_2$","H","$H_1$","$H_2$","I","$I_1$","L","M","N","$N_1$","X","Y","$Y_1$","$Y_2$","$Y_3$","Z"]
        HiSymKptLv = np.array([[   0,   0,   0],
                               [  nu,  nu,omega],
                               [1-nu,1-nu,1-omega],
                               [  nu,nu-1,omega],
                               [zeta,zeta, eta],
                               [1-zeta,-zeta,1-eta],
                               [-zeta,-zeta,1-eta],
                               [ rho,1-rho,1/2],
                               [1-rho,rho-1, 1/2],
                               [ 1/2, 1/2, 1/2],
                               [ 1/2,   0, 1/2],
                               [ 1/2,   0,   0],
                               [   0,-1/2,   0],
                               [ 1/2,-1/2,   0],
                               [  mu,  mu,delta],
                               [1-mu, -mu,-delta],
                               [ -mu, -mu,-delta],
                               [  mu,mu-1,delta],
                               [   0,   0, 1/2]])
        HiSymLine = [["${\Gamma}$","Y","F","L","I"],["$I_1$","Z","H","$F_1$"],
                     ["$H_1$","$Y_1$","X","${\Gamma}$","N"],["M","${\Gamma}$"]]
    elif LatName == "TRI1a" or LatName == "TRI2a" or LatName == "TRIa":
        HiSymKptName = ["${\Gamma}$","L","M","N","R","X","Y","Z"]
        HiSymKptLv = np.array([[   0,   0,   0],
                               [ 1/2, 1/2,   0],
                               [   0, 1/2, 1/2],
                               [ 1/2,   0, 1/2],
                               [ 1/2, 1/2, 1/2],
                               [ 1/2,   0,   0],
                               [   0, 1/2,   0],
                               [   0,   0, 1/2]])
        HiSymLine = [["X","${\Gamma}$","Y"],["L","${\Gamma}$","Z"],["N","${\Gamma}$","M"],
                     ["R","${\Gamma}$"]]
    elif LatName == "TRI1b" or LatName == "TRI2b" or LatName == "TRIb":
        HiSymKptName = ["${\Gamma}$","L","M","N","R","X","Y","Z"]
        HiSymKptLv = np.array([[   0,   0,   0],
                               [ 1/2,-1/2,   0],
                               [   0,   0, 1/2],
                               [-1/2,-1/2, 1/2],
                               [   0,-1/2, 1/2],
                               [   0,-1/2,   0],
                               [ 1/2,   0,   0],
                               [-1/2,   0, 1/2]])
        HiSymLine = [["X","${\Gamma}$","Y"],["L","${\Gamma}$","Z"],["N","${\Gamma}$","M"],
                     ["R","${\Gamma}$"]]
    return HiSymKptName, HiSymKptLv, HiSymLine

def get_brillouin_zone_2d(cell):
    
    cell = np.asarray(cell, dtype=float)
    assert cell.shape == (2, 2)

    px, py = np.tensordot(cell, np.mgrid[-1:2, -1:2], axes=[0,0])
    points = np.c_[px.ravel(), py.ravel()]
    
    from scipy.spatial import Voronoi
    vor = Voronoi(points)
    
    bz_ridges = []
    bz_vertices = []


    for pid, rid in zip(vor.ridge_points, vor.ridge_vertices):
        if(pid[0] == 4 or pid[1] == 4):
            bz_ridges.append(vor.vertices[np.r_[rid, [rid[0]]]])
            bz_vertices += rid

    bz_vertices = list(set(bz_vertices))

    return vor.vertices[bz_vertices], bz_ridges
    
def get_brillouin_zone_3d(cell):
    """
    Ref: http://staff.ustc.edu.cn/~zqj/posts/howto-plot-brillouin-zone/
    
    Generate the Brillouin Zone of a given cell. The BZ is the Wigner-Seitz cell
    of the reciprocal lattice, which can be constructed by Voronoi decomposition
    to the reciprocal lattice.  A Voronoi diagram is a subdivision of the space
    into the nearest neighborhoods of a given set of points. 

    https://en.wikipedia.org/wiki/Wigner%E2%80%93Seitz_cell
    https://docs.scipy.org/doc/scipy/reference/tutorial/spatial.html#voronoi-diagrams
    """

    cell = np.asarray(cell, dtype=float)
    assert cell.shape == (3, 3)

    px, py, pz = np.tensordot(cell, np.mgrid[-1:2, -1:2, -1:2], axes=[0, 0])
    points = np.c_[px.ravel(), py.ravel(), pz.ravel()]

    from scipy.spatial import Voronoi
    vor = Voronoi(points)

    bz_facets = []
    bz_ridges = []
    bz_vertices = []

    # for rid in vor.ridge_vertices:
    #     if( np.all(np.array(rid) >= 0) ):
    #         bz_ridges.append(vor.vertices[np.r_[rid, [rid[0]]]])
    #         bz_facets.append(vor.vertices[rid])

    for pid, rid in zip(vor.ridge_points, vor.ridge_vertices):
        # WHY 13 ????
        # The Voronoi ridges/facets are perpendicular to the lines drawn between the
        # input points. The 14th input point is [0, 0, 0].
        if(pid[0] == 13 or pid[1] == 13):
            bz_ridges.append(vor.vertices[np.r_[rid, [rid[0]]]])
            bz_facets.append(vor.vertices[rid])
            bz_vertices += rid

    bz_vertices = list(set(bz_vertices))

    return vor.vertices[bz_vertices], bz_ridges, bz_facets

if __name__ == '__main__':
    plt.close("all")
    # BraLat = ["cP","cF","cI","tP","tI","oP","oF","oI","oS","hP","hR","mP","mS","aP"]
    # for i in range(len(BraLat)):
    #     print(BraLat[i])
    #     PlotBrillouinZone(3,BraLat[i])
    BraLat = ["tP","hP","oP","oC","mP"]
    for i in range(len(BraLat)):
        print(BraLat[i])
        PlotBrillouinZone(2,BraLat[i])
    # PlotBrillouinZone(2,"tP")