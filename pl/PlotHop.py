import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
plt.rcParams['figure.dpi']= 300
ColorBank = plt.rcParams["axes.prop_cycle"]

def PlotHoppingTerm(ParaIn,ParaNbr,ParaRel,Name="",HopInd=[]):
    
    '''
    This function is to plot hopping terms as arrows.
    The symmetry-related hopping terms are painted the same color
    The integers in input variable HopInd will be considered as the wanted orders of hopping.
    If not input HopInd is found, all the calculated hopping terms will show in 2D material, 
        while at most 7th order will show in 3D material including onsite energy.
    '''
    
    # Parameters
    AtName       = ParaIn["AtomName"]
    AtTypeInd    = ParaIn["AtomTypeIndex"]
    Lv           = ParaIn["LatticeVector"]
    LvAt         = ParaNbr["LvAt"]
    LvAtXyz      = ParaNbr["LvAtXyz"]
    LvAtAtOO     = ParaRel["LvAtAtOO"]
    
    # Classify different atoms
    LvAtType = LvAt[:,-1]
    NumType = np.max(LvAtType)
    IndType = np.array([np.where(LvAtType==i)[0] for i in range(1,NumType+1)])
    NumAtType = np.max(AtTypeInd) + 1
    IndAtType = [np.where(AtTypeInd==i)[0] for i in range(NumAtType)]
    
    # Reduce hopping terms from AtOrb-AtOrb to At-At
    LvAtAtInd = []
    for LvAtAtOOi in LvAtAtOO:
        LvAtAti = np.unique(LvAtAtOOi[:,:5],axis=0)
        LvAtAtIndi = []
        for j in range(len(LvAtAti)):
            n1, n2, n3, iAt, jAt = LvAtAti[j]
            LvAt1 = np.array([  0,  0,  0, iAt])
            LvAt2 = np.array([ n1, n2, n3, jAt])
            Ind1  = FindIndex(LvAt1,LvAt)
            Ind2  = FindIndex(LvAt2,LvAt)
            LvAtAtIndi.append([Ind1,Ind2])
        LvAtAtInd.append(LvAtAtIndi)
    
    # Plot atoms
    fig = plt.figure()
    x = LvAtXyz[:,0]; y = LvAtXyz[:,1]; z = LvAtXyz[:,2]
    if np.linalg.norm(z):
        # 3D
        if len(HopInd) == 0:
            HopInd = np.arange(1,7+1)
        ax = fig.add_subplot(111, projection='3d')
        ax.set_box_aspect([1, 1, 1])
        # plot cell edges
        # http://www.mathchina.com/bbs/forum.php?mod=viewthread&tid=2053333
        A = np.array([0,0,0])
        B = Lv[0]
        C = Lv[0] + Lv[1]
        D = Lv[1]
        E = Lv[2]
        F = Lv[0] + Lv[2]
        G = Lv[0] + Lv[1] + Lv[2]
        H = Lv[1] + Lv[2]
        EdgeX0 = np.array([A[0],B[0],C[0],B[0],F[0],G[0],F[0],E[0],H[0],E[0],A[0],D[0],H[0],G[0],C[0],D[0]])
        EdgeY0 = np.array([A[1],B[1],C[1],B[1],F[1],G[1],F[1],E[1],H[1],E[1],A[1],D[1],H[1],G[1],C[1],D[1]])
        EdgeZ0 = np.array([A[2],B[2],C[2],B[2],F[2],G[2],F[2],E[2],H[2],E[2],A[2],D[2],H[2],G[2],C[2],D[2]])
        plt.plot(EdgeX0,EdgeY0,EdgeZ0,"-",c="k")
        handles = []
        # Plot hopping terms
        LvAtAtIndPlt = []
        for i in range(len(LvAtAtInd)):
            if i+1 not in HopInd:
                continue
            LvAtAtIndi = LvAtAtInd[i]
            ic = np.mod(9-i,10)
            if LvAtAtIndi[0][0] == LvAtAtIndi[0][1]: # Onsite hopping
                for j in range(len(LvAtAtIndi)):
                    Ind0 = LvAtAtIndi[j][0]
                    x0, y0, z0 = LvAtXyz[Ind0]
                    ax.scatter(x0,y0,z0, color="C%d"%ic, s = 100)
                    LvAtAtIndPlt.append(Ind0)
                handles.append(Line2D([0], [0], color="C%d"%ic, ls='', marker='.', label='$t_{%d}$'%(i+1)))
            else:                                    # Offsite hopping
                for j in range(len(LvAtAtIndi)):
                    Ind1, Ind2 = LvAtAtIndi[j]
                    x1, y1, z1 = LvAtXyz[Ind1]
                    x2, y2, z2 = LvAtXyz[Ind2]
                    LvAtAtIndPlt.append(Ind1)
                    LvAtAtIndPlt.append(Ind2)
                    plt.quiver(x1, y1, z1, x2-x1, y2-y1, z2-z1, color="C%d"%ic, alpha=0.8)
                handles.append(Line2D([0], [0], color="C%d"%ic, ls='-', label='$t_{%d}$'%(i+1)))
        LvAtAtIndPlt = np.unique(LvAtAtIndPlt)
        # add atoms in 000
        Ind0 = np.where(np.linalg.norm(LvAt[:,:3],axis=1)==0)[0]
        LvAtAtIndPlt = list(set(Ind0) | set(LvAtAtIndPlt))
        # Plot atoms
        for i in range(NumAtType):
            Ind = []
            for j in IndAtType[i]:
                Ind += IndType[j].tolist()
            Ind = list(set(Ind) & set(LvAtAtIndPlt))
            ax.scatter(x[Ind],y[Ind],z[Ind], color="C%d"%i)
            handles.append(Line2D([0], [0], color="C%d"%i, ls='', marker='.', label=AtName[i]))
        xyzmax = (np.max(abs(LvAtXyz))+0.5*np.max(abs(Lv)))
        ax.set_xlim([-xyzmax,xyzmax])
        ax.set_ylim([-xyzmax,xyzmax])
        ax.set_zlim([-xyzmax,xyzmax])
        plt.axis("off")
        handles = handles[len(handles)-NumAtType:] + handles[:len(handles)-NumAtType]
        plt.legend(handles=handles)
    else:
        # 2D
        if len(HopInd) == 0:
            HopInd = np.arange(1,len(LvAtAtOO)+1)
        Lv2D = np.unique(LvAt[:,:2],axis=0)
        # plot cell edges
        # n1a1+n2a2 -> (n1+1)a1+n2a2 -> (n1+1)a1+(n2+1)a2 -> n1a1+(n2+1)a2 -> n1a1+n2a2
        EdgeX0 = np.array([0,Lv[0,0],Lv[0,0]+Lv[1,0],Lv[1,0],0])
        EdgeY0 = np.array([0,Lv[0,1],Lv[0,1]+Lv[1,1],Lv[1,1],0])
        for n1,n2 in Lv2D:
            EdgeX = EdgeX0 + n1*Lv[0,0] + n2*Lv[1,0]
            EdgeY = EdgeY0 + n1*Lv[0,1] + n2*Lv[1,1]
            plt.plot(EdgeX,EdgeY,"-",c="gray")
        plt.plot(EdgeX0,EdgeY0,"-",c="k") # 0 -> a1 -> a1+a2 -> a2 -> 0
        # Plot atoms
        handles = []
        for i in range(NumAtType):
            Ind = []
            for j in IndAtType[i]:
                Ind += IndType[j].tolist()
            linei, = plt.plot(x[Ind],y[Ind],"o",label=AtName[i])
            handles.append(linei)
        # Plot hopping terms
        for i in range(len(LvAtAtInd)):
            if i+1 not in HopInd:
                continue
            LvAtAtIndi = LvAtAtInd[i]
            ic = np.mod(9-i,10)
            if LvAtAtIndi[0][0] == LvAtAtIndi[0][1]: # Onsite hopping
                for j in range(len(LvAtAtIndi)):
                    Ind0 = LvAtAtIndi[j][0]
                    x0, y0, _ = LvAtXyz[Ind0]
                    plt.plot(x0,y0,".",color="C%d"%ic,label=AtName[i])
                handles.append(Line2D([0], [0], color="C%d"%ic, ls='', marker='.', label='$t_{%d}$'%(i+1)))
            else:                                    # Offsite hopping
                for j in range(len(LvAtAtIndi)):
                    Ind1, Ind2 = LvAtAtIndi[j]
                    x1, y1, _ = LvAtXyz[Ind1]
                    x2, y2, _ = LvAtXyz[Ind2]
                    plt.annotate('',xytext = (x1, y1),xy = (x2, y2),arrowprops = {'arrowstyle': '->', 'connectionstyle': 'arc3,rad=0.1', 'edgecolor' : "C%d"%ic, 'alpha' : 0.8, 'label' : '$t_%d$'%(i+1)})
                handles.append(Line2D([0], [0], color="C%d"%ic, ls='-', label='$t_{%d}$'%(i+1)))
        plt.axis("equal")
        plt.axis("off")
        plt.legend(handles=handles)
    plt.show()
    
    # Output
    ParaHopPlt = {"FigureHop": fig,
                  }
    
    return ParaHopPlt

def FindIndex(x,X,tol = 1e-3):
    dX = np.sum(abs(X - x), axis=1)
    min_dX = min(dX)
    if min_dX < tol:
        ind = np.where(dX == min_dX)[0][0]
        return ind
    else:
        return -1