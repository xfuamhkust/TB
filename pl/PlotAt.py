import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['figure.dpi']= 300

def PlotAtoms(ParaIn,ParaNbr,Name=""):
    
    '''
    This function is to plot hopping terms.
    Each atom are displayed as a point in 2D or 3D figure.
    The same atoms are painted the same color.
    Edges of primitive cells are showed with black (000) or gray(n1,n2,n3!=000), as well as lattice vectors (RG(B) arrows) at origin.
    '''
    
    # Parameters
    AtName       = ParaIn["AtomName"]
    AtTypeInd    = ParaIn["AtomTypeIndex"]
    Lv           = ParaIn["LatticeVector"]
    LvAt         = ParaNbr["LvAt"]
    LvAtXyz      = ParaNbr["LvAtXyz"]
    
    # Classify different atoms
    LvAtType = LvAt[:,-1]
    NumType = np.max(LvAtType)
    IndType = np.array([np.where(LvAtType==i)[0] for i in range(1,NumType+1)])
    NumAtType = np.max(AtTypeInd) + 1
    IndAtType = [np.where(AtTypeInd==i)[0] for i in range(NumAtType)]
    
    # Plot atoms
    fig = plt.figure()
    x = LvAtXyz[:,0]; y = LvAtXyz[:,1]; z = LvAtXyz[:,2]
    if np.linalg.norm(z):
        # 3D
        Lv3D = np.unique(LvAt[:,:3],axis=0)
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
        ax = fig.add_subplot(111, projection='3d')
        ax.set_box_aspect([1, 1, 1])
        for n1,n2,n3 in Lv3D:
            EdgeX = EdgeX0 + n1*Lv[0,0] + n2*Lv[1,0] + n3*Lv[2,0]
            EdgeY = EdgeY0 + n1*Lv[0,1] + n2*Lv[1,1] + n3*Lv[2,1]
            EdgeZ = EdgeZ0 + n1*Lv[0,2] + n2*Lv[1,2] + n3*Lv[2,2]
            plt.plot(EdgeX,EdgeY,EdgeZ,"-",c="gray")
        plt.plot(EdgeX0,EdgeY0,EdgeZ0,"-",c="k")
        plt.quiver(0., 0., 0., Lv[0,0], Lv[0,1], Lv[0,2], color='red'  )
        plt.quiver(0., 0., 0., Lv[1,0], Lv[1,1], Lv[1,2], color='green')
        plt.quiver(0., 0., 0., Lv[2,0], Lv[2,1], Lv[2,2], color='blue' )
        plt.plot([0],[0],"k.")
        # Plot atoms
        for i in range(NumAtType):
            Ind = []
            for j in IndAtType[i]:
                Ind += IndType[j].tolist()
            ax.scatter(x[Ind],y[Ind],z[Ind],label=AtName[i])
        xyzmax = (np.max(abs(LvAtXyz))+0.5*np.max(abs(Lv)))
        ax.set_xlim([-xyzmax,xyzmax])
        ax.set_ylim([-xyzmax,xyzmax])
        ax.set_zlim([-xyzmax,xyzmax])
        plt.axis("off")
        plt.legend()
    else:
        # 2D
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
        # plot lattice vectors
        plt.annotate('',ha = 'center', va = 'bottom',xytext = (0., 0.),xy = (Lv[0,0], Lv[0,1]),arrowprops = { 'edgecolor' : 'red',   'facecolor' : 'red',   'shrink' : 0.05, 'alpha' : 0.8})
        plt.annotate('',ha = 'center', va = 'bottom',xytext = (0., 0.),xy = (Lv[1,0], Lv[1,1]),arrowprops = { 'edgecolor' : 'green', 'facecolor' : 'green', 'shrink' : 0.05, 'alpha' : 0.8})
        plt.plot([0],[0],"k.")
        # Plot atoms
        for i in range(NumAtType):
            Ind = []
            for j in IndAtType[i]:
                Ind += IndType[j].tolist()
            plt.plot(x[Ind],y[Ind],"o",label=AtName[i])
        plt.axis("equal")
        plt.axis("off")
        plt.legend()
    plt.show()
    
    # Output
    ParaAtPlt = {"FigureAt": fig,
                  }
    
    return ParaAtPlt