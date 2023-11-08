import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['figure.dpi']= 300

def PlotAtomStrip(ParaStrip,Name=""):
    
    '''
    This function is to plot 2D Strip.
    '''
    
    # Parameters
    Dir       = ParaStrip["Dir"]
    Wid       = ParaStrip["Wid"]
    LvS       = ParaStrip["Lv"]
    AtName    = ParaStrip["AtName"]
    AtNumS    = ParaStrip["AtNum"]
    AtLvS     = ParaStrip["AtLv"]
    AtXyzS    = AtLvS @ LvS
    AtTypeNum = len(AtNumS)
    
    # Direction
    if Dir == 1:
        # x = a1, y = a2
        Lvx =  LvS[0]
        Lvy =  LvS[1]
    elif Dir == 2:
        # x = -a2, y = a1
        Lvx = -LvS[1]
        Lvy =  LvS[0]
    
    # Plot
    fig = plt.figure()
    Len = np.max([Wid*3,5])
    # Plot cell edges
    # n1a1+n2a2 -> (n1+1)a1+n2a2 -> (n1+1)a1+(n2+1)a2 -> n1a1+(n2+1)a2 -> n1a1+n2a2
    EdgeX0 = np.array([0,LvS[0,0],LvS[0,0]+LvS[1,0],LvS[1,0],0])
    EdgeY0 = np.array([0,LvS[0,1],LvS[0,1]+LvS[1,1],LvS[1,1],0])
    for iW in range(Wid):
        for iL in range(Len):
            EdgeX = EdgeX0 + iW*Lvx[0] + iL*Lvy[0]
            EdgeY = EdgeY0 + iW*Lvx[1] + iL*Lvy[1]
            plt.plot(EdgeX,EdgeY,"-",c="gray")
    # Plot atoms
    count0 = 0
    for iAt in range(AtTypeNum):
        AtPloti = np.zeros((AtNumS[iAt]*Wid*Len,3))
        count1 = 0
        for jAt in range(AtNumS[iAt]):
            for iW in range(Wid):
                for iL in range(Len):
                    AtPloti[count1] = AtXyzS[count0] + iW*Lvx + iL*Lvy
                    count1 += 1
            count0 += 1
        plt.plot(AtPloti[:,0],AtPloti[:,1],"o",label=AtName[iAt])
        
    plt.axis("equal")
    plt.axis("off")
    plt.legend()
    plt.show()
    
    # Output
    ParaAtPlt = {"FigureAt": fig,
                  }
    
    return ParaAtPlt