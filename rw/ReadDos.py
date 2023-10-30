def ReadDosPara(FileName):
    f = [line.strip().split() for line in open(FileName,"r").readlines()]
    Method = f[0][0]
    if Method == "Diag" or Method == "Inv":
        Emin, Emax, nE, eta = float(f[1][0]), float(f[1][1]), int(f[1][2]), float(f[1][3]),
        nk0    = int(f[2][0])
        ParaDosIn = {"Method": Method,
                     "nE":     nE,
                     "Emin":   Emin,
                     "Emax":   Emax,
                     "eta":    eta,
                     "nk0":    nk0,
                     }
    elif Method == "KPM":
        Emin, Emax, nE = float(f[1][0]), float(f[1][1]), int(f[1][2]),
        nk0    = int(f[2][0])
        nsr    = int(f[3][0]) # number of series
        ParaDosIn = {"Method": Method,
                     "nE":     nE,
                     "Emin":   Emin,
                     "Emax":   Emax,
                     "nk0":    nk0,
                     "nsr":    nsr,
                     }
    return ParaDosIn