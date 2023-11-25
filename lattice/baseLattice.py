import numpy as np
from lattice.assembleMixinModules import mixinModules




class lattice(*mixinModules):
    #base class for a lattice
    def __init__(self):
        self.keywords=["Name", "Dim", "Spin", "Nbr",
          "LatType", "LatVec", "AtTpNum","Bases",
          "AtomSite", "supercellSize","supercellVacancy","supercellSubstitution",
          "supercellInterstitial"]

        self.OrbGrp = np.array(["1S","2S","2P","2P","2P",
                   "3S","3P","3P","3P","3D","3D","3D","3D","3D",
                   "4S","4P","4P","4P","4D","4D","4D","4D","4D",
                   "4F","4F","4F","4F","4F","4F","4F",
                   "5S","5P","5P","5P","5D","5D","5D","5D","5D",
                   "5F","5F","5F","5F","5F","5F","5F",
                   "6S","6P","6P","6P","6D","6D","6D","6D","6D",
                   "6F","6F","6F","6F","6F","6F","6F",
                   "7S","7P","7P","7P","7D","7D","7D","7D","7D",
                   "7F","7F","7F","7F","7F","7F","7F",
                   "S","P","P","P","D","D","D","D","D",
                   "F","F","F","F","F","F","F"])
        self.OrbIdv = np.array(["1s","2s","2px","2py","2pz",
                   "3s","3px","3py","3pz","3dxy","3dyz","3dzx","3dx2-y2","3dz2",
                   "4s","4px","4py","4pz","4dxy","4dyz","4dzx","4dx2-y2","4dz2",
                   "4fz3","4fxz2","4fyz2","4fxyz","4fz(x2-y2)","4fx(x2-3y2)","4fy(3x2-y2)",
                   "5s","5px","5py","5pz","5dxy","5dyz","5dzx","5dx2-y2","5dz2",
                   "5fz3","5fxz2","5fyz2","5fxyz","5fz(x2-y2)","5fx(x2-3y2)","5fy(3x2-y2)",
                   "6s","6px","6py","6pz","6dxy","6dyz","6dzx","6dx2-y2","6dz2",
                   "6fz3","6fxz2","6fyz2","6fxyz","6fz(x2-y2)","6fx(x2-3y2)","6fy(3x2-y2)",
                   "7s","7px","7py","7pz","7dxy","7dyz","7dzx","7dx2-y2","7dz2",
                   "7fz3","7fxz2","7fyz2","7fxyz","7fz(x2-y2)","7fx(x2-3y2)","7fy(3x2-y2)",
                   "s","px","py","pz","dxy","dyz","dzx","dx2-y2","dz2",
                   "fz3","fxz2","fyz2","fxyz","fz(x2-y2)","fx(x2-3y2)","fy(3x2-y2)"])



