from lattice.baseLattice import baseLattice
import numpy as np

class superLattice(baseLattice):
    # class for a superlattice
    def __init__(self):
        super().__init__()

    def constructSuperLattice(self):
        """
        construct parameters for superLattice from baseLattice
        :return:
        """
        self.superLatticeParaIn={}
        self.superLatticeParaIn["Name"]=self.ParaIn["Name"]+"_supercell"

        Dim = self.ParaIn["Dimension"]
        self.superLatticeParaIn["Dimension"]=Dim

        self.superLatticeParaIn["Spin"]=self.ParaIn["Spin"]

        n=self.ParaIn["supercell"]["supercellNbr"]
        self.superLatticeParaIn["NeighborNumber"]=[n, n, 0 if Dim == 2 else n]

        self.superLatticeParaIn["Lattice type"]="super"

        #construct atom site, atom number, atom type, bases, lattice vectors

        #for now we only consider the case where base lattice is primitive
        #TODO: generalize the transformation matrix from baseLattice to superLattice, now this matrix is diagonal diag(N0,N1,N2)

        #construct lattice vector from baseLattice and supercellSize
        N0,N1,N2=self.ParaIn["supercell"]["supercellSize"]
        baseLv=self.ParaIn["LatticeVector"]
        vec0=N0*baseLv[0,:]
        vec1=N1*baseLv[1,:]
        vec2=N2*baseLv[2,:]
        self.superLatticeParaIn["LatticeVector"]=np.array([vec0,vec1,vec2])

        #

        #TODO: for now we assume that when the interstitial atoms are the same as in the baseLattice, they also have the same orbitals
        #TODO: for now we make the restriction that the substuted atom is different from the new atom
        #duplicate the baseLattice's atoms Nj times along direction j, j=0,1,2

        #construct set cellRowAtm
        cellRowAtm=set()
        # for n0 in range(0,N0):
        #     for n1 in range(0,N1):
        #         for n2 in range(0,N2):



        return