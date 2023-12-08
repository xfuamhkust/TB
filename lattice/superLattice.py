from lattice.baseLattice import baseLattice
import numpy as np
import copy
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
        #TODO: for now we make the restriction that the substituted atom is different from the new atom
        #duplicate the baseLattice's atoms Nj times along direction j, j=0,1,2

        #construct set cellRowAtm
        cellRowAtm=set()

        for n0 in range(0,N0):
            for n1 in range(0,N1):
                for n2 in range(0,N2):
                    for j,ind in enumerate(self.ParaIn["AtomTypeIndex"]):
                        atm=self.ParaIn["AtomName"][ind]
                        cellRowAtm.add((n0,n1,n2,j,atm))
        # delete elements from cellRowAtm according to supercellVacancy
        for row in self.ParaIn["supercell"]["supercellVacancy"]:
            n0n1n2, k, atm = row
            n0, n1, n2 = n0n1n2
            toDelete=(n0,n1,n2,self.map2RowNum(k,atm),atm)
            cellRowAtm.remove(toDelete)
        # construct superLattice from supercellSubstitution
        for row in self.ParaIn["supercell"]["supercellSubstitution"]:
            n0n1n2, k, atmOld, atmNewAndOrbs = row
            n0, n1, n2 = n0n1n2
            atmNew = atmNewAndOrbs[0]
            j=self.map2RowNum(k,atmOld)
            toDelete=(n0,n1,n2,j,atmOld)
            cellRowAtm.remove(toDelete)
            toAdd=(n0,n1,n2,j,atmNew)
            cellRowAtm.add(toAdd)

        # construct supercellLattice from supercellInterstitial
        self.superLatticeParaIn["AtomSite"]=copy.deepcopy(self.ParaIn["AtomSite"])
        for row in self.ParaIn["supercell"]["supercellInterstitial"]:
            n0n1n2, s0s1s2, atm, orbs = row
            n0, n1, n2 = n0n1n2
            s0, s1, s2 = s0s1s2
            self.superLatticeParaIn["AtomSite"]=np.vstack([self.superLatticeParaIn["AtomSite"],[s0,s1,s2]])
            l=len(self.superLatticeParaIn["AtomSite"])
            toAdd=(n0,n1,n2,l-1,atm)
            cellRowAtm.add(toAdd)







        return

    def map2RowNum(self,k,atm):
        """

        :param k: ordinal of atm among the same atms in a baseLattice
        :param atm: atom name
        :return: row of atm in self.ParaIn["AtomSite"], in order to construct superLattice's parameters from
        supercellVacancy and supercellSubstitution
        """

        #index of atm  in all atom names
        ind=self.ParaIn["AtomName"].index(atm)

        #posistions of ind in self.ParaIn["AtomTypeIndex"], which is also the row number of atm in self.ParaIn["AtomSite"]
        positions=[i for i,val in enumerate(self.ParaIn["AtomTypeIndex"]) if i==ind]
        return positions[k]

