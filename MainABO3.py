import sys
import numpy as np
import matplotlib.pyplot as plt
plt.close("all")

# Main
from cd.SymGroup import GetSpaceGroup
from cd.NbrAtom  import FindNeighbor
from cd.SymAtom  import FindAtomSymmetry
from cd.HopRel   import FindRelation
from cd.HopVal   import GetHoppingValue
from cd.HmtReal  import GetHamiltonianReal
from cd.KPoint   import GetKPoint
from cd.HmtK     import GetHamiltonianK
from cd.Eig      import GetEigenSolution

# Read & Write
from rw.ReadTBIN import ReadInput
from rw.ReadHop  import ReadHopping
from rw.ReadRel  import ReadRelation
from rw.ReadHmt  import ReadHamiltonianReal
from rw.ReadKpt  import ReadKPoint
from rw.WriteRel import WriteRelation
from rw.WriteHmt import WriteHamiltonianReal

# Plot
from pl.PlotEB   import PlotEnergyBand
from pl.PlotAt   import PlotAtoms

# Check
# from ck.CheckWRR import CheckWriteReadRelation
# from ck.CheckSym import CheckSymmetry

if __name__ == '__main__':
    
    Name = "ABO3"
    
    # From input to hopping relations
    ParaIn     = ReadInput("data/TBIN_" + Name + ".txt")
    ParaSym    = GetSpaceGroup(ParaIn)
    # CheckSymmetry(ParaSym)
    ParaNbr    = FindNeighbor(ParaIn)
    PlotAtoms(ParaNbr,Name)
    ParaSymAt  = FindAtomSymmetry(ParaIn,ParaSym,ParaNbr)
    ParaRel    = FindRelation(ParaIn,ParaSym,ParaNbr,ParaSymAt)
    WriteRelation(ParaIn,ParaRel)
    sys.exit(0)
    # From hopping relations and Hamiltonian in real space
    HopValIn   = ReadHopping( "data//HopValIN_" + Name + ".txt")
    ParaRel    = ReadRelation("data//" + Name + "//HopRel.txt")
    # CheckWriteReadRelation(ParaRel,ParaRel1)
    HopValClas = GetHoppingValue(HopValIn, ParaRel)
    ParaHmtR   = GetHamiltonianReal(ParaRel,HopValClas)
    WriteHamiltonianReal(ParaIn, ParaHmtR)
    
    # From Hamiltonian in real space to Hamiltonian in k space
    ParaHmtR   = ReadHamiltonianReal(Name)
    ParaKptIn  = ReadKPoint("data//KptIN_" + Name + ".txt")
    ParaKpt    = GetKPoint(ParaKptIn)
    ParaHmtK   = GetHamiltonianK(ParaHmtR,ParaKpt)
    ParaEig    = GetEigenSolution(ParaHmtK)
    ParaEigPlt = PlotEnergyBand(ParaKpt,ParaEig,Name)
    
    # print(ParaIn.keys())