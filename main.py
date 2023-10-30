import sys
import numpy as np
import matplotlib.pyplot as plt
import os
import pathlib
from datetime import datetime
import sympy as smp

plt.close("all")
# Main
# from cd.SymGroup import GetSpaceGroupPrimitive
from cd.SymGroup import paraInPrim2Conv,conv2SGN,getParaSym
from cd.NbrAtom  import FindNeighbor
from cd.SymAtom  import FindAtomSymmetry
from cd.HopRel   import FindRelation
from cd.HopVal   import GetHoppingValue
from cd.HmtReal  import GetHamiltonianReal
from cd.KPoint   import GetKPoint
from cd.HmtK     import GetHamiltonianK
from cd.HmtK     import HkSp2Np
from cd.Eig      import GetEigenSolution
from cd.Berry    import GetBerryK
from cd.Berry    import GetAHNC
from cd.Dos      import GetDos
# Read & Write
from rw.ReadTBIN import ReadInput
from rw.ReadHop  import ReadHopping
from rw.ReadRel  import ReadRelation
from rw.ReadHmt  import ReadHamiltonianReal
from rw.ReadKpt  import ReadKPoint
from rw.WriteRel import WriteRelation
from rw.WriteHmt import WriteHamiltonianReal
from rw.WriteAna import WriteHamiltonianAnalytic
from rw.WriteAna import writeHkAnalytic
from rw.printHk  import Hk2html
from rw.ReadBer  import ReadBerryPara
from rw.ReadDos  import ReadDosPara
# Plot
from pl.PlotEB   import PlotEnergyBand
from pl.PlotEB   import plotEigsFromHkSp
from pl.PlotBC   import PlotBerryCurvature
from pl.PlotDos  import PlotDOS
# from pl.PlotAt   import PlotAtoms
# from pl.PlotHop  import PlotHoppingTerm
# from pl.PlotBZ  import PlotBrillouinZone
# Check
from ck.CheckOrbCpl    import CheckOrbitalCompleteness as CheckOrbital
from ck.CheckHmtSym    import CheckHamiltonianSymmetry as CheckH
from ck.CheckEigValSym import CheckEnergySymmetry as CheckE

#This is the demo program for Auh-17-2023
#1. It reads config file containing information of the crystal,
# the reader may give either a primitive cell or a conventional cell,
# the program will compute primitive from conventional, or conventional from primitive.
#2. When the

# material = 'data/ABO3/primitive_TBIN_ABO3.txt'
material = 'data/Graphene/primitive_TBIN_Graphene.txt'
# material = 'data/h-BN/primitive_TBIN_h-BN.txt'
# material = 'data/NaCl/primitive_TBIN_NaCl.txt'
# material = 'data/Si/primitive_TBIN_Si.txt'
lenParams=len(sys.argv)
if lenParams<2:
    sys.argv.append(material)
else:
    sys.argv[1] = material

lenParams=len(sys.argv)
if lenParams!=2:
    raise RuntimeError("Wrong number of arguments.")

inConfigName=sys.argv[1]
#check if given path exist
if not os.path.exists(inConfigName):
    raise NotADirectoryError("path does not exist.")
#check if given file exists
if not os.path.isfile(inConfigName):
    raise FileNotFoundError("Not a file.")
# print(inConfigName)
# inConfigFolder=pathlib.Path(inConfigName).parent
# print(inConfigFolder)
#read info
ParaIn=ReadInput(inConfigName)
Name=ParaIn["Name"]


'''################### Determine space group of crystal ####################'''
# ParaSym = GetSpaceGroupPrimitive(ParaIn)
atmUnderConvVector,atmIndsConv=paraInPrim2Conv(ParaIn)
SGN, originBilbao=conv2SGN(ParaIn,atmUnderConvVector,atmIndsConv)
ParaSym=getParaSym(ParaIn,SGN, originBilbao)
ParaIn["origin Bilbao"]=ParaSym["origin Bilbao"]


'''##################### Complete Electronic orbitals ######################'''
ParaIn,ParaSym = CheckOrbital(ParaIn,ParaSym)


'''################### Relate hopping terms by symmetry ####################'''
ParaNbr    = FindNeighbor(ParaIn)
# PlotAtoms(ParaIn,ParaNbr,Name)
tFindingRelationStart=datetime.now()
ParaSymAt  = FindAtomSymmetry(ParaIn,ParaSym,ParaNbr)
ParaRel    = FindRelation(ParaIn,ParaSym,ParaNbr,ParaSymAt)
# PlotHoppingTerm(ParaIn,ParaNbr,ParaRel,Name,[5,6])
tFindingRelationEnd=datetime.now()
print("Finding symmetry relations: ",tFindingRelationEnd-tFindingRelationStart)
WriteRelation(ParaIn,ParaRel)


'''################### Calculate real space Hamiltonian ####################'''
# From hopping relations to Hamiltonian in real space
HopValIn   = ReadHopping(ParaIn["Folder"]+"/HopValIN_"+ParaIn["Name"]+".txt")
ParaRel    = ReadRelation(ParaIn["Folder"]+ "/HopRel.txt")
HopValClas = GetHoppingValue(HopValIn, ParaRel)
ParaHmtR   = GetHamiltonianReal(ParaRel,HopValClas)
WriteHamiltonianReal(ParaIn, ParaHmtR)


# '''###### Compute analytic forms of real space and k-space Hamiltonian #####'''
# # Write analytic form of Hamiltonian
# ParaHmtR   = ReadHamiltonianReal(ParaIn)
# ParaKptIn  = ReadKPoint(ParaIn["Folder"]+"/KptIN_" + Name + ".txt")
# ParaKpt    = GetKPoint(ParaKptIn)
# WriteHamiltonianAnalytic(ParaIn,ParaRel)#Hamiltonian in real space
# # Calculate energy bands from analytic form
# HkMatSp,k1,k2,k3,tValsSpAll=writeHkAnalytic(ParaIn,ParaRel,ParaHmtR)# H(k): Hamiltonian in k-space
# Hk2html(ParaIn,HkMatSp,ParaIn["Folder"])#write H(k) to html file
# pltSpEigs=plotEigsFromHkSp(HkMatSp,k1,k2,k3,tValsSpAll,ParaKpt,HopValIn,ParaRel,ParaIn["Folder"])
# # HkTmp=HkSp2Np(HkMatSp,k1,k2,k3,tValsSpAll,1,2,3,HopValIn)


'''######################## Calculate Energy Bands #########################'''
# From Hamiltonian in real space to Hamiltonian in k space
ParaHmtR   = ReadHamiltonianReal(ParaIn)
ParaKptIn  = ReadKPoint(ParaIn["Folder"]+"/KptIN_" + Name + ".txt")
ParaKpt    = GetKPoint(ParaKptIn)
ParaHmtK   = GetHamiltonianK(ParaHmtR,ParaKpt)
# Solve the k-space Hamiltonian to get energy bands
tEigStart=datetime.now()
ParaEig    = GetEigenSolution(ParaHmtK)
ParaEigPlt = PlotEnergyBand(ParaKpt,ParaEig,ParaIn["Folder"])
tEigEnd=datetime.now()
print(Name+ " enerygy band: ",tEigEnd-tEigStart)


'''######### Check symmetry of Energy Bands and k-space Hamiltonian ########'''
CheckE(ParaIn,ParaSym,ParaSymAt,ParaHmtR,0)
CheckH(ParaIn,ParaSym,ParaSymAt,ParaHmtR)


# '''####################### Calculate Berry Curvature #######################'''
# # Read real space Hamiltonian and k points
# ParaHmtR   = ReadHamiltonianReal(ParaIn)
# ParaKptIn  = ReadKPoint(ParaIn["Folder"]+"/KptIN_" + Name + ".txt")
# ParaKpt    = GetKPoint(ParaKptIn)
# # Calculate Berry Curvature
# ParaBerK = GetBerryK(ParaIn,ParaHmtR,ParaKpt)
# ParaBerPlt = PlotBerryCurvature(ParaKpt,ParaBerK,ParaIn["Folder"])
# '''This part needs an isolated test!'''


# '''########### Calculate intrinsic anomalous Hall conductivity and 
#                          intrinsic anomalous Nernst conductivity ###########'''
# # Read real space Hamiltonian and input parameters (Temperature, Fermi Energy, nk)
# ParaHmtR   = ReadHamiltonianReal(ParaIn)
# ParaBerIn  = ReadBerryPara(ParaIn["Folder"]+ "/BerIN_" + Name + ".txt")
# # Calculate intrinsic AHC & ANC
# ParaAHNC = GetAHNC(ParaIn,ParaHmtR,ParaBerIn)
# '''This part needs an isolated test!'''


'''###################### Calculate density of states ######################'''
# Read real space Hamiltonian and input parameters (method,nE,nk)
ParaHmtR   = ReadHamiltonianReal(ParaIn)
ParaDosIn  = ReadDosPara(ParaIn["Folder"]+ "/DosIN_" + Name + ".txt")
# Calculate DOS
ParaDos = GetDos(ParaIn,ParaHmtR,ParaDosIn)
ParaDosPlt = PlotDOS(ParaDos,ParaIn["Folder"])
