import sys
import numpy as np
import matplotlib.pyplot as plt
import os
import pathlib
from datetime import datetime
import sympy as smp

plt.close("all")
# lattice class
# from lattice.baseLattice import baseLattice
from lattice.superLattice import superLattice
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
from cd.Strip    import GetStrip
from cd.EigStrip import GetEigStrip
from cd.Qc       import GetQc
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
from rw.ReadQc   import ReadQcPara
# Plot
from pl.PlotEB   import PlotEnergyBand
from pl.PlotEB   import plotEigsFromHkSp
from pl.PlotBC   import PlotBerryCurvature
from pl.PlotDos  import PlotDOS
from pl.PlotStrip import PlotAtomStrip
from pl.PlotEB1D import PlotEnergyBand1D
from pl.PlotQc   import PlotQuanCond
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

material = 'data/ABO3/supercell_TBIN_ABO3.txt'
# material = 'data/Graphene/primitive_TBIN_Graphene.txt'
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
superCrystal=superLattice()
superCrystal.ParaIn=ReadInput(inConfigName)
superCrystal.checkSupercellInfoSanity()
superCrystal.constructSuperLattice()
# Name=ParaIn["Name"]


'''################### Determine space group of crystal ####################'''
# ParaSym = GetSpaceGroupPrimitive(ParaIn)
atmUnderConvVector,atmIndsConv=paraInPrim2Conv(superCrystal.ParaIn)
SGN, originBilbao=conv2SGN(superCrystal.ParaIn,atmUnderConvVector,atmIndsConv)
superCrystal.ParaSym=getParaSym(superCrystal.ParaIn,SGN, originBilbao)
superCrystal.ParaIn["origin Bilbao"]=superCrystal.ParaSym["origin Bilbao"]


'''##################### Complete Electronic orbitals ######################'''
# CheckOrbital(baseCrystal.ParaIn,baseCrystal.ParaSym)


'''################### Relate hopping terms by symmetry ####################'''
# baseCrystal.ParaNbr    = FindNeighbor(baseCrystal.ParaIn)
# PlotAtoms(baseCrystal.ParaIn,baseCrystal.ParaNbr,baseCrystal.ParaIn["Name"])
# tFindingRelationStart=datetime.now()
# baseCrystal.ParaSymAt  = FindAtomSymmetry(baseCrystal.ParaIn,baseCrystal.ParaSym,baseCrystal.ParaNbr)
# baseCrystal.ParaRel    = FindRelation(baseCrystal.ParaIn,baseCrystal.ParaSym,baseCrystal.ParaNbr,baseCrystal.ParaSymAt)
# PlotHoppingTerm(baseCrystal.ParaIn,baseCrystal.ParaNbr,baseCrystal.ParaRel,baseCrystal.ParaIn["Name"],[5,6])
# tFindingRelationEnd=datetime.now()
# print("Finding symmetry relations: ",tFindingRelationEnd-tFindingRelationStart)
# WriteRelation(baseCrystal.ParaIn,baseCrystal.ParaRel)


'''################### Calculate real space Hamiltonian ####################'''
# # From hopping relations to Hamiltonian in real space
# HopValIn   = ReadHopping(baseCrystal.ParaIn["Folder"]+"/HopValIN_"+baseCrystal.ParaIn["Name"]+".txt")
# # ParaRel    = ReadRelation(baseCrystal.ParaIn["Folder"]+ "/HopRel.txt")
# baseCrystal.HopValClas = GetHoppingValue(HopValIn, baseCrystal.ParaRel)
# ParaHmtR   = GetHamiltonianReal(baseCrystal.ParaRel,baseCrystal.HopValClas)
# WriteHamiltonianReal(baseCrystal.ParaIn, ParaHmtR)


# '''###### Compute analytic forms of real space and k-space Hamiltonian #####'''
# # Write analytic form of Hamiltonian
# ParaHmtR   = ReadHamiltonianReal(baseCrystal.ParaIn)
# ParaKptIn  = ReadKPoint(baseCrystal.ParaIn["Folder"]+"/KptIN_" + baseCrystal.ParaIn["Name"] + ".txt")
# ParaKpt    = GetKPoint(ParaKptIn)
# WriteHamiltonianAnalytic(baseCrystal.ParaIn,baseCrystal.ParaRel)#Hamiltonian in real space
# # Calculate energy bands from analytic form
# HkMatSp,k1,k2,k3,tValsSpAll=writeHkAnalytic(baseCrystal.ParaIn,baseCrystal.ParaRel,ParaHmtR)# H(k): Hamiltonian in k-space
# Hk2html(baseCrystal.ParaIn,HkMatSp,baseCrystal.ParaIn["Folder"])#write H(k) to html file
# pltSpEigs=plotEigsFromHkSp(HkMatSp,k1,k2,k3,tValsSpAll,ParaKpt,HopValIn,baseCrystal.ParaRel,baseCrystal.ParaIn["Folder"])
# HkTmp=HkSp2Np(HkMatSp,k1,k2,k3,tValsSpAll,1,2,3,HopValIn)


'''######################## Calculate Energy Bands #########################'''
# From Hamiltonian in real space to Hamiltonian in k space
# ParaHmtR   = ReadHamiltonianReal(baseCrystal.ParaIn)
# ParaKptIn  = ReadKPoint(baseCrystal.ParaIn["Folder"]+"/KptIN_" + baseCrystal.ParaIn["Name"]+ ".txt")
# ParaKpt    = GetKPoint(ParaKptIn)
# ParaHmtK   = GetHamiltonianK(ParaHmtR,ParaKpt)
# # Solve the k-space Hamiltonian to get energy bands
# tEigStart=datetime.now()
# ParaEig    = GetEigenSolution(ParaHmtK)
# ParaEigPlt = PlotEnergyBand(ParaKpt,ParaEig,baseCrystal.ParaIn["Folder"])
# tEigEnd=datetime.now()
# print(baseCrystal.ParaIn["Name"]+ " enerygy band: ",tEigEnd-tEigStart)


'''######### Check symmetry of Energy Bands and k-space Hamiltonian ########'''
# CheckE(baseCrystal.ParaIn,baseCrystal.ParaSym,baseCrystal.ParaSymAt,ParaHmtR,0)
# CheckH(baseCrystal.ParaIn,baseCrystal.ParaSym,baseCrystal.ParaSymAt,ParaHmtR)


# '''####################### Calculate Berry Curvature #######################'''
# # Read real space Hamiltonian and k points
# ParaHmtR   = ReadHamiltonianReal(baseCrystal.ParaIn)
# ParaKptIn  = ReadKPoint(baseCrystal.ParaIn["Folder"]+"/KptIN_" + baseCrystal.ParaIn["Name"] + ".txt")
# ParaKpt    = GetKPoint(ParaKptIn)
# # Calculate Berry Curvature
# ParaBerK = GetBerryK(baseCrystal.ParaIn,ParaHmtR,ParaKpt)
# ParaBerPlt = PlotBerryCurvature(ParaKpt,ParaBerK,baseCrystal.ParaIn["Folder"])
# '''This part needs an isolated test!'''


# '''########### Calculate intrinsic anomalous Hall conductivity and 
#                          intrinsic anomalous Nernst conductivity ###########'''
# # Read real space Hamiltonian and input parameters (Temperature, Fermi Energy, nk)
# ParaHmtR   = ReadHamiltonianReal(baseCrystal.ParaIn)
# ParaBerIn  = ReadBerryPara(baseCrystal.ParaIn["Folder"]+ "/BerIN_" + baseCrystal.ParaIn["Name"] + ".txt")
# # Calculate intrinsic AHC & ANC
# ParaAHNC = GetAHNC(baseCrystal.ParaIn,ParaHmtR,ParaBerIn)
# '''This part needs an isolated test!'''


# '''###################### Calculate density of states ######################'''
# # Read real space Hamiltonian and input parameters (method,nE,nk)
# ParaHmtR   = ReadHamiltonianReal(baseCrystal.ParaIn)
# ParaDosIn  = ReadDosPara(baseCrystal.ParaIn["Folder"]+ "/DosIN_" + baseCrystal.ParaIn["Name"] + ".txt")
# # Calculate DOS
# ParaDos = GetDos(baseCrystal.ParaIn,ParaHmtR,ParaDosIn)
# ParaDosPlt = PlotDOS(ParaDos,baseCrystal.ParaIn["Folder"])

'''######################## Calculate info of strip ########################'''
# Read real space Hamiltonian and input parameters (method,nE,nk)
# ParaHmtR   = ReadHamiltonianReal(baseCrystal.ParaIn)
# ParaQcIn   = ReadQcPara(baseCrystal.ParaIn["Folder"]+ "/QcIN_" + baseCrystal.ParaIn["Name"] + ".txt")
# # Info of strip
# ParaStrip  = GetStrip(baseCrystal.ParaIn,ParaHmtR,ParaQcIn)
# # Plot strip
# PlotAtomStrip(ParaStrip,baseCrystal.ParaIn["Name"])
# # Calculate energy band
# ParaEigStrip = GetEigStrip(ParaStrip)
# PlotEnergyBand1D(ParaEigStrip,baseCrystal.ParaIn["Folder"])
# # Calculate quantum conductance
# ParaQc     = GetQc(ParaQcIn,ParaStrip)
# PlotQuanCond(ParaQc,baseCrystal.ParaIn["Folder"])
'''This part is only tested with graphene & BN!'''