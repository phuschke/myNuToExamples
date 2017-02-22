import subprocess

matrixMeshFile              =   "3d_3_point_bending_cuboid.msh"
fiberMeshFile               =   "trusses.msh"
tensileStrength             =   "1.2"
compressiveStrength         =   "300"
fractureEnergy              =   "0.01"
nonlocalRadius              =   "0.2"
youngsModulusFiber          =   "1.75e1"
crossSectionFiber           =   "1.0"
circumferenceFiber          =   "0.628"
interfaceNormalStiffness    =   "1e6"
alpha                       =   "1"
maxBondStress               =   "3.e1"
residualBondStress          =   "1.e1"
slipAtMaxBondStress         =   "0.1"
slipAtResidualBondStress    =   "1"
subDirectory                =   "1337"
args = ("./3d_reference", matrixMeshFile, fiberMeshFile, tensileStrength, compressiveStrength, fractureEnergy, nonlocalRadius, youngsModulusFiber, crossSectionFiber, circumferenceFiber, interfaceNormalStiffness, alpha, maxBondStress, residualBondStress, slipAtMaxBondStress, slipAtResidualBondStress, subDirectory)

proc = subprocess.Popen([args], shell=True,stdin=None, stdout=None, stderr=None, close_fds=True)




