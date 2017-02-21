import subprocess
import matplotlib.pyplot as plt

matrixMeshFile              =   "3d_3_point_bending_cuboid.msh"
fiberMeshFile               =   "trusses.msh"
tensileStrength             =   "3"
compressiveStrength         =   "300"
fractureEnergy              =   "0.01"
nonlocalRadius              =   "0.2"
youngsModulusFiber          =   "1.75e4"
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

popen = subprocess.Popen(args, stdout=subprocess.PIPE)
popen.wait()
output = popen.stdout.read()
print output


path = "/home/phuschke/results/results_3_point_bending_reference/" + subDirectory
with open(path + '/displacements.dat') as f:
    content = f.readlines()

content = [x.strip() for x in content]
content = [x.split() for x in content]
disp = [x[1] for x in content]


with open(path + '/force.dat') as f:
    content = f.readlines()

content = [x.strip() for x in content]
content = [x.split() for x in content]
force = [x[1] for x in content]

print(content[1][0])

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.plot(disp, force,'-bo')
plt.xlabel('displacement')
plt.ylabel('force')
#plt.show()
plt.savefig(path + '/forceDisp.png')



########################################################################################################################

youngsModulusFiber          =   "2.0e4"
subDirectory                =   "1"
args = ("./3d_reference", matrixMeshFile, fiberMeshFile, tensileStrength, compressiveStrength, fractureEnergy, nonlocalRadius, youngsModulusFiber, crossSectionFiber, circumferenceFiber, interfaceNormalStiffness, alpha, maxBondStress, residualBondStress, slipAtMaxBondStress, slipAtResidualBondStress, subDirectory)

popen = subprocess.Popen(args, stdout=subprocess.PIPE)
popen.wait()
output = popen.stdout.read()
print output


path = "/home/phuschke/results/results_3_point_bending_reference/" + subDirectory
with open(path + '/displacements.dat') as f:
    content = f.readlines()

content = [x.strip() for x in content]
content = [x.split() for x in content]
disp = [x[1] for x in content]


with open(path + '/force.dat') as f:
    content = f.readlines()

content = [x.strip() for x in content]
content = [x.split() for x in content]
force = [x[1] for x in content]

print(content[1][0])

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.plot(disp, force,'-bo')
plt.xlabel('displacement')
plt.ylabel('force')
#plt.show()
plt.savefig(path + '/forceDisp.png')

########################################################################################################################

youngsModulusFiber          =   "1.5e4"
subDirectory                =   "2"
args = ("./3d_reference", matrixMeshFile, fiberMeshFile, tensileStrength, compressiveStrength, fractureEnergy, nonlocalRadius, youngsModulusFiber, crossSectionFiber, circumferenceFiber, interfaceNormalStiffness, alpha, maxBondStress, residualBondStress, slipAtMaxBondStress, slipAtResidualBondStress, subDirectory)

popen = subprocess.Popen(args, stdout=subprocess.PIPE)
popen.wait()
output = popen.stdout.read()
print output


path = "/home/phuschke/results/results_3_point_bending_reference/" + subDirectory
with open(path + '/displacements.dat') as f:
    content = f.readlines()

content = [x.strip() for x in content]
content = [x.split() for x in content]
disp = [x[1] for x in content]


with open(path + '/force.dat') as f:
    content = f.readlines()

content = [x.strip() for x in content]
content = [x.split() for x in content]
force = [x[1] for x in content]

print(content[1][0])

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.plot(disp, force,'-bo')
plt.xlabel('displacement')
plt.ylabel('force')
#plt.show()
plt.savefig(path + '/forceDisp.png')
