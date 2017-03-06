
matrixMeshFile="3d_3_point_bending_cuboid.msh"
fiberMeshFile="trusses.msh"
tensileStrength="1.2"
compressiveStrength="300"
fractureEnergy="0.05"
nonlocalRadius="0.02"
youngsModulusFiber="1.75e1"
crossSectionFiber="1.0"
circumferenceFiber="0.628"
interfaceNormalStiffness="1e6"
alpha="1"
maxBondStress="3.e1"
residualBondStress="1.e1"
slipAtMaxBondStress="0.1"
slipAtResidualBondStress="1"
subDirectory="15"
writeOutputToFile="15_bla.txt"

./3d_reference $matrixMeshFile $fiberMeshFile $tensileStrength $compressiveStrength $fractureEnergy $nonlocalRadius $youngsModulusFiber $crossSectionFiber $circumferenceFiber $interfaceNormalStiffness $alpha $maxBondStress $residualBondStress $slipAtMaxBondStress $slipAtResidualBondStress $subDirectory > $writeOutputToFile 2>&1 &

fractureEnergy="0.02"
nonlocalRadius="0.02"
subDirectory="16"
writeOutputToFile="16_bla.txt"

./3d_reference $matrixMeshFile $fiberMeshFile $tensileStrength $compressiveStrength $fractureEnergy $nonlocalRadius $youngsModulusFiber $crossSectionFiber $circumferenceFiber $interfaceNormalStiffness $alpha $maxBondStress $residualBondStress $slipAtMaxBondStress $slipAtResidualBondStress $subDirectory > $writeOutputToFile 2>&1 &

fractureEnergy="0.05"
nonlocalRadius="0.01"
subDirectory="17"
writeOutputToFile="17_bla.txt"

./3d_reference $matrixMeshFile $fiberMeshFile $tensileStrength $compressiveStrength $fractureEnergy $nonlocalRadius $youngsModulusFiber $crossSectionFiber $circumferenceFiber $interfaceNormalStiffness $alpha $maxBondStress $residualBondStress $slipAtMaxBondStress $slipAtResidualBondStress $subDirectory > $writeOutputToFile 2>&1 &

fractureEnergy="0.02"
nonlocalRadius="0.01"
subDirectory="18"
writeOutputToFile="18_bla.txt"

./3d_reference $matrixMeshFile $fiberMeshFile $tensileStrength $compressiveStrength $fractureEnergy $nonlocalRadius $youngsModulusFiber $crossSectionFiber $circumferenceFiber $interfaceNormalStiffness $alpha $maxBondStress $residualBondStress $slipAtMaxBondStress $slipAtResidualBondStress $subDirectory > $writeOutputToFile 2>&1 &


fractureEnergy="0.01"
nonlocalRadius="0.01"
subDirectory="19"
writeOutputToFile="19_bla.txt"

./3d_reference $matrixMeshFile $fiberMeshFile $tensileStrength $compressiveStrength $fractureEnergy $nonlocalRadius $youngsModulusFiber $crossSectionFiber $circumferenceFiber $interfaceNormalStiffness $alpha $maxBondStress $residualBondStress $slipAtMaxBondStress $slipAtResidualBondStress $subDirectory > $writeOutputToFile 2>&1 &
