#include "nuto/math/FullMatrix.h"
#include "nuto/mechanics/structures/unstructured/Structure.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStrain3D.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStress3D.h"
#include "nuto/mechanics/constitutive/mechanics/GradientDamageEngineeringStress.h"
#include "nuto/mechanics/MechanicsException.h"
#include "nuto/mechanics/constitutive/ConstitutiveStaticDataBase.h"
#include "nuto/mechanics/constitutive/mechanics/ConstitutiveStaticDataGradientDamage1D.h"

#include "nuto/math/SparseMatrixCSRVector2General.h"
#include "nuto/math/SparseDirectSolverMUMPS.h"

#include "nuto/mechanics/timeIntegration/NewmarkDirect.h"
#include "nuto/base/Logger.h"

#include <boost-1_55/boost/filesystem.hpp>
#include <boost-1_55/boost/lexical_cast.hpp>
#include <boost-1_55/boost/progress.hpp>

#include <iostream>
#include <fstream>
#include <string>

const unsigned int dimension = 2;

class Parameters
{
public:

    static const int mDimension = dimension;

    static const bool mPerformLineSearch = true;
    static const bool mAutomaticTimeStepping = true;

    static constexpr double mMatrixYoungsModulus = 4.0e5;   // concrete
    static constexpr double mMatrixPoissonsRatio = 0.2;
    static constexpr double mMatrixThickness = 1;

    static constexpr double mBondYoungsModulus = 1.0e3;   // bond
    static constexpr double mBondPoissonsRatio = 0.2;
    static constexpr double mBondThickness = 10;

    static constexpr double mFibreYoungsModulus = 2.1e6;   // steel
    static constexpr double mFibrePoissonsRatio = 0.2;
    static constexpr double mFibreCrossSection = 0.01;

    static constexpr double mSpringStiffnessSoft = 1.0e5;   // spring soft
    static constexpr double mSpringStiffnessRigid = 1.0e7;   // spring rigid

    static constexpr double mNonlocalRadius = 4;
    static constexpr double mTensileStrength = 3;
    static constexpr double mCompressiveStrength = 30;
    static constexpr double mFractureEnergy = 0.01;
    static constexpr double mTimeStep = 1e-3;
    static constexpr double mMinTimeStep = 1e-5;
    static constexpr double mMaxTimeStep = 1e-1;
    static constexpr double mToleranceForce = 1e-6;
    static constexpr double mSimulationTime = 1.0;
    static constexpr double mLoad = 1;

    static const boost::filesystem::path mOutputPath;
    static const boost::filesystem::path mMeshFilePath;

    static const NuTo::FullVector<double, dimension> mDirectionX;
    static const NuTo::FullVector<double, dimension> mDirectionY;
};

const boost::filesystem::path Parameters::mOutputPath("/home/phuschke/2d_bond_stress_slip/");
const boost::filesystem::path Parameters::mMeshFilePath("/home/phuschke/develop/nuto/myNutoExamples/2d_bond_stress_slip.msh");

const NuTo::FullVector<double, dimension> Parameters::mDirectionX = NuTo::FullVector<double, 2>::UnitX();
const NuTo::FullVector<double, dimension> Parameters::mDirectionY = NuTo::FullVector<double, 2>::UnitY();

int main(int argc, char* argv[])
{

    std::cout << "**********************************************" << std::endl;
    std::cout << "                  Structure                   " << std::endl;
    std::cout << "**********************************************" << std::endl;

    NuTo::Structure myStructure(Parameters::mDimension);
    myStructure.SetVerboseLevel(10);
    myStructure.SetShowTime(false);

    std::cout << "**********************************************" << std::endl;
    std::cout << "                  Integration Scheme          " << std::endl;
    std::cout << "**********************************************" << std::endl;

    NuTo::NewmarkDirect myIntegrationScheme(&myStructure);
    myIntegrationScheme.SetTimeStep(Parameters::mTimeStep);
    myIntegrationScheme.SetMinTimeStep(Parameters::mMinTimeStep);
    myIntegrationScheme.SetMaxTimeStep(Parameters::mMaxTimeStep);
    myIntegrationScheme.SetToleranceForce(Parameters::mToleranceForce);
    myIntegrationScheme.SetAutomaticTimeStepping(Parameters::mAutomaticTimeStepping);
    myIntegrationScheme.SetPerformLineSearch(Parameters::mPerformLineSearch);

    std::cout << "**********************************************" << std::endl;
    std::cout << "                  Section                     " << std::endl;
    std::cout << "**********************************************" << std::endl;

    int springSection = myStructure.SectionCreate(NuTo::Section::SPRING);

    int matrixSection = myStructure.SectionCreate(NuTo::Section::PLANE_STRESS);
    myStructure.SectionSetThickness(matrixSection, Parameters::mMatrixThickness);

    int fibreSection = myStructure.SectionCreate(NuTo::Section::TRUSS);
    myStructure.SectionSetArea(fibreSection, Parameters::mFibreCrossSection);

    std::cout << "**********************************************" << std::endl;
    std::cout << "                  Material                    " << std::endl;
    std::cout << "**********************************************" << std::endl;

    int springMaterialRigid = myStructure.ConstitutiveLawCreate(NuTo::Constitutive::eConstitutiveType::LINEAR_SPRING);
    myStructure.ConstitutiveLawSetParameterDouble(springMaterialRigid, NuTo::Constitutive::eConstitutiveParameter::SPRING_STIFFNESS, Parameters::mSpringStiffnessRigid);
    myStructure.ConstitutiveLawSetParameterFullVectorDouble(springMaterialRigid, NuTo::Constitutive::eConstitutiveParameter::SPRING_DIRECTION, Parameters::mDirectionY);

    int springMaterialSoft = myStructure.ConstitutiveLawCreate(NuTo::Constitutive::eConstitutiveType::LINEAR_SPRING);
    myStructure.ConstitutiveLawSetParameterDouble(springMaterialSoft, NuTo::Constitutive::eConstitutiveParameter::SPRING_STIFFNESS, Parameters::mSpringStiffnessSoft);
    myStructure.ConstitutiveLawSetParameterFullVectorDouble(springMaterialSoft, NuTo::Constitutive::eConstitutiveParameter::SPRING_DIRECTION, Parameters::mDirectionX);

    int matrixMaterial = myStructure.ConstitutiveLawCreate(NuTo::Constitutive::eConstitutiveType::LINEAR_ELASTIC_ENGINEERING_STRESS);
    myStructure.ConstitutiveLawSetParameterDouble(matrixMaterial, NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, Parameters::mMatrixYoungsModulus);
    myStructure.ConstitutiveLawSetParameterDouble(matrixMaterial, NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO, Parameters::mMatrixPoissonsRatio);

    int fibreMaterial = myStructure.ConstitutiveLawCreate(NuTo::Constitutive::eConstitutiveType::LINEAR_ELASTIC_ENGINEERING_STRESS);
    myStructure.ConstitutiveLawSetParameterDouble(fibreMaterial, NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, Parameters::mFibreYoungsModulus);
    myStructure.ConstitutiveLawSetParameterDouble(fibreMaterial, NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO, Parameters::mFibrePoissonsRatio);

    std::cout << "**********************************************" << std::endl;
    std::cout << "                  Interpolation Type          " << std::endl;
    std::cout << "**********************************************" << std::endl;

    int matrixInterpolationType = myStructure.InterpolationTypeCreate(NuTo::Interpolation::eShapeType::QUAD2D);
    myStructure.InterpolationTypeAdd(matrixInterpolationType, NuTo::Node::COORDINATES, NuTo::Interpolation::EQUIDISTANT1);
    myStructure.InterpolationTypeAdd(matrixInterpolationType, NuTo::Node::DISPLACEMENTS, NuTo::Interpolation::EQUIDISTANT1);

    int fibreInterpolationType = myStructure.InterpolationTypeCreate(NuTo::Interpolation::eShapeType::TRUSSXD);
    myStructure.InterpolationTypeAdd(fibreInterpolationType, NuTo::Node::COORDINATES, NuTo::Interpolation::EQUIDISTANT1);
    myStructure.InterpolationTypeAdd(fibreInterpolationType, NuTo::Node::DISPLACEMENTS, NuTo::Interpolation::EQUIDISTANT1);

    int springInterpolationType = myStructure.InterpolationTypeCreate(NuTo::Interpolation::eShapeType::SPRING);
    myStructure.InterpolationTypeAdd(springInterpolationType, NuTo::Node::COORDINATES, NuTo::Interpolation::EQUIDISTANT1);
    myStructure.InterpolationTypeAdd(springInterpolationType, NuTo::Node::DISPLACEMENTS, NuTo::Interpolation::EQUIDISTANT1);

    std::cout << "**********************************************" << std::endl;
    std::cout << "                  Nodes                       " << std::endl;
    std::cout << "**********************************************" << std::endl;

    NuTo::FullVector<double, Eigen::Dynamic> nodeCoords(Parameters::mDimension);
    std::set<NuTo::Node::eAttributes> dofs;
    dofs.insert(NuTo::Node::COORDINATES);
    dofs.insert(NuTo::Node::DISPLACEMENTS);

    nodeCoords[0] = 0.0;
    nodeCoords[1] = 0.0;
    int nodeMatrix00 = myStructure.NodeCreate(nodeCoords, dofs);

    nodeCoords[0] = 1.0;
    nodeCoords[1] = 0.0;
    int nodeMatrix01 = myStructure.NodeCreate(nodeCoords, dofs);

    nodeCoords[0] = 1.0;
    nodeCoords[1] = 1.0;
    int nodeMatrix02 = myStructure.NodeCreate(nodeCoords, dofs);

    nodeCoords[0] = 0.0;
    nodeCoords[1] = 1.0;
    int nodeMatrix03 = myStructure.NodeCreate(nodeCoords, dofs);

    nodeCoords[0] = 0.0;
    nodeCoords[1] = -1.0;
    int nodeFibre00 = myStructure.NodeCreate(nodeCoords, dofs);

    nodeCoords[0] = 1.0;
    nodeCoords[1] = -1.0;
    int nodeFibre01 = myStructure.NodeCreate(nodeCoords, dofs);

    std::cout << "**********************************************" << std::endl;
    std::cout << "                  Matrix Geometry             " << std::endl;
    std::cout << "**********************************************" << std::endl;

    NuTo::FullVector<int, Eigen::Dynamic> nodeIndicesMatrix(4);

    nodeIndicesMatrix[0] = nodeMatrix00;
    nodeIndicesMatrix[1] = nodeMatrix01;
    nodeIndicesMatrix[2] = nodeMatrix02;
    nodeIndicesMatrix[3] = nodeMatrix03;
    int elementMatirx00 = myStructure.ElementCreate(matrixInterpolationType, nodeIndicesMatrix, NuTo::ElementData::eElementDataType::CONSTITUTIVELAWIP, NuTo::IpData::eIpDataType::STATICDATA);

    myStructure.ElementSetSection(elementMatirx00, matrixSection);
    myStructure.ElementSetConstitutiveLaw(elementMatirx00, matrixMaterial);

    std::cout << "**********************************************" << std::endl;
    std::cout << "                  Fibre Geometry              " << std::endl;
    std::cout << "**********************************************" << std::endl;

    NuTo::FullVector<int, Eigen::Dynamic> nodeIndicesFibre(2);

    nodeIndicesFibre[0] = nodeFibre00;
    nodeIndicesFibre[1] = nodeFibre01;
    int elementFibre00 = myStructure.ElementCreate(fibreInterpolationType, nodeIndicesFibre, NuTo::ElementData::eElementDataType::CONSTITUTIVELAWIP, NuTo::IpData::eIpDataType::NOIPDATA);

    myStructure.ElementSetSection(elementFibre00, fibreSection);
    myStructure.ElementSetConstitutiveLaw(elementFibre00, fibreMaterial);

    myStructure.NodeBuildGlobalDofs();
    myStructure.Info();
    myStructure.CheckCoefficientMatrix_0(1e-6, true);

    std::cout << "**********************************************" << std::endl;
    std::cout << "                  Spring Geometry Rigid       " << std::endl;
    std::cout << "**********************************************" << std::endl;

    NuTo::FullVector<int, Eigen::Dynamic> nodeIndicesSpring(2);

    // springs rigid
    nodeIndicesSpring[0] = nodeFibre00;
    nodeIndicesSpring[1] = nodeMatrix00;
    int elementSpring00 = myStructure.ElementCreate(springInterpolationType, nodeIndicesSpring, NuTo::ElementData::eElementDataType::CONSTITUTIVELAWIP, NuTo::IpData::eIpDataType::NOIPDATA);

    nodeIndicesSpring[0] = nodeFibre01;
    nodeIndicesSpring[1] = nodeMatrix01;
    int elementSpring01 = myStructure.ElementCreate(springInterpolationType, nodeIndicesSpring, NuTo::ElementData::eElementDataType::CONSTITUTIVELAWIP, NuTo::IpData::eIpDataType::NOIPDATA);

    int groupElementSpringRigid = myStructure.GroupCreate(NuTo::Groups::eGroupId::Elements);
    myStructure.GroupAddElement(groupElementSpringRigid, elementSpring00);
    myStructure.GroupAddElement(groupElementSpringRigid, elementSpring01);

    myStructure.ElementGroupSetSection(groupElementSpringRigid, springSection);
    myStructure.ElementGroupSetConstitutiveLaw(groupElementSpringRigid, springMaterialRigid);

    std::cout << "**********************************************" << std::endl;
    std::cout << "                  Spring Geometry Soft        " << std::endl;
    std::cout << "**********************************************" << std::endl;

    // springs soft
    nodeIndicesSpring[0] = nodeFibre00;
    nodeIndicesSpring[1] = nodeMatrix00;
    int elementSpring02 = myStructure.ElementCreate(springInterpolationType, nodeIndicesSpring, NuTo::ElementData::eElementDataType::CONSTITUTIVELAWIP, NuTo::IpData::eIpDataType::NOIPDATA);

    nodeIndicesSpring[0] = nodeFibre01;
    nodeIndicesSpring[1] = nodeMatrix01;
    int elementSpring03 = myStructure.ElementCreate(springInterpolationType, nodeIndicesSpring, NuTo::ElementData::eElementDataType::CONSTITUTIVELAWIP, NuTo::IpData::eIpDataType::NOIPDATA);

    int groupElementSpringSoft = myStructure.GroupCreate(NuTo::Groups::eGroupId::Elements);
    myStructure.GroupAddElement(groupElementSpringSoft, elementSpring02);
    myStructure.GroupAddElement(groupElementSpringSoft, elementSpring03);

    myStructure.ElementGroupSetSection(groupElementSpringSoft, springSection);
    myStructure.ElementGroupSetConstitutiveLaw(groupElementSpringSoft, springMaterialSoft);

    std::cout << "**********************************************" << std::endl;
    std::cout << "                  Boundary Conditions         " << std::endl;
    std::cout << "**********************************************" << std::endl;

    int groupNodeBCLeft = myStructure.GroupCreate(NuTo::Groups::eGroupId::Nodes);
    myStructure.GroupAddNodeCoordinateRange(groupNodeBCLeft, 0, -1e-6, +1e-6);

    int groupNodeBCBottom = myStructure.GroupCreate(NuTo::Groups::eGroupId::Nodes);
    myStructure.GroupAddNodeCoordinateRange(groupNodeBCBottom, 1, -1-1e-6, -1+1e-6);

    myStructure.ConstraintLinearSetDisplacementNodeGroup(groupNodeBCLeft, Parameters::mDirectionX, 0);
    myStructure.ConstraintLinearSetDisplacementNodeGroup(groupNodeBCBottom, Parameters::mDirectionY, 0);

    std::cout << "**********************************************" << std::endl;
    std::cout << "                  Loads                       " << std::endl;
    std::cout << "**********************************************" << std::endl;

    myStructure.LoadCreateNodeForce(0, nodeFibre01, Parameters::mDirectionX, 1e4);

    std::cout << "**********************************************" << std::endl;
    std::cout << "                  Visualisation               " << std::endl;
    std::cout << "**********************************************" << std::endl;

    myStructure.AddVisualizationComponentDisplacements();
    myStructure.AddVisualizationComponentSection();
    myStructure.AddVisualizationComponentConstitutive();
    myStructure.AddVisualizationComponentElement();

    std::cout << "**********************************************" << std::endl;
    std::cout << "                  Solver                      " << std::endl;
    std::cout << "**********************************************" << std::endl;

    myStructure.NodeBuildGlobalDofs();
    myStructure.CheckCoefficientMatrix_0(1e-6, true);
    myStructure.Info();

    // build global stiffness matrix and equivalent load vector which correspond to prescribed boundary values
    NuTo::SparseMatrixCSRVector2General<double> stiffnessMatrixVec;
    NuTo::FullVector<double, Eigen::Dynamic> dispForceVector;
    myStructure.CalculateMaximumIndependentSets();
    myStructure.BuildGlobalCoefficientMatrix0(stiffnessMatrixVec, dispForceVector);
    NuTo::SparseMatrixCSRGeneral<double> stiffnessMatrix(stiffnessMatrixVec);

    std::cout << "stiffnessMatrix.Info() \n" << std::endl;
    stiffnessMatrix.Info();
    stiffnessMatrix += stiffnessMatrix;
    std::cout << "stiffnessMatrix.Info() \n" << std::endl;

    stiffnessMatrix.Info();

    // build global external load vector
    NuTo::FullVector<double, Eigen::Dynamic> extForceVector;
    myStructure.BuildGlobalExternalLoadVector(0, extForceVector);

    std::cout << "extForceVector \n" << extForceVector << std::endl;
    std::cout << "dispForceVector \n" << dispForceVector << std::endl;
    // calculate right hand side
    NuTo::FullVector<double, Eigen::Dynamic> rhsVector = dispForceVector + extForceVector;

    std::cout << "rhsVector \n" << rhsVector << std::endl;
    std::cout << "?????????????????????????????????????????????????????????????????" << std::endl;
    std::cout << "stiffnessMatrix.Info() \n" << std::endl;
    stiffnessMatrix.Info();
    // solve
    NuTo::SparseDirectSolverMUMPS mySolver;
    NuTo::FullVector<double, Eigen::Dynamic> displacementVector;
    stiffnessMatrix.SetOneBasedIndexing();
    mySolver.Solve(stiffnessMatrix, rhsVector, displacementVector);

    // write displacements to node
    myStructure.NodeMergeActiveDofValues(displacementVector);

    // calculate residual
    NuTo::FullVector<double, Eigen::Dynamic> intForceVector;
    myStructure.BuildGlobalGradientInternalPotentialVector(intForceVector);
    NuTo::FullVector<double, Eigen::Dynamic> residualVector = extForceVector - intForceVector;
    std::cout << "residual: " << residualVector.Norm() << std::endl;

    myStructure.ExportVtkDataFileElements(Parameters::mOutputPath.string() + "Truss2D.vtk");

    myStructure.InterpolationTypeInfo(matrixInterpolationType);
    myStructure.IntegrationTypeInfo(10);

    myStructure.CheckCoefficientMatrix_0(1e-6,true);
    std::cout << "**********************************************" << std::endl;
    std::cout << "                  End                         " << std::endl;
    std::cout << "**********************************************" << std::endl;
}

