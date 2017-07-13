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

    static constexpr double mMatrixYoungsModulus = 4.0e4; // concrete
    static constexpr double mMatrixPoissonsRatio = 0.2;
    static constexpr double mMatrixThickness = 10;

    static constexpr double mBondYoungsModulus = 1.0e3; // bond
    static constexpr double mBondPoissonsRatio = 0.2;
    static constexpr double mBondThickness = 10;

    static constexpr double mFibreYoungsModulus = 1.0e3; // steel
    static constexpr double mFibrePoissonsRatio = 0.2;
    static constexpr double mFibreCrossSection = 0.01;

    static constexpr double mSpringStiffnessSoft = 1.0e0; // spring soft
    static constexpr double mSpringStiffnessRigid = 1.0e0; // spring rigid

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

const boost::filesystem::path Parameters::mOutputPath("/home/phuschke/2d_bond_stress_slip_2_trusses/");
const boost::filesystem::path
        Parameters::mMeshFilePath("/home/phuschke/develop/nuto/myNutoExamples/2d_bond_stress_slip.msh");

const NuTo::FullVector<double, dimension> Parameters::mDirectionX = NuTo::FullVector<double, 2>::UnitX();
const NuTo::FullVector<double, dimension> Parameters::mDirectionY = NuTo::FullVector<double, 2>::UnitY();

int main(int argc, char* argv[])
{
    boost::filesystem::remove_all(Parameters::mOutputPath);
    boost::filesystem::create_directory(Parameters::mOutputPath);

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
    myStructure.ConstitutiveLawSetParameterDouble(springMaterialRigid,
                                                  NuTo::Constitutive::eConstitutiveParameter::SPRING_STIFFNESS,
                                                  Parameters::mSpringStiffnessRigid);
    myStructure.ConstitutiveLawSetParameterFullVectorDouble(
            springMaterialRigid, NuTo::Constitutive::eConstitutiveParameter::SPRING_DIRECTION, Parameters::mDirectionY);

    int springMaterialSoft = myStructure.ConstitutiveLawCreate(NuTo::Constitutive::eConstitutiveType::LINEAR_SPRING);
    myStructure.ConstitutiveLawSetParameterDouble(springMaterialSoft,
                                                  NuTo::Constitutive::eConstitutiveParameter::SPRING_STIFFNESS,
                                                  Parameters::mSpringStiffnessSoft);
    myStructure.ConstitutiveLawSetParameterFullVectorDouble(
            springMaterialSoft, NuTo::Constitutive::eConstitutiveParameter::SPRING_DIRECTION, Parameters::mDirectionX);

    int matrixMaterial =
            myStructure.ConstitutiveLawCreate(NuTo::Constitutive::eConstitutiveType::LINEAR_ELASTIC_ENGINEERING_STRESS);
    myStructure.ConstitutiveLawSetParameterDouble(matrixMaterial,
                                                  NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS,
                                                  Parameters::mMatrixYoungsModulus);
    myStructure.ConstitutiveLawSetParameterDouble(matrixMaterial,
                                                  NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO,
                                                  Parameters::mMatrixPoissonsRatio);

    int fibreMaterial =
            myStructure.ConstitutiveLawCreate(NuTo::Constitutive::eConstitutiveType::LINEAR_ELASTIC_ENGINEERING_STRESS);
    myStructure.ConstitutiveLawSetParameterDouble(
            fibreMaterial, NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, Parameters::mFibreYoungsModulus);
    myStructure.ConstitutiveLawSetParameterDouble(
            fibreMaterial, NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO, Parameters::mFibrePoissonsRatio);

    std::cout << "**********************************************" << std::endl;
    std::cout << "                  Interpolation Type          " << std::endl;
    std::cout << "**********************************************" << std::endl;

    int springInterpolationType = myStructure.InterpolationTypeCreate(NuTo::Interpolation::eShapeType::SPRING);
    myStructure.InterpolationTypeAdd(springInterpolationType, NuTo::Node::COORDINATES,
                                     NuTo::Interpolation::EQUIDISTANT1);
    myStructure.InterpolationTypeAdd(springInterpolationType, NuTo::Node::DISPLACEMENTS,
                                     NuTo::Interpolation::EQUIDISTANT1);

    int fibreInterpolationType = myStructure.InterpolationTypeCreate(NuTo::Interpolation::eShapeType::TRUSSXD);
    myStructure.InterpolationTypeAdd(fibreInterpolationType, NuTo::Node::COORDINATES,
                                     NuTo::Interpolation::EQUIDISTANT1);
    myStructure.InterpolationTypeAdd(fibreInterpolationType, NuTo::Node::DISPLACEMENTS,
                                     NuTo::Interpolation::EQUIDISTANT1);
    myStructure.InterpolationTypeSetIntegrationType(
            fibreInterpolationType, NuTo::IntegrationType::IntegrationType1D2NGauss2Ip, NuTo::IpData::STATICDATA);

    std::cout << "**********************************************" << std::endl;
    std::cout << "                  Fibre Geometry              " << std::endl;
    std::cout << "**********************************************" << std::endl;

    NuTo::FullVector<double, Eigen::Dynamic> nodeCoords(Parameters::mDimension);
    NuTo::FullVector<int, Eigen::Dynamic> nodeIndicesFibre(2);

    std::set<NuTo::Node::eAttributes> dofs;
    dofs.insert(NuTo::Node::COORDINATES);
    dofs.insert(NuTo::Node::DISPLACEMENTS);

    nodeCoords[0] = 0.0;
    nodeCoords[1] = 0.0;
    int node00 = myStructure.NodeCreate(nodeCoords, dofs);

    nodeCoords[0] = 10.0;
    nodeCoords[1] = 0.0;
    int node01 = myStructure.NodeCreate(nodeCoords, dofs);

    nodeCoords[0] = 0.0;
    nodeCoords[1] = 0.0;
    int node02 = myStructure.NodeCreate(nodeCoords, dofs);

    nodeCoords[0] = 10.0;
    nodeCoords[1] = 0.0;
    int node03 = myStructure.NodeCreate(nodeCoords, dofs);

    // trusses
    nodeIndicesFibre[0] = node00;
    nodeIndicesFibre[1] = node01;
    int elementFibre00 = myStructure.ElementCreate(fibreInterpolationType, nodeIndicesFibre,
                                                   NuTo::ElementData::eElementDataType::CONSTITUTIVELAWIP,
                                                   NuTo::IpData::eIpDataType::NOIPDATA);

    nodeIndicesFibre[0] = node02;
    nodeIndicesFibre[1] = node03;
    int elementFibre01 = myStructure.ElementCreate(fibreInterpolationType, nodeIndicesFibre,
                                                   NuTo::ElementData::eElementDataType::CONSTITUTIVELAWIP,
                                                   NuTo::IpData::eIpDataType::NOIPDATA);

    // springs rigid
    nodeIndicesFibre[0] = node00;
    nodeIndicesFibre[1] = node02;
    int elementSpring00 = myStructure.ElementCreate(springInterpolationType, nodeIndicesFibre,
                                                    NuTo::ElementData::eElementDataType::CONSTITUTIVELAWIP,
                                                    NuTo::IpData::eIpDataType::NOIPDATA);

    nodeIndicesFibre[0] = node01;
    nodeIndicesFibre[1] = node03;
    int elementSpring01 = myStructure.ElementCreate(springInterpolationType, nodeIndicesFibre,
                                                    NuTo::ElementData::eElementDataType::CONSTITUTIVELAWIP,
                                                    NuTo::IpData::eIpDataType::NOIPDATA);

    // springs rigid
    nodeIndicesFibre[0] = node00;
    nodeIndicesFibre[1] = node02;
    int elementSpring02 = myStructure.ElementCreate(springInterpolationType, nodeIndicesFibre,
                                                    NuTo::ElementData::eElementDataType::CONSTITUTIVELAWIP,
                                                    NuTo::IpData::eIpDataType::NOIPDATA);

    nodeIndicesFibre[0] = node01;
    nodeIndicesFibre[1] = node03;
    int elementSpring03 = myStructure.ElementCreate(springInterpolationType, nodeIndicesFibre,
                                                    NuTo::ElementData::eElementDataType::CONSTITUTIVELAWIP,
                                                    NuTo::IpData::eIpDataType::NOIPDATA);

    int groupElementFibre = myStructure.GroupCreate(NuTo::Groups::eGroupId::Elements);
    myStructure.GroupAddElement(groupElementFibre, elementFibre00);
    myStructure.GroupAddElement(groupElementFibre, elementFibre01);
    myStructure.ElementGroupSetSection(groupElementFibre, fibreSection);
    myStructure.ElementGroupSetConstitutiveLaw(groupElementFibre, fibreMaterial);

    int groupElementSpringRigid = myStructure.GroupCreate(NuTo::Groups::eGroupId::Elements);
    myStructure.GroupAddElement(groupElementSpringRigid, elementSpring00);
    myStructure.GroupAddElement(groupElementSpringRigid, elementSpring01);
    myStructure.ElementGroupSetSection(groupElementSpringRigid, springSection);
    myStructure.ElementGroupSetConstitutiveLaw(groupElementSpringRigid, springMaterialRigid);

    int groupElementSpringSoft = myStructure.GroupCreate(NuTo::Groups::eGroupId::Elements);
    myStructure.GroupAddElement(groupElementSpringSoft, elementSpring02);
    myStructure.GroupAddElement(groupElementSpringSoft, elementSpring03);
    myStructure.ElementGroupSetSection(groupElementSpringSoft, springSection);
    myStructure.ElementGroupSetConstitutiveLaw(groupElementSpringSoft, springMaterialSoft);

    myStructure.Info();
    myStructure.NodeBuildGlobalDofs();

    //    myStructure.CheckCoefficientMatrix_0(1e-6, true);

    std::cout << "**********************************************" << std::endl;
    std::cout << "                  Boundary Conditions         " << std::endl;
    std::cout << "**********************************************" << std::endl;

    myStructure.ConstraintLinearSetDisplacementNode(node00, Parameters::mDirectionX, 0);
    myStructure.ConstraintLinearSetDisplacementNode(node00, Parameters::mDirectionY, 0);
    myStructure.ConstraintLinearSetDisplacementNode(node01, Parameters::mDirectionY, 0);

    std::cout << "**********************************************" << std::endl;
    std::cout << "                  Loads                       " << std::endl;
    std::cout << "**********************************************" << std::endl;
    NuTo::FullVector<double, dimension> loadDirection = NuTo::FullVector<double, 2>::Ones();
    int loadCase = 0;
    myStructure.LoadCreateNodeForce(loadCase, node03, loadDirection, 10);

    std::cout << "**********************************************" << std::endl;
    std::cout << "                  Visualisation               " << std::endl;
    std::cout << "**********************************************" << std::endl;

    myStructure.AddVisualizationComponentDisplacements();
    myStructure.AddVisualizationComponentSection();
    myStructure.AddVisualizationComponentConstitutive();

    std::cout << "**********************************************" << std::endl;
    std::cout << "                  Solver                      " << std::endl;
    std::cout << "**********************************************" << std::endl;

    myStructure.Info();
    //    myStructure.NodeBuildGlobalDofs();
    // myStructure.CheckCoefficientMatrix_0(1e-6, true);
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

    myStructure.Info();

    myStructure.ExportVtkDataFileElements(Parameters::mOutputPath.string() + "Truss2D.vtk");

    myStructure.IntegrationTypeInfo(10);
    std::cout << "**********************************************" << std::endl;
    std::cout << "                  End                         " << std::endl;
    std::cout << "**********************************************" << std::endl;
}
