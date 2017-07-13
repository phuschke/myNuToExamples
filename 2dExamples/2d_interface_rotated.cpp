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

constexpr unsigned int dimension = 2;

class Parameters
{
public:
    static const int mDimension = dimension;

    static const bool mPerformLineSearch = true;
    static const bool mAutomaticTimeStepping = false;

    static constexpr double mInterfaceNormalStiffness = 100;
    static constexpr double mInterfaceTangentialStiffness = 1;

    static constexpr double mInterfaceThickness = 10;


    static constexpr double mNonlocalRadius = 4;
    static constexpr double mTensileStrength = 3;
    static constexpr double mCompressiveStrength = 30;
    static constexpr double mFractureEnergy = 0.01;
    static constexpr double mTimeStep = 0.1;
    static constexpr double mMinTimeStep = 1e-4;
    static constexpr double mMaxTimeStep = 1.;
    static constexpr double mToleranceForce = 1e-6;
    static constexpr double mSimulationTime = 1.0;
    static constexpr double mLoad = 1.;

    static const boost::filesystem::path mOutputPath;
    static const boost::filesystem::path mMeshFilePath;

    static const NuTo::FullVector<double, dimension> mDirectionX;
    static const NuTo::FullVector<double, dimension> mDirectionY;
};

const boost::filesystem::path Parameters::mOutputPath("/home/phuschke/2d_interface_goodman_rotated/");

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
    myIntegrationScheme.SetResultDirectory(Parameters::mOutputPath.string(), true);

    std::cout << "**********************************************" << std::endl;
    std::cout << "                  Section                     " << std::endl;
    std::cout << "**********************************************" << std::endl;

    int interfaceSection = myStructure.SectionCreate(NuTo::Section::PLANE_STRESS);
    myStructure.SectionSetThickness(interfaceSection, Parameters::mInterfaceThickness);

    std::cout << "**********************************************" << std::endl;
    std::cout << "                  Material                    " << std::endl;
    std::cout << "**********************************************" << std::endl;

    int interfaceMaterial = myStructure.ConstitutiveLawCreate(NuTo::Constitutive::eConstitutiveType::INTERFACE_GOODMAN);
    myStructure.ConstitutiveLawSetParameterDouble(interfaceMaterial,
                                                  NuTo::Constitutive::eConstitutiveParameter::NORMAL_STIFFNESS,
                                                  Parameters::mInterfaceNormalStiffness);
    myStructure.ConstitutiveLawSetParameterDouble(interfaceMaterial,
                                                  NuTo::Constitutive::eConstitutiveParameter::TANGENTIAL_STIFFNESS,
                                                  Parameters::mInterfaceTangentialStiffness);

    std::cout << "**********************************************" << std::endl;
    std::cout << "                  Interpolation Type          " << std::endl;
    std::cout << "**********************************************" << std::endl;

    int matrixInterpolationType = myStructure.InterpolationTypeCreate(NuTo::Interpolation::eShapeType::INTERFACE);
    myStructure.InterpolationTypeAdd(matrixInterpolationType, NuTo::Node::COORDINATES,
                                     NuTo::Interpolation::EQUIDISTANT2);
    myStructure.InterpolationTypeAdd(matrixInterpolationType, NuTo::Node::DISPLACEMENTS,
                                     NuTo::Interpolation::EQUIDISTANT2);

    std::cout << "**********************************************" << std::endl;
    std::cout << "                  Nodes                       " << std::endl;
    std::cout << "**********************************************" << std::endl;

    NuTo::FullVector<double, Eigen::Dynamic> nodeCoords(Parameters::mDimension);
    std::set<NuTo::Node::eAttributes> dofs;
    dofs.insert(NuTo::Node::COORDINATES);
    dofs.insert(NuTo::Node::DISPLACEMENTS);

    nodeCoords[0] = 0.0;
    nodeCoords[1] = 0.0;
    int nodeInterface00 = myStructure.NodeCreate(nodeCoords, dofs);

    nodeCoords[0] = 1.0;
    nodeCoords[1] = 1.0;
    int nodeInterface01 = myStructure.NodeCreate(nodeCoords, dofs);

    nodeCoords[0] = 2.0;
    nodeCoords[1] = 2.0;
    int nodeInterface02 = myStructure.NodeCreate(nodeCoords, dofs);

    nodeCoords[0] = 0.0;
    nodeCoords[1] = 4.0;
    int nodeInterface03 = myStructure.NodeCreate(nodeCoords, dofs);

    nodeCoords[0] = -1.0;
    nodeCoords[1] = 3.0;
    int nodeInterface04 = myStructure.NodeCreate(nodeCoords, dofs);

    nodeCoords[0] = -2.0;
    nodeCoords[1] = 2.0;
    int nodeInterface05 = myStructure.NodeCreate(nodeCoords, dofs);

    std::cout << "**********************************************" << std::endl;
    std::cout << "                  Interface Geometry          " << std::endl;
    std::cout << "**********************************************" << std::endl;

    NuTo::FullVector<int, Eigen::Dynamic> nodeIndicesMatrix(6);

    nodeIndicesMatrix[0] = nodeInterface00;
    nodeIndicesMatrix[1] = nodeInterface01;
    nodeIndicesMatrix[2] = nodeInterface02;
    nodeIndicesMatrix[3] = nodeInterface03;
    nodeIndicesMatrix[4] = nodeInterface04;
    nodeIndicesMatrix[5] = nodeInterface05;
    int elementMatirx00 = myStructure.ElementCreate(matrixInterpolationType, nodeIndicesMatrix,
                                                    NuTo::ElementData::eElementDataType::CONSTITUTIVELAWIP,
                                                    NuTo::IpData::eIpDataType::NOIPDATA);

    myStructure.ElementSetSection(elementMatirx00, interfaceSection);
    myStructure.ElementSetConstitutiveLaw(elementMatirx00, interfaceMaterial);

    std::cout << "**********************************************" << std::endl;
    std::cout << "                  Boundary Conditions         " << std::endl;
    std::cout << "**********************************************" << std::endl;

    NuTo::FullVector<double, dimension> directionAligned = NuTo::FullVector<double, dimension>::Ones();
    NuTo::FullVector<double, dimension> directionOrthogonal;
    directionOrthogonal(0, 0) = -1.0;
    directionOrthogonal(1, 0) = 1.0;


    int groupNodeBC = myStructure.GroupCreate(NuTo::Groups::eGroupId::Nodes);
    myStructure.GroupAddNode(groupNodeBC, nodeInterface00);
    myStructure.GroupAddNode(groupNodeBC, nodeInterface01);
    myStructure.GroupAddNode(groupNodeBC, nodeInterface02);

    myStructure.ConstraintLinearSetDisplacementNodeGroup(groupNodeBC, directionAligned, 0);
    myStructure.ConstraintLinearSetDisplacementNodeGroup(groupNodeBC, directionOrthogonal, 0);

    std::cout << "**********************************************" << std::endl;
    std::cout << "                  Loads                       " << std::endl;
    std::cout << "**********************************************" << std::endl;
    int loadCase = 0;
    myStructure.LoadCreateNodeForce(loadCase, nodeInterface03, directionAligned, 1e4);

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

    myStructure.Info();
    myStructure.NodeBuildGlobalDofs();
    myStructure.CalculateMaximumIndependentSets();
    myStructure.CheckCoefficientMatrix_0(1e-6, true);
    myStructure.Info();


    NuTo::FullMatrix<double, -1, -1> dispRHS(2, 2);
    dispRHS << 0, 0, 1, 1;


    myIntegrationScheme.SetTimeDependentLoadCase(loadCase, dispRHS);
    myIntegrationScheme.Solve(Parameters::mSimulationTime);

    myStructure.CheckCoefficientMatrix_0(1e-6, true);
    std::cout << "**********************************************" << std::endl;
    std::cout << "                  End                         " << std::endl;
    std::cout << "**********************************************" << std::endl;
}
