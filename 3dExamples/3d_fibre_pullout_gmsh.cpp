#include "nuto/math/FullMatrix.h"
#include "nuto/mechanics/structures/unstructured/Structure.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStrain3D.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStress3D.h"
#include "nuto/mechanics/constitutive/mechanics/GradientDamageEngineeringStress.h"
#include "nuto/mechanics/MechanicsException.h"
#include "nuto/mechanics/constitutive/ConstitutiveStaticDataBase.h"
#include "nuto/mechanics/constitutive/mechanics/ConstitutiveStaticDataGradientDamage1D.h"

#include "nuto/math/SparseMatrixCSRVector2General.h"

#include "nuto/mechanics/timeIntegration/NewmarkDirect.h"
#include "nuto/base/Logger.h"

#include <boost-1_55/boost/filesystem.hpp>
#include <boost-1_55/boost/lexical_cast.hpp>
#include <boost-1_55/boost/progress.hpp>

#include <iostream>
#include <fstream>
#include <string>

class Parameters
{
public:

    static const int mDimension = 3;

    static const bool mPerformLineSearch = true;
    static const bool mAutomaticTimeStepping = true;

    static constexpr double mMatrixYoungsModulus = 4.0e4;   // concrete
    static constexpr double mMatrixPoissonsRatio = 0.2;

    static constexpr double mFibreYoungsModulus = 2.1e9;
    static constexpr double mFibrePoissonsRatio = 0.2;
    static constexpr double mFibreCrossSection = 0.01;

    static constexpr double mNonlocalRadius = 4;
    static constexpr double mTensileStrength = 3;
    static constexpr double mCompressiveStrength = 30;
    static constexpr double mFractureEnergy = 0.01;
    static constexpr double mTimeStep = 1e-3;
    static constexpr double mMinTimeStep = 1e-5;
    static constexpr double mMaxTimeStep = 1e-1;
    static constexpr double mToleranceForce = 1e-6;
    static constexpr double mSimulationTime = 1.0;
    static constexpr double mLoad = 1.0;

    static const boost::filesystem::path mOutputPath;
    static const boost::filesystem::path mMeshFilePath;

    static const NuTo::FullVector<double, 3> mDirectionX;
    static const NuTo::FullVector<double, 3> mDirectionY;
    static const NuTo::FullVector<double, 3> mDirectionZ;
};

const boost::filesystem::path Parameters::mOutputPath("/home/phuschke/3d_fibre_pullout_result/");
const boost::filesystem::path Parameters::mMeshFilePath("/home/phuschke/develop/nuto/myNutoExamples/3d_fibre_pullout_mesh.msh");

const NuTo::FullVector<double, 3> Parameters::mDirectionX = NuTo::FullVector<double, 3>::UnitX();
const NuTo::FullVector<double, 3> Parameters::mDirectionY = NuTo::FullVector<double, 3>::UnitY();
const NuTo::FullVector<double, 3> Parameters::mDirectionZ = NuTo::FullVector<double, 3>::UnitZ();

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

    int matrixSection = myStructure.SectionCreate(NuTo::Section::VOLUME);

    int fibreSection = myStructure.SectionCreate(NuTo::Section::TRUSS);
    myStructure.SectionSetArea(fibreSection, Parameters::mFibreCrossSection);

    std::cout << "**********************************************" << std::endl;
    std::cout << "                  Material                    " << std::endl;
    std::cout << "**********************************************" << std::endl;

    int matrixMaterial = myStructure.ConstitutiveLawCreate(NuTo::Constitutive::eConstitutiveType::LINEAR_ELASTIC_ENGINEERING_STRESS);
    myStructure.ConstitutiveLawSetYoungsModulus(matrixMaterial, Parameters::mMatrixYoungsModulus);
    myStructure.ConstitutiveLawSetPoissonsRatio(matrixMaterial, Parameters::mMatrixPoissonsRatio);

    int fibreMaterial = myStructure.ConstitutiveLawCreate(NuTo::Constitutive::eConstitutiveType::LINEAR_ELASTIC_ENGINEERING_STRESS);
    myStructure.ConstitutiveLawSetYoungsModulus(fibreMaterial, Parameters::mFibreYoungsModulus);
    myStructure.ConstitutiveLawSetPoissonsRatio(fibreMaterial, Parameters::mFibrePoissonsRatio);

    std::cout << "**********************************************" << std::endl;
    std::cout << "                  Interpolation Type          " << std::endl;
    std::cout << "**********************************************" << std::endl;

    int matrixInterpolationType = myStructure.InterpolationTypeCreate(NuTo::Interpolation::eShapeType::BRICK3D);
    myStructure.InterpolationTypeAdd(matrixInterpolationType, NuTo::Node::COORDINATES, NuTo::Interpolation::EQUIDISTANT1);
    myStructure.InterpolationTypeAdd(matrixInterpolationType, NuTo::Node::DISPLACEMENTS, NuTo::Interpolation::EQUIDISTANT1);
    myStructure.InterpolationTypeSetIntegrationType(matrixInterpolationType, NuTo::IntegrationType::IntegrationType3D8NGauss2x2x2Ip, NuTo::IpData::STATICDATA);

    int fibreInterpolationType = myStructure.InterpolationTypeCreate(NuTo::Interpolation::eShapeType::TRUSSXD);
    myStructure.InterpolationTypeAdd(fibreInterpolationType, NuTo::Node::COORDINATES, NuTo::Interpolation::EQUIDISTANT1);
    myStructure.InterpolationTypeAdd(fibreInterpolationType, NuTo::Node::DISPLACEMENTS, NuTo::Interpolation::EQUIDISTANT1);
    myStructure.InterpolationTypeSetIntegrationType(fibreInterpolationType, NuTo::IntegrationType::IntegrationType1D2NGauss1Ip, NuTo::IpData::STATICDATA);

    std::cout << "**********************************************" << std::endl;
    std::cout << "                  Matrix Geometry             " << std::endl;
    std::cout << "**********************************************" << std::endl;

    NuTo::FullMatrix<int, Eigen::Dynamic, Eigen::Dynamic> createdGroupIds = myStructure.ImportFromGmsh(Parameters::mMeshFilePath.string(), NuTo::ElementData::CONSTITUTIVELAWIP, NuTo::IpData::eIpDataType::STATICDATA);
    int groupIdMatrix = createdGroupIds.GetValue(0, 0);

    myStructure.ElementGroupSetInterpolationType(groupIdMatrix, matrixInterpolationType);

    myStructure.ElementGroupSetSection(groupIdMatrix, matrixSection);
    myStructure.ElementGroupSetConstitutiveLaw(groupIdMatrix, matrixMaterial);

    std::cout << "**********************************************" << std::endl;
    std::cout << "                  Fibre Geometry              " << std::endl;
    std::cout << "**********************************************" << std::endl;

    myStructure.ElementTotalConvertToInterpolationType(1.e-6, 10);
    myStructure.NodeBuildGlobalDofs();
    myStructure.Info();

    NuTo::FullVector<int, Eigen::Dynamic> nodeIndicesFibre(2);
    NuTo::FullVector<double, Eigen::Dynamic> nodeCoords(Parameters::mDimension);

    nodeCoords[0] = 0.0;
    nodeCoords[1] = 0.0;
    nodeCoords[2] = 0.0;
    int node00 = myStructure.NodeGetIdAtCoordinate(nodeCoords, 1e-6);

    nodeCoords[0] = 0.0;
    nodeCoords[1] = 0.0;
    nodeCoords[2] = 1.0;
    int node01 = myStructure.NodeGetIdAtCoordinate(nodeCoords, 1e-6);

    nodeIndicesFibre[0] = node00;
    nodeIndicesFibre[1] = node01;
    int elementFibre01 = myStructure.ElementCreate(fibreInterpolationType, nodeIndicesFibre, NuTo::ElementData::eElementDataType::CONSTITUTIVELAWIP, NuTo::IpData::eIpDataType::NOIPDATA);

    int groupElementFibre = myStructure.GroupCreate(NuTo::Groups::eGroupId::Elements);
//    myStructure.GroupAddElement(groupElementFibre, elementFibre00);
    myStructure.GroupAddElement(groupElementFibre, elementFibre01);
    myStructure.ElementGroupSetSection(groupElementFibre, fibreSection);
    myStructure.ElementGroupSetConstitutiveLaw(groupElementFibre, fibreMaterial);

    myStructure.ElementTotalConvertToInterpolationType(1.e-6, 10);

    myStructure.NodeBuildGlobalDofs();
    myStructure.Info();
    myStructure.CheckCoefficientMatrix_0(1e-6, true);

    std::cout << "**********************************************" << std::endl;
    std::cout << "                  Boundary Conditions         " << std::endl;
    std::cout << "**********************************************" << std::endl;

    int groupNodeBCLeft = myStructure.GroupCreate(NuTo::Groups::eGroupId::Nodes);
    myStructure.GroupAddNodeCoordinateRange(groupNodeBCLeft, 0, -1e-6, +1e-6);

    int groupNodeBCFront = myStructure.GroupCreate(NuTo::Groups::eGroupId::Nodes);
    myStructure.GroupAddNodeCoordinateRange(groupNodeBCFront, 1, -1e-6, +1e-6);


    int groupNodeBCBottom = myStructure.GroupCreate(NuTo::Groups::eGroupId::Nodes);
    myStructure.GroupAddNodeCoordinateRange(groupNodeBCBottom, 2, -1e-6, +1e-6);


    myStructure.ConstraintLinearSetDisplacementNodeGroup(groupNodeBCLeft, Parameters::mDirectionX, 0);
    myStructure.ConstraintLinearSetDisplacementNodeGroup(groupNodeBCFront, Parameters::mDirectionY, 0);
    myStructure.ConstraintLinearSetDisplacementNodeGroup(groupNodeBCBottom, Parameters::mDirectionZ, 0);
    std::cout << "**********************************************" << std::endl;
    std::cout << "                  Loads                       " << std::endl;
    std::cout << "**********************************************" << std::endl;

    int loadCase = 0;
    int groupNodeLoad = myStructure.GroupCreate(NuTo::Groups::eGroupId::Nodes);
    myStructure.GroupAddNodeCoordinateRange(groupNodeLoad, 0, 1 - 1e-6, 1 + 1e-6);
    //myStructure.LoadCreateNodeGroupForce(loadCase, groupNodeLoad, Parameters::mDirectionX, 1);
    int timeDependentConstraint = myStructure.ConstraintLinearSetDisplacementNodeGroup(groupNodeLoad, Parameters::mDirectionX, 1);

    std::cout << "**********************************************" << std::endl;
    std::cout << "                  Visualisation               " << std::endl;
    std::cout << "**********************************************" << std::endl;

    myStructure.AddVisualizationComponentDisplacements();
    myStructure.AddVisualizationComponentEngineeringStrain();
    myStructure.AddVisualizationComponentEngineeringStress();
    myStructure.AddVisualizationComponentSection();
    myStructure.AddVisualizationComponentConstitutive();

    std::cout << "**********************************************" << std::endl;
    std::cout << "                  Solver                      " << std::endl;
    std::cout << "**********************************************" << std::endl;

    myStructure.NodeBuildGlobalDofs();
    myStructure.CheckCoefficientMatrix_0(1e-6, true);
    myStructure.Info();

    boost::filesystem::remove_all(Parameters::mOutputPath);
    boost::filesystem::create_directory(Parameters::mOutputPath);

    myStructure.NodeBuildGlobalDofs();
    myStructure.CalculateMaximumIndependentSets();

    NuTo::FullMatrix<double, 2, 2> timeDependentLoad;
    timeDependentLoad(0, 0) = 0;
    timeDependentLoad(1, 0) = Parameters::mSimulationTime;
    timeDependentLoad(0, 1) = 0;
    timeDependentLoad(1, 1) = Parameters::mLoad;

    //myIntegrationScheme.SetTimeDependentLoadCase(loadCase, timeDependentLoad);
    myIntegrationScheme.SetTimeDependentConstraint(timeDependentConstraint, timeDependentLoad);

    myIntegrationScheme.SetResultDirectory(Parameters::mOutputPath.string(), false);

    myIntegrationScheme.Solve(Parameters::mSimulationTime);

    std::cout << "**********************************************" << std::endl;
    std::cout << "                  End                         " << std::endl;
    std::cout << "**********************************************" << std::endl;
}

