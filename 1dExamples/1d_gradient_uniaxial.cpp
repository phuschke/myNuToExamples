#include "nuto/mechanics/structures/unstructured/Structure.h"
#include "nuto/mechanics/timeIntegration/NewmarkDirect.h"
#include "nuto/mechanics/sections/SectionTruss.h"
#include <boost/filesystem.hpp>

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <chrono>

class Parameters
{
public:

    static const int mDimension = 1;

    static const bool mPerformLineSearch = true;
    static const bool mAutomaticTimeStepping = true;

    static constexpr double mMatrixYoungsModulus = 3.0e4;   // concrete
    static constexpr double mMatrixPoissonsRatio = 0.2;
    static constexpr double mMatrixThickness = 0.1;
    static constexpr double mMatrixNonlocalRadius = 2;
    static constexpr double mMatrixTensileStrength = 3;
    static constexpr double mMatrixCompressiveStrength = 30;
    static constexpr double mMatrixFractureEnergy = 0.1;


    static constexpr double mCrossSection = 0.1;
    static constexpr double mTimeStep = 1e-2;
    static constexpr double mMinTimeStep = 1e-7;
    static constexpr double mMaxTimeStep = 1e-1;
    static constexpr double mToleranceForce = 1e-4;
    static constexpr double mSimulationTime = 1.0;
    static constexpr double mLoad =0.2;

    static const boost::filesystem::path mOutputPath;

    static const NuTo::FullVector<double, 1> mDirectionX;

};

const boost::filesystem::path Parameters::mOutputPath("/home/phuschke/1d_gradient/");

const NuTo::FullVector<double, 1> Parameters::mDirectionX = NuTo::FullVector<double, 1>::UnitX();





int main(int argc, char* argv[])
{

    const int numElements = 100;
    const int dispOrder = 2;
    const int nlOrder = 1;
    const int ipOrder = 2;

    const double length = 100;

    std::cout << "***********************************" << std::endl;
    std::cout << "**      Structure                **" << std::endl;
    std::cout << "***********************************" << std::endl;

    NuTo::Structure myStructure(Parameters::mDimension);

    char resultDirectoryName[200];
    sprintf(resultDirectoryName, "/home/phuschke/1d_gradient_results_ele_%04d_dispOrder_%02d_nlOrder_%02d_ip_%02d/", numElements, dispOrder, nlOrder, ipOrder);
    boost::filesystem::path resultDirectory(resultDirectoryName);
    boost::filesystem::remove_all(resultDirectory);
    boost::filesystem::create_directory(resultDirectory);

    std::cout << "***********************************" << std::endl;
    std::cout << "**      Integration Scheme       **" << std::endl;
    std::cout << "***********************************" << std::endl;

    NuTo::NewmarkDirect myIntegrationScheme(&myStructure);
    myIntegrationScheme.SetTimeStep(1e-6 * Parameters::mSimulationTime);
    myIntegrationScheme.SetMinTimeStep(1e-7 * Parameters::mSimulationTime);
    myIntegrationScheme.SetMaxTimeStep(1e-3 * Parameters::mSimulationTime);
    myIntegrationScheme.SetToleranceForce(1e-4);
    myIntegrationScheme.SetAutomaticTimeStepping(true);
    myIntegrationScheme.SetVerboseLevel(0);
    myIntegrationScheme.SetShowTime(false);
    myIntegrationScheme.SetPerformLineSearch(false);
    myIntegrationScheme.SetCheckCoefficientMatrix(false);

    std::cout << "***********************************" << std::endl;
    std::cout << "**      Section                  **" << std::endl;
    std::cout << "***********************************" << std::endl;

    double alpha = 0.1;
    double xWeakSpot = 50.;
    double lWeakSpot = 10.;
    double areaParameters[4];
    areaParameters[0] = xWeakSpot;
    areaParameters[1] = lWeakSpot;
    areaParameters[2] = alpha;
    areaParameters[3] = 1;

    int mySection = myStructure.SectionCreate(NuTo::Section::TRUSS);
    myStructure.SectionSetArea(mySection, Parameters::mCrossSection);

    myStructure.SectionGetSectionPtr(mySection)->AsSectionTruss()->SetAreaParameters(areaParameters);

    std::cout << "***********************************" << std::endl;
    std::cout << "**      Material                 **" << std::endl;
    std::cout << "***********************************" << std::endl;

    int matrixMaterial = myStructure.ConstitutiveLawCreate(NuTo::Constitutive::eConstitutiveType::GRADIENT_DAMAGE_ENGINEERING_STRESS);
    myStructure.ConstitutiveLawSetParameterDouble(matrixMaterial, NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, Parameters::mMatrixYoungsModulus);
    myStructure.ConstitutiveLawSetParameterDouble(matrixMaterial, NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO, Parameters::mMatrixPoissonsRatio);
    myStructure.ConstitutiveLawSetParameterDouble(matrixMaterial, NuTo::Constitutive::eConstitutiveParameter::TENSILE_STRENGTH, Parameters::mMatrixTensileStrength);
    myStructure.ConstitutiveLawSetParameterDouble(matrixMaterial, NuTo::Constitutive::eConstitutiveParameter::COMPRESSIVE_STRENGTH, Parameters::mMatrixCompressiveStrength);
    myStructure.ConstitutiveLawSetParameterDouble(matrixMaterial, NuTo::Constitutive::eConstitutiveParameter::NONLOCAL_RADIUS, Parameters::mMatrixNonlocalRadius);
    myStructure.ConstitutiveLawSetParameterDouble(matrixMaterial, NuTo::Constitutive::eConstitutiveParameter::FRACTURE_ENERGY, Parameters::mMatrixFractureEnergy);

    //**********************************************
    //          Geometry
    //**********************************************

    //create nodes
    int numNodesX = numElements + 1;
    double deltaX = length / numElements;

    int nodeNum = 0;
    for (int countX = 0; countX < numNodesX; countX++)
    {
        NuTo::FullVector<double, Eigen::Dynamic> coordinates(1);
        coordinates(0) = countX * deltaX;
        myStructure.NodeCreate(nodeNum, coordinates);
        nodeNum++;
    }

    int myInterpolationType = myStructure.InterpolationTypeCreate(NuTo::Interpolation::eShapeType::TRUSS1D);
    myStructure.InterpolationTypeAdd(myInterpolationType, NuTo::Node::COORDINATES, NuTo::Interpolation::EQUIDISTANT1);

    switch (dispOrder)
    {
    case 2:
        myStructure.InterpolationTypeAdd(myInterpolationType, NuTo::Node::DISPLACEMENTS, NuTo::Interpolation::EQUIDISTANT2);
        break;
    case 3:
        myStructure.InterpolationTypeAdd(myInterpolationType, NuTo::Node::DISPLACEMENTS, NuTo::Interpolation::EQUIDISTANT3);
        break;
    case 4:
        myStructure.InterpolationTypeAdd(myInterpolationType, NuTo::Node::DISPLACEMENTS, NuTo::Interpolation::EQUIDISTANT4);
        break;
    default:
        return 1;
    }

    switch (nlOrder)
    {
    case 1:
        myStructure.InterpolationTypeAdd(myInterpolationType, NuTo::Node::NONLOCALEQSTRAIN, NuTo::Interpolation::EQUIDISTANT1);
        break;
    case 2:
        myStructure.InterpolationTypeAdd(myInterpolationType, NuTo::Node::NONLOCALEQSTRAIN, NuTo::Interpolation::EQUIDISTANT2);
        break;
    case 3:
        myStructure.InterpolationTypeAdd(myInterpolationType, NuTo::Node::NONLOCALEQSTRAIN, NuTo::Interpolation::EQUIDISTANT3);
        break;
    default:
        return 1;
    }

    switch (ipOrder)
    {
    case 2:
        myStructure.InterpolationTypeSetIntegrationType(myInterpolationType, NuTo::IntegrationType::IntegrationType1D2NGauss2Ip, NuTo::IpData::STATICDATA);
        break;
    case 3:
        myStructure.InterpolationTypeSetIntegrationType(myInterpolationType, NuTo::IntegrationType::IntegrationType1D2NGauss3Ip, NuTo::IpData::STATICDATA);
        break;
    case 4:
        myStructure.InterpolationTypeSetIntegrationType(myInterpolationType, NuTo::IntegrationType::IntegrationType1D2NGauss4Ip, NuTo::IpData::STATICDATA);
        break;
    case 5:
        myStructure.InterpolationTypeSetIntegrationType(myInterpolationType, NuTo::IntegrationType::IntegrationType1D2NGauss5Ip, NuTo::IpData::STATICDATA);
        break;
    default:
        return 1;
    }

    //create elements
    for (int countX = 0; countX < numElements; countX++)
    {
        NuTo::FullVector<int, Eigen::Dynamic> nodes(2);
        nodes(0) = countX;
        nodes(1) = countX + 1;
        myStructure.ElementCreate(myInterpolationType, nodes, NuTo::ElementData::eElementDataType::CONSTITUTIVELAWIP, NuTo::IpData::eIpDataType::STATICDATA);
    }

    myStructure.ElementTotalConvertToInterpolationType(1.e-6, 10);

    myStructure.ElementTotalSetConstitutiveLaw(matrixMaterial);
    myStructure.ElementTotalSetSection(mySection);

    int weakNodes = myStructure.GroupCreate(NuTo::Groups::Nodes);
    myStructure.GroupAddNodeCoordinateRange(weakNodes, 0, 45 - 1e-6, 55 + 1e-6);
    int weakElemets = myStructure.GroupCreate(NuTo::Groups::Elements);
    myStructure.GroupAddElementsFromNodes(weakElemets, weakNodes, true);

    myStructure.IntegrationTypeInfo(10);

    //**********************************************
    //          Boundary Conditions
    //**********************************************

    // left boundary
    int nodesLeft = myStructure.GroupCreate(NuTo::Groups::Nodes);
    myStructure.GroupAddNodeCoordinateRange(nodesLeft, 0, -1e-6, 1e-6);
    myStructure.ConstraintLinearSetDisplacementNodeGroup(nodesLeft, Parameters::mDirectionX, 0);

    //**********************************************
    //          Loads
    //**********************************************

    // right load
    int nodesRight = myStructure.GroupCreate(NuTo::Groups::Nodes);
    myStructure.GroupAddNodeCoordinateRange(nodesRight, 0, length - 1e-6, length + 1e-6);
    int load = myStructure.ConstraintLinearSetDisplacementNodeGroup(nodesRight, Parameters::mDirectionX, 0);

    //**********************************************
    //          Visualisation
    //**********************************************

    myStructure.AddVisualizationComponentDisplacements();
    myStructure.AddVisualizationComponentDisplacements(weakElemets);

//    myStructure.AddVisualizationComponentEngineeringStrain();
//    myStructure.AddVisualizationComponentEngineeringStress();
//    myStructure.AddVisualizationComponentDamage();
//    myStructure.AddVisualizationComponentNonlocalEqStrain();
//    myStructure.AddVisualizationComponentSection();
//    myStructure.AddVisualizationComponentConstitutive();

    //**********************************************
    //          Solver
    //**********************************************

    myStructure.NodeBuildGlobalDofs();
    myStructure.CalculateMaximumIndependentSets();

    NuTo::FullMatrix<double, 2, 2> timeDependentLoad;
    timeDependentLoad(0, 0) = 0;
    timeDependentLoad(1, 0) = Parameters::mSimulationTime;
    timeDependentLoad(0, 1) = 0;
    timeDependentLoad(1, 1) = Parameters::mLoad;

    myIntegrationScheme.AddResultGroupNodeForce("force", nodesRight);
    myIntegrationScheme.AddResultNodeDisplacements("displacements", myStructure.GroupGetMemberIds(nodesRight).GetValue(0, 0));

    myIntegrationScheme.SetTimeDependentConstraint(load, timeDependentLoad);

    myIntegrationScheme.SetResultDirectory(resultDirectory.string(), false);
    myStructure.SetVerboseLevel(10);
    myStructure.Info();

    myIntegrationScheme.Solve(Parameters::mSimulationTime);

    std::cout << " ===> End <=== " << std::endl;
}

