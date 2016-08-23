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

int main(int argc, char* argv[])
{

//    if (argc != 4)
//    {
//        std::cout << "input arguments: numElements, interpolation order, integration order" << std::endl;
//        return 1;
//    }
    const int numElements = 40;//std::stoi(argv[1]);
    const int order = 4;//std::stoi(argv[2]);
    const int ipOrder = 4;//std::stoi(argv[3]);

    const int dimension = 1;
    const double youngsModulus = 4500;
    const double poissonsRatio = 0.0;
    const double crossSection = 55;
    const double displacement = 10;

    const double length = 100;
    //directionX
    NuTo::FullVector<double, dimension> directionX;
    directionX[0] = 1;

    //**********************************************
    //          Structure
    //**********************************************

    NuTo::Structure myStructure(dimension);

    char resultDirectoryName[200];
    sprintf(resultDirectoryName, "/home/phuschke/1D_Displacement_results_ele_%04d_disp_order_%02d_ip_%02d/", numElements, order, ipOrder);
    boost::filesystem::path resultDirectory(resultDirectoryName);
    boost::filesystem::remove_all(resultDirectory);
    boost::filesystem::create_directory(resultDirectory);

    //**********************************************
    //         Integration Scheme
    //**********************************************

    double simulationTime = 1;
    NuTo::NewmarkDirect myIntegrationScheme(&myStructure);
    myIntegrationScheme.SetTimeStep(1e-4 * simulationTime);
    myIntegrationScheme.SetMinTimeStep(1e-5 * simulationTime);
    myIntegrationScheme.SetMaxTimeStep(1e-2 * simulationTime);
    myIntegrationScheme.SetToleranceForce(1e-6);
    myIntegrationScheme.SetAutomaticTimeStepping(true);
    myIntegrationScheme.SetVerboseLevel(0);
    myIntegrationScheme.SetShowTime(false);
    myIntegrationScheme.SetPerformLineSearch(false);
    myIntegrationScheme.SetCheckCoefficientMatrix(false);

    //**********************************************
    //          Section
    //**********************************************

    int mySection = myStructure.SectionCreate(NuTo::Section::TRUSS);
    myStructure.SectionSetArea(mySection, crossSection);

    int weakSection = myStructure.SectionCreate(NuTo::Section::TRUSS);
    myStructure.SectionSetArea(weakSection, 50);

    //**********************************************
    //          Material
    //**********************************************

    int myMaterial = myStructure.ConstitutiveLawCreate(NuTo::Constitutive::eConstitutiveType::LINEAR_ELASTIC_ENGINEERING_STRESS);
    myStructure.ConstitutiveLawSetPoissonsRatio(myMaterial, poissonsRatio);
    myStructure.ConstitutiveLawSetYoungsModulus(myMaterial, youngsModulus);


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

    // triangular elements
    if (order == 2)
    {
        myStructure.InterpolationTypeAdd(myInterpolationType, NuTo::Node::DISPLACEMENTS, NuTo::Interpolation::EQUIDISTANT2);
    } else if (order == 3)
    {
        myStructure.InterpolationTypeAdd(myInterpolationType, NuTo::Node::DISPLACEMENTS, NuTo::Interpolation::EQUIDISTANT3);
    } else if (order == 4)
    {
        myStructure.InterpolationTypeAdd(myInterpolationType, NuTo::Node::DISPLACEMENTS, NuTo::Interpolation::EQUIDISTANT4);
    } else
    {
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

    myStructure.ElementTotalSetConstitutiveLaw(myMaterial);
    myStructure.ElementTotalSetSection(mySection);

    int weakNodes = myStructure.GroupCreate(NuTo::Groups::Nodes);
    myStructure.GroupAddNodeCoordinateRange(weakNodes, 0, 45 - 1e-6, 55 + 1e-6);
    int weakElemets = myStructure.GroupCreate(NuTo::Groups::Elements);
    myStructure.GroupAddElementsFromNodes(weakElemets, weakNodes, true);
    myStructure.ElementGroupSetSection(weakElemets, weakSection);
    //myStructure.ElementGroupSetConstitutiveLaw(weakElemets, weakMaterial);
    myStructure.IntegrationTypeInfo(10);

    //**********************************************
    //          Boundary Conditions
    //**********************************************

    // left boundary
    int nodesLeft = myStructure.GroupCreate(NuTo::Groups::Nodes);
    myStructure.GroupAddNodeCoordinateRange(nodesLeft, 0, -1e-6, 1e-6);
    myStructure.ConstraintLinearSetDisplacementNodeGroup(nodesLeft, directionX, 0);

    //**********************************************
    //          Loads
    //**********************************************

    // right load
    int nodesRight = myStructure.GroupCreate(NuTo::Groups::Nodes);
    myStructure.GroupAddNodeCoordinateRange(nodesRight, 0, length - 1e-6, length + 1e-6);
    int load = myStructure.ConstraintLinearSetDisplacementNodeGroup(nodesRight, directionX, 0);

    //**********************************************
    //          Visualisation
    //**********************************************

    myStructure.AddVisualizationComponentDisplacements();
    myStructure.AddVisualizationComponentEngineeringStrain();
    myStructure.AddVisualizationComponentEngineeringStress();
    myStructure.AddVisualizationComponentSection();
    myStructure.AddVisualizationComponentConstitutive();
    //**********************************************
    //          Solver
    //**********************************************

    myStructure.NodeBuildGlobalDofs();
    myStructure.CalculateMaximumIndependentSets();

    NuTo::FullMatrix<double, 2, 2> dispRHS;
    dispRHS << 0, 0, simulationTime, displacement;

    myIntegrationScheme.AddResultGroupNodeForce("myforce", nodesRight);
    myIntegrationScheme.AddResultNodeDisplacements("mydisplacements", myStructure.GroupGetMemberIds(nodesRight).GetValue(0, 0));

    myIntegrationScheme.SetTimeDependentConstraint(load, dispRHS);

    myIntegrationScheme.SetResultDirectory(resultDirectory.string(), false);
    myStructure.SetVerboseLevel(10);
    myStructure.Info();

    myIntegrationScheme.Solve(simulationTime);
    myStructure.ExportVtkDataFileElements(resultDirectory.string() + "out.vtk");

    // output stress at ip
    boost::filesystem::path stressFile(resultDirectory.string() + "stressFile.dat");
    std::ofstream stressOutput;
    stressOutput.open(stressFile.string());

    NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> stress;
    NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> coords;

    myStructure.Info();
    for (int ele = 0; ele < myStructure.GetNumElements(); ++ele)
    {
        myStructure.ElementGetEngineeringStress(ele, stress);
        myStructure.ElementGetIntegrationPointCoordinates(ele, coords);

        int numIP = myStructure.ElementGetElementPtr(ele)->GetNumIntegrationPoints();
        for (int ip = 0; ip < numIP; ++ip)
        {
            stressOutput << stress.GetValue(0, ip) << "\t" << coords.GetValue(0, ip) << std::endl;
        }

    }
    stressOutput.close();



    char forceDispName[200];
    sprintf(forceDispName, "/home/phuschke/develop/nuto/myNutoExamples/GradientEnhancedDamageProject/latex/data/forceDisp1Ddisplacement_%03d_%02d_%02d.dat", numElements, order, ipOrder);
    std::string command = "paste " + resultDirectory.string() + "myforce.dat " + resultDirectory.string() + "mydisplacements.dat > " + forceDispName;
    system(command.c_str());

    char stressName[200];
    sprintf(stressName, "/home/phuschke/develop/nuto/myNutoExamples/GradientEnhancedDamageProject/latex/data/stress1Ddisplacement_%03d_%02d_%02d.dat", numElements, order, ipOrder);
    command = "mv " + resultDirectory.string() + "stressFile.dat " + stressName;
    system(command.c_str());
    std::cout << " ===> End 1DGradient <=== " << std::endl;
}

