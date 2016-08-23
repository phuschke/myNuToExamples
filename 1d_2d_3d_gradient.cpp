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
#include <vector>

const double youngsModulus = 40000;
const double poissonsRatio = 0.0;
const double nonlocalRadius = 10;
const double tensileStrength = 2;
const double compressiveStrength = 20;
const double fractureEnergy = 0.5;
const double displacement = 1;
const double lengthX = 200;
const double lengthY = 10;
const double lengthZ = 10;
const double simulationTime = 1;

int assignMaterial(NuTo::Structure& myStructure, const double factor)
{
    NuTo::FullVector<double, Eigen::Dynamic> myDamageLaw(1);
    myDamageLaw(0) = NuTo::Constitutive::eDamageLawType::ISOTROPIC_EXPONENTIAL_SOFTENING;

    int myMaterial = myStructure.ConstitutiveLawCreate(NuTo::Constitutive::eConstitutiveType::GRADIENT_DAMAGE_ENGINEERING_STRESS);

    myStructure.ConstitutiveLawSetPoissonsRatio(myMaterial, poissonsRatio);
    myStructure.ConstitutiveLawSetYoungsModulus(myMaterial, factor * youngsModulus);
    myStructure.ConstitutiveLawSetNonlocalRadius(myMaterial, nonlocalRadius);
    myStructure.ConstitutiveLawSetTensileStrength(myMaterial, factor * tensileStrength);
    myStructure.ConstitutiveLawSetFractureEnergy(myMaterial, fractureEnergy);
    myStructure.ConstitutiveLawSetDamageLaw(myMaterial, myDamageLaw);
    myStructure.ConstitutiveLawSetCompressiveStrength(myMaterial, factor * compressiveStrength);

    return myMaterial;
}

void setIntegrationScheme(NuTo::NewmarkDirect& myIntegrationScheme)
{
    myIntegrationScheme.SetTimeStep(1e-2 * simulationTime);
    myIntegrationScheme.SetMinTimeStep(1e-5 * simulationTime);
    myIntegrationScheme.SetMaxTimeStep(1e-1 * simulationTime);
    myIntegrationScheme.SetMinTimeStepPlot(1e-2);
    myIntegrationScheme.SetToleranceForce(1e-6);
    myIntegrationScheme.SetAutomaticTimeStepping(true);
    myIntegrationScheme.SetVerboseLevel(0);
    myIntegrationScheme.SetShowTime(false);
    myIntegrationScheme.SetPerformLineSearch(false);
    myIntegrationScheme.SetCheckCoefficientMatrix(false);
}

void addVisualization(NuTo::Structure& myStructure)
{
    myStructure.AddVisualizationComponentDisplacements();
    myStructure.AddVisualizationComponentEngineeringStrain();
    myStructure.AddVisualizationComponentEngineeringStress();
    myStructure.AddVisualizationComponentDamage();
    myStructure.AddVisualizationComponentNonlocalEqStrain();
    myStructure.AddVisualizationComponentConstitutive();
}

void run1D()
{
    boost::filesystem::path resultDirectory("/home/phuschke/1D2D3D_Gradient_results_1D/");
    boost::filesystem::remove_all(resultDirectory);
    boost::filesystem::create_directory(resultDirectory);

    const int dimension = 1;
    NuTo::FullVector<double, dimension> directionX(
    { 1 });

    //**********************************************
    //          Structure 1D
    //**********************************************
    NuTo::Structure myStructure(dimension);

    //**********************************************
    //          Material 1D
    //**********************************************
    int myMaterial = assignMaterial(myStructure, 1.);
    int weakMaterial = assignMaterial(myStructure, 0.9);

    //**********************************************
    //          Section 1D
    //**********************************************
    int mySection = myStructure.SectionCreate(NuTo::Section::TRUSS);
    myStructure.SectionSetArea(mySection, lengthZ * lengthY);

    //**********************************************
    //          Geometry 1D
    //**********************************************

    //create nodes
    int numElements = 20;
    int numNodesX = numElements + 1;
    double deltaX = lengthX / numElements;

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
    myStructure.InterpolationTypeAdd(myInterpolationType, NuTo::Node::DISPLACEMENTS, NuTo::Interpolation::EQUIDISTANT2);
    myStructure.InterpolationTypeAdd(myInterpolationType, NuTo::Node::NONLOCALEQSTRAIN, NuTo::Interpolation::EQUIDISTANT1);
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

    NuTo::NewmarkDirect myIntegrationScheme(&myStructure);
    setIntegrationScheme(myIntegrationScheme);
    int weak_Nodes = myStructure.GroupCreate(NuTo::Groups::Nodes);
    myStructure.GroupAddNodeCoordinateRange(weak_Nodes, 0, 0.5 * lengthX - 0.1, 0.5 * lengthX + 0.1);
    int weak_Elements = myStructure.GroupCreate(NuTo::Groups::Elements);
    myStructure.GroupAddElementsFromNodes(weak_Elements, weak_Nodes, false);
    myStructure.ElementGroupSetConstitutiveLaw(weak_Elements, weakMaterial);

    //**********************************************
    //          BC 1D
    //**********************************************
    int nodesLeft = myStructure.GroupCreate(NuTo::Groups::Nodes);
    myStructure.GroupAddNodeCoordinateRange(nodesLeft, 0, 0, 1e-6);
    myStructure.ConstraintLinearSetDisplacementNodeGroup(nodesLeft, directionX, 0);

    //**********************************************
    //          Loads 1D
    //**********************************************
    int nodesRight = myStructure.GroupCreate(NuTo::Groups::Nodes);
    myStructure.GroupAddNodeCoordinateRange(nodesRight, 0, lengthX - 1e-6, lengthX + 1e-6);
    int load = myStructure.ConstraintLinearSetDisplacementNodeGroup(nodesRight, directionX, 0);

    //**********************************************
    //          Visualisation 1D
    //**********************************************

    addVisualization(myStructure);

    //**********************************************
    //          Solve 1D
    //**********************************************

    myStructure.NodeBuildGlobalDofs();
    myStructure.CalculateMaximumIndependentSets();

    NuTo::FullMatrix<double, 2, 2> dispRHS;
    dispRHS << 0, 0, simulationTime, displacement;

    myIntegrationScheme.SetTimeDependentConstraint(load, dispRHS);

    myIntegrationScheme.SetResultDirectory(resultDirectory.string(), false);
    myStructure.SetVerboseLevel(10);
    myStructure.Info();

    myIntegrationScheme.Solve(simulationTime);
    myStructure.ExportVtkDataFileElements(resultDirectory.string() + "out.vtk");

    std::cout << "End 1D" << std::endl;
}

void run2D()
{
    boost::filesystem::path resultDirectory("/home/phuschke/1D2D3D_Gradient_results_2D/");
    boost::filesystem::remove_all(resultDirectory);
    boost::filesystem::create_directory(resultDirectory);

    const int dimension = 2;
    NuTo::FullVector<double, dimension> directionX(
    { 1, 0 });
    //directionY
    NuTo::FullVector<double, dimension> directionY(
    { 0, 1 });

    //**********************************************
    //          Structure 2D
    //**********************************************
    NuTo::Structure myStructure(dimension);

    //**********************************************
    //          Material 2D
    //**********************************************
    int myMaterial = assignMaterial(myStructure, 1.);
    int weakMaterial = assignMaterial(myStructure, 0.9);

    //**********************************************
    //          Section 2D
    //**********************************************
    int mySection = myStructure.SectionCreate(NuTo::Section::PLANE_STRESS);
    myStructure.SectionSetThickness(mySection, lengthZ);

    //**********************************************
    //          Geometry 2D
    //**********************************************
    boost::filesystem::path meshFile("/home/phuschke/develop/nuto/myNutoExamples/beam2D.msh");
    NuTo::FullMatrix<int, Eigen::Dynamic, Eigen::Dynamic> createdGroupIds = myStructure.ImportFromGmsh(meshFile.string(), NuTo::ElementData::CONSTITUTIVELAWIP, NuTo::IpData::eIpDataType::STATICDATA);
    int groupId = createdGroupIds.GetValue(0, 0);

    int myInterpolationType = myStructure.InterpolationTypeCreate(NuTo::Interpolation::eShapeType::QUAD2D);
    myStructure.InterpolationTypeAdd(myInterpolationType, NuTo::Node::COORDINATES, NuTo::Interpolation::EQUIDISTANT1);
    myStructure.InterpolationTypeAdd(myInterpolationType, NuTo::Node::DISPLACEMENTS, NuTo::Interpolation::EQUIDISTANT2);
    myStructure.InterpolationTypeAdd(myInterpolationType, NuTo::Node::NONLOCALEQSTRAIN, NuTo::Interpolation::EQUIDISTANT1);
    myStructure.ElementGroupSetInterpolationType(groupId, myInterpolationType);
    myStructure.ElementTotalConvertToInterpolationType(1.e-6, 10);
    myStructure.InterpolationTypeSetIntegrationType(myInterpolationType, NuTo::IntegrationType::IntegrationType2D4NGauss4Ip, NuTo::IpData::STATICDATA);

    myStructure.ElementTotalSetConstitutiveLaw(myMaterial);
    myStructure.ElementTotalSetSection(mySection);

    NuTo::NewmarkDirect myIntegrationScheme(&myStructure);
    setIntegrationScheme(myIntegrationScheme);
    int weak_Nodes = myStructure.GroupCreate(NuTo::Groups::Nodes);
    myStructure.GroupAddNodeCoordinateRange(weak_Nodes, 0, 0.5 * lengthX - 0.1, 0.5 * lengthX + 0.1);
    int weak_Elements = myStructure.GroupCreate(NuTo::Groups::Elements);
    myStructure.GroupAddElementsFromNodes(weak_Elements, weak_Nodes, false);
    myStructure.ElementGroupSetConstitutiveLaw(weak_Elements, weakMaterial);

    //**********************************************
    //          BC 2D
    //**********************************************
    int nodesLeft = myStructure.GroupCreate(NuTo::Groups::Nodes);
    myStructure.GroupAddNodeCoordinateRange(nodesLeft, 0, 0, 1e-6);
    myStructure.ConstraintLinearSetDisplacementNodeGroup(nodesLeft, directionX, 0);

    // fix single node
    NuTo::FullVector<double, 2> center;
    center[0] = 0;
    center[1] = 0;
    int node_0_0_0 = myStructure.GroupCreate(NuTo::Groups::Nodes);
    myStructure.GroupAddNodeRadiusRange(node_0_0_0, center, 0, 1e-6);
    myStructure.ConstraintLinearSetDisplacementNodeGroup(node_0_0_0, directionY, 0);

    // fix single node
    center[0] = lengthX;
    center[1] = 0;
    int node_lenX_0_0 = myStructure.GroupCreate(NuTo::Groups::Nodes);
    myStructure.GroupAddNodeRadiusRange(node_lenX_0_0, center, 0, 1e-6);
    myStructure.ConstraintLinearSetDisplacementNodeGroup(node_lenX_0_0, directionY, 0);

    //**********************************************
    //          Loads 2D
    //**********************************************
    int nodesRight = myStructure.GroupCreate(NuTo::Groups::Nodes);
    myStructure.GroupAddNodeCoordinateRange(nodesRight, 0, lengthX - 1e-6, lengthX + 1e-6);
    int load = myStructure.ConstraintLinearSetDisplacementNodeGroup(nodesRight, directionX, 0);

    //**********************************************
    //          Visualisation 2D
    //**********************************************

    addVisualization(myStructure);

    //**********************************************
    //          Solve 2D
    //**********************************************

    myStructure.NodeBuildGlobalDofs();
    myStructure.CalculateMaximumIndependentSets();

    NuTo::FullMatrix<double, 2, 2> dispRHS;
    dispRHS << 0, 0, simulationTime, displacement;

    myIntegrationScheme.SetTimeDependentConstraint(load, dispRHS);

    myIntegrationScheme.SetResultDirectory(resultDirectory.string(), false);
    myStructure.SetVerboseLevel(10);
    myStructure.Info();

    myIntegrationScheme.Solve(simulationTime);
    myStructure.ExportVtkDataFileElements(resultDirectory.string() + "out.vtk");

    std::cout << "End 2D" << std::endl;

}

void run3D()
{

    boost::filesystem::path resultDirectory("/home/phuschke/1D2D3D_Gradient_results_3D/");
    boost::filesystem::remove_all(resultDirectory);
    boost::filesystem::create_directory(resultDirectory);

    const int dimension = 3;
    NuTo::FullVector<double, dimension> directionX(
    { 1, 0, 0 });
    //directionY
    NuTo::FullVector<double, dimension> directionY(
    { 0, 1, 0 });
    //directionZ
    NuTo::FullVector<double, dimension> directionZ(
    { 0, 0, 1 });

    //**********************************************
    //          Structure 3D
    //**********************************************
    NuTo::Structure myStructure(dimension);

    //**********************************************
    //          Material 3D
    //**********************************************
    int myMaterial = assignMaterial(myStructure, 1.);
    int weakMaterial = assignMaterial(myStructure, 0.9);

    //**********************************************
    //          Section 3D
    //**********************************************
    int mySection = myStructure.SectionCreate(NuTo::Section::VOLUME);

    //**********************************************
    //          Geometry 3D
    //**********************************************
    boost::filesystem::path meshFile("/home/phuschke/develop/nuto/myNutoExamples/beam3D.msh");
    NuTo::FullMatrix<int, Eigen::Dynamic, Eigen::Dynamic> createdGroupIds = myStructure.ImportFromGmsh(meshFile.string(), NuTo::ElementData::CONSTITUTIVELAWIP, NuTo::IpData::eIpDataType::STATICDATA);
    int groupId = createdGroupIds.GetValue(0, 0);

    int myInterpolationType = myStructure.InterpolationTypeCreate(NuTo::Interpolation::eShapeType::BRICK3D);
    myStructure.InterpolationTypeAdd(myInterpolationType, NuTo::Node::COORDINATES, NuTo::Interpolation::EQUIDISTANT1);
    myStructure.InterpolationTypeAdd(myInterpolationType, NuTo::Node::DISPLACEMENTS, NuTo::Interpolation::EQUIDISTANT2);
    myStructure.InterpolationTypeAdd(myInterpolationType, NuTo::Node::NONLOCALEQSTRAIN, NuTo::Interpolation::EQUIDISTANT1);
    myStructure.ElementGroupSetInterpolationType(groupId, myInterpolationType);
    myStructure.ElementTotalConvertToInterpolationType(1.e-6, 10);
    myStructure.InterpolationTypeSetIntegrationType(myInterpolationType, NuTo::IntegrationType::IntegrationType3D8NGauss2x2x2Ip, NuTo::IpData::STATICDATA);

    myStructure.ElementTotalSetConstitutiveLaw(myMaterial);
    myStructure.ElementTotalSetSection(mySection);

    NuTo::NewmarkDirect myIntegrationScheme(&myStructure);
    setIntegrationScheme(myIntegrationScheme);
    int weak_Nodes = myStructure.GroupCreate(NuTo::Groups::Nodes);
    myStructure.GroupAddNodeCoordinateRange(weak_Nodes, 0, 0.5 * lengthX - 0.1, 0.5 * lengthX + 0.1);
    int weak_Elements = myStructure.GroupCreate(NuTo::Groups::Elements);
    myStructure.GroupAddElementsFromNodes(weak_Elements, weak_Nodes, false);
    myStructure.ElementGroupSetConstitutiveLaw(weak_Elements, weakMaterial);

    //**********************************************
    //          BC 3D
    //**********************************************
    int nodesLeft = myStructure.GroupCreate(NuTo::Groups::Nodes);
    myStructure.GroupAddNodeCoordinateRange(nodesLeft, 0, 0, 1e-6);
    myStructure.ConstraintLinearSetDisplacementNodeGroup(nodesLeft, directionX, 0);

    int nodesBottom = myStructure.GroupCreate(NuTo::Groups::Nodes);
    myStructure.GroupAddNodeCoordinateRange(nodesBottom, 1, 0, 1e-6);
    myStructure.ConstraintLinearSetDisplacementNodeGroup(nodesBottom, directionY, 0);

    int nodesFront = myStructure.GroupCreate(NuTo::Groups::Nodes);
    myStructure.GroupAddNodeCoordinateRange(nodesFront, 2, 0, 1e-6);
    myStructure.ConstraintLinearSetDisplacementNodeGroup(nodesFront, directionZ, 0);

    //**********************************************
    //          Loads 3D
    //**********************************************
    int nodesRight = myStructure.GroupCreate(NuTo::Groups::Nodes);
    myStructure.GroupAddNodeCoordinateRange(nodesRight, 0, lengthX - 1e-6, lengthX + 1e-6);
    int load = myStructure.ConstraintLinearSetDisplacementNodeGroup(nodesRight, directionX, 0);

    //**********************************************
    //          Visualisation 3D
    //**********************************************

    addVisualization(myStructure);

    //**********************************************
    //          Solve 3D
    //**********************************************

    myStructure.NodeBuildGlobalDofs();
    myStructure.CalculateMaximumIndependentSets();

    NuTo::FullMatrix<double, 2, 2> dispRHS;
    dispRHS << 0, 0, simulationTime, displacement;

    myIntegrationScheme.SetTimeDependentConstraint(load, dispRHS);

    myIntegrationScheme.SetResultDirectory(resultDirectory.string(), false);
    myStructure.SetVerboseLevel(10);
    myStructure.Info();

    myIntegrationScheme.Solve(simulationTime);
    myStructure.ExportVtkDataFileElements(resultDirectory.string() + "out.vtk");

    std::cout << "End 3D" << std::endl;
}

int main(int argc, char* argv[])
{

    run3D();
    run2D();
    run1D();

    std::cout << " ===> End all <=== " << std::endl;
}

