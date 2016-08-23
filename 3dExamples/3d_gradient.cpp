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

int main(int argc, char* argv[])
{

    const int dimension               = 3;
    const double youngsModulus        = 40000;
    const double poissonsRatio        = 0.0;
    const double nonlocalRadius       = 1;
    const double tensileStrength      = 2;
    const double compressiveStrength  = 20;
    const double fractureEnergy       = 1;
    const double displacement         = 0.2;
    const double lengthX              = 200;
    //directionX
    NuTo::FullVector<double, dimension> directionX(
    { 1, 0, 0 });
    //directionY
    NuTo::FullVector<double, dimension> directionY(
    { 0, 1, 0 });
    //directionZ
    NuTo::FullVector<double, dimension> directionZ(
    { 0, 0, 1 });

    //**********************************************
    //          Structure
    //**********************************************

    NuTo::Structure myStructure(dimension);

    boost::filesystem::path resultDirectory("/home/phuschke/3D_Gradient_results/");
    boost::filesystem::remove_all(resultDirectory);
    boost::filesystem::create_directory(resultDirectory);
    boost::filesystem::path meshFile("/home/phuschke/develop/nuto/myNutoExamples/beam3D.msh");

    //**********************************************
    //         Integration Scheme
    //**********************************************

    double simulationTime = 1;
    NuTo::NewmarkDirect myIntegrationScheme(&myStructure);
    myIntegrationScheme.SetTimeStep(1e-2 * simulationTime);
    myIntegrationScheme.SetMinTimeStep(1e-5 * simulationTime);
    myIntegrationScheme.SetMaxTimeStep(1e-1 * simulationTime);
    myIntegrationScheme.SetToleranceForce(1e-6);
    myIntegrationScheme.SetAutomaticTimeStepping(true);
    myIntegrationScheme.SetVerboseLevel(0);
    myIntegrationScheme.SetShowTime(false);
    myIntegrationScheme.SetPerformLineSearch(false);
    myIntegrationScheme.SetCheckCoefficientMatrix(false);

    //**********************************************
    //          Section
    //**********************************************

    int mySection = myStructure.SectionCreate(NuTo::Section::VOLUME);

    //**********************************************
    //          Material
    //**********************************************

    NuTo::FullVector<double, Eigen::Dynamic> myDamageLaw(1);
    myDamageLaw(0) = NuTo::Constitutive::eDamageLawType::ISOTROPIC_EXPONENTIAL_SOFTENING;

    int myMaterial = myStructure.ConstitutiveLawCreate(NuTo::Constitutive::eConstitutiveType::GRADIENT_DAMAGE_ENGINEERING_STRESS);
    myStructure.ConstitutiveLawSetPoissonsRatio(myMaterial, poissonsRatio);
    myStructure.ConstitutiveLawSetYoungsModulus(myMaterial, youngsModulus);
    myStructure.ConstitutiveLawSetNonlocalRadius(myMaterial, nonlocalRadius);
    myStructure.ConstitutiveLawSetTensileStrength(myMaterial, tensileStrength);
    myStructure.ConstitutiveLawSetFractureEnergy(myMaterial, fractureEnergy);
    myStructure.ConstitutiveLawSetDamageLaw(myMaterial, myDamageLaw);
    myStructure.ConstitutiveLawSetCompressiveStrength(myMaterial, compressiveStrength);

    int myMaterialweak = myStructure.ConstitutiveLawCreate(NuTo::Constitutive::eConstitutiveType::GRADIENT_DAMAGE_ENGINEERING_STRESS);
    myStructure.ConstitutiveLawSetPoissonsRatio(myMaterialweak, poissonsRatio);
    myStructure.ConstitutiveLawSetYoungsModulus(myMaterialweak, youngsModulus*0.5);
    myStructure.ConstitutiveLawSetNonlocalRadius(myMaterialweak, nonlocalRadius);
    myStructure.ConstitutiveLawSetTensileStrength(myMaterialweak, tensileStrength*0.5);
    myStructure.ConstitutiveLawSetFractureEnergy(myMaterialweak, fractureEnergy);
    myStructure.ConstitutiveLawSetDamageLaw(myMaterialweak, myDamageLaw);
    myStructure.ConstitutiveLawSetCompressiveStrength(myMaterialweak, compressiveStrength*0.5);

    //**********************************************
    //          Geometry
    //**********************************************

    NuTo::FullMatrix<int, Eigen::Dynamic, Eigen::Dynamic> createdGroupIds = myStructure.ImportFromGmsh(meshFile.string(), NuTo::ElementData::CONSTITUTIVELAWIP, NuTo::IpData::eIpDataType::STATICDATA);
    int groupId = createdGroupIds.GetValue(0, 0);

    // triangular elements
    int myInterpolationType = myStructure.InterpolationTypeCreate(NuTo::Interpolation::eShapeType::BRICK3D);
    myStructure.InterpolationTypeAdd(myInterpolationType, NuTo::Node::COORDINATES, NuTo::Interpolation::EQUIDISTANT1);
    myStructure.InterpolationTypeAdd(myInterpolationType, NuTo::Node::DISPLACEMENTS, NuTo::Interpolation::EQUIDISTANT2);
    myStructure.InterpolationTypeAdd(myInterpolationType, NuTo::Node::NONLOCALEQSTRAIN, NuTo::Interpolation::EQUIDISTANT1);
    myStructure.ElementGroupSetInterpolationType(groupId, myInterpolationType);
    myStructure.ElementTotalConvertToInterpolationType(1.e-6, 10);
    myStructure.InterpolationTypeSetIntegrationType(myInterpolationType, NuTo::IntegrationType::IntegrationType3D8NGauss2x2x2Ip, NuTo::IpData::STATICDATA);

    myStructure.ElementTotalSetConstitutiveLaw(myMaterial);
    myStructure.ElementTotalSetSection(mySection);

    int weak_Nodes = myStructure.GroupCreate(NuTo::Groups::Nodes);
    myStructure.GroupAddNodeCoordinateRange(weak_Nodes, 0 , 0.5*lengthX-0.1, 0.5*lengthX+0.1);
    int weak_Elements = myStructure.GroupCreate(NuTo::Groups::Elements);
    myStructure.GroupAddElementsFromNodes(weak_Elements, weak_Nodes, false);
    myStructure.ElementGroupSetConstitutiveLaw(weak_Elements,myMaterialweak);

//    NuTo::FullVector<int, Eigen::Dynamic> weak_Element_Vector = myStructure.GroupGetMemberIds(weak_Elements);
//    for (int i = 0; i < weak_Element_Vector.size(); ++i)
//    {
//        int ele = weak_Element_Vector[i];
//        int numIP = myStructure.ElementGetElementPtr(ele)->GetNumIntegrationPoints();
//
//
//        for (int theIP = 0; theIP < numIP; ++theIP)
//        {
//            myStructure.ElementGetElementPtr(ele)->GetStaticData(theIP)->AsGradientDamage1D()->SetKappa(7e-5);
//        }
//
//
//    }

    //**********************************************
    //          Boundary Conditions
    //**********************************************

    // left boundary
    int nodesLeft = myStructure.GroupCreate(NuTo::Groups::Nodes);
    myStructure.GroupAddNodeCoordinateRange(nodesLeft, 0, 0, 1e-6);
    myStructure.ConstraintLinearSetDisplacementNodeGroup(nodesLeft, directionX, 0);

    // fix single node
    NuTo::FullVector<double, 3> center;
    center.setZero();
    int node_0_0_0 = myStructure.GroupCreate(NuTo::Groups::Nodes);
    myStructure.GroupAddNodeRadiusRange(node_0_0_0,center,0,1e-6);
    myStructure.ConstraintLinearSetDisplacementNodeGroup(node_0_0_0, directionY, 0);
    myStructure.ConstraintLinearSetDisplacementNodeGroup(node_0_0_0, directionZ, 0);

    // fix single node
    center[0]  = lengthX;
    center[1]  = 0;
    center[2]  = 0;
    int node_lenX_0_0 = myStructure.GroupCreate(NuTo::Groups::Nodes);
    myStructure.GroupAddNodeRadiusRange(node_lenX_0_0,center,0,1e-6);
    myStructure.ConstraintLinearSetDisplacementNodeGroup(node_lenX_0_0, directionY, 0);
    myStructure.ConstraintLinearSetDisplacementNodeGroup(node_lenX_0_0, directionZ, 0);

    //**********************************************
    //          Loads
    //**********************************************

    // right load
    int nodesRight = myStructure.GroupCreate(NuTo::Groups::Nodes);
    myStructure.GroupAddNodeCoordinateRange(nodesRight, 0, lengthX - 1e-6, lengthX + 1e-6);
    int load = myStructure.ConstraintLinearSetDisplacementNodeGroup(nodesRight, directionX, 0);

    //**********************************************
    //          Visualisation
    //**********************************************

    myStructure.AddVisualizationComponentDisplacements();
    myStructure.AddVisualizationComponentEngineeringStrain();
    myStructure.AddVisualizationComponentEngineeringStress();
    myStructure.AddVisualizationComponentDamage();
    myStructure.AddVisualizationComponentNonlocalEqStrain();
    myStructure.AddVisualizationComponentConstitutive();
    //**********************************************
    //          Solver
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

    std::cout << " ===> End 3DGradient <=== " << std::endl;
}

