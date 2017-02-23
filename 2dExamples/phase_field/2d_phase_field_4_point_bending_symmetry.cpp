#include "nuto/mechanics/structures/unstructured/Structure.h"
#include "nuto/mechanics/timeIntegration/NewmarkDirect.h"

#include <boost/filesystem.hpp>
#include <iostream>
#include <fstream>
#include <string>
#include <chrono>

using std::cout;
using std::endl;

constexpr unsigned int dimension = 2;
class Parameters
{
public:

    static const int mDimension = dimension;

    static const bool mPerformLineSearch = true;
    static const bool mAutomaticTimeStepping = true;

    static constexpr double mMatrixYoungsModulus = 4.0e4;   // concrete
    static constexpr double mMatrixPoissonsRatio = 0.2;
    static constexpr double mMatrixThickness = 50;
    static constexpr double mMatrixNonlocalRadius = 4;
    static constexpr double mMatrixTensileStrength = 3;
    static constexpr double mMatrixCompressiveStrength = 30;
    static constexpr double mMatrixFractureEnergy = 0.01;


    static constexpr double mTimeStep = 1e-2;
    static constexpr double mMinTimeStep = 1e-5;
    static constexpr double mMaxTimeStep = 1e-1;
    static constexpr double mToleranceForce = 1e-6;
    static constexpr double mSimulationTime = 1.0;
    static constexpr double mLoad = -10*0.04;

    static boost::filesystem::path mOutputPath;
    static const boost::filesystem::path mMeshFilePathMatrix;

    static const NuTo::FullVector<double, dimension> mDirectionX;
    static const NuTo::FullVector<double, dimension> mDirectionY;
};

boost::filesystem::path Parameters::mOutputPath;
const boost::filesystem::path Parameters::mMeshFilePathMatrix("/home/phuschke/meshFiles/2d/fourPointBending.msh");


const NuTo::FullVector<double, dimension> Parameters::mDirectionX = NuTo::FullVector<double, dimension>::UnitX();
const NuTo::FullVector<double, dimension> Parameters::mDirectionY = NuTo::FullVector<double, dimension>::UnitY();


int main(int argc, char* argv[])
{

    if (argc != 5)
    {
        std::cout << "input arguments: displacement order, damage order, integration order, mesh size" << std::endl;
        return EXIT_FAILURE;
    }

    // material
    const int           dispOrder   = std::stoi(argv[1]);
    const int           nlOrder     = std::stoi(argv[2]);
    const int           ipOrder     = std::stoi(argv[3]);
    const int           meshSize    = std::stoi(argv[4]);


    cout << "**********************************************" << endl;
    cout << "**  strucutre                               **" << endl;
    cout << "**********************************************" << endl;

    NuTo::Structure myStructure(dimension);
    myStructure.SetShowTime(false);
    char resultDirectoryName[200];
    sprintf(resultDirectoryName, "/home/phuschke/2d_phase_field_4_point_bending/dispOrder_%02d_damageOrder_%02d_ip_%02d_meshSize_%02d/", dispOrder, nlOrder, ipOrder, meshSize);
    boost::filesystem::path resultDirectory(resultDirectoryName);
    boost::filesystem::remove_all(resultDirectory);
    boost::filesystem::create_directory(resultDirectory);
    Parameters::mOutputPath = resultDirectory;

    cout << "**********************************************" << endl;
    cout << "**  integration sheme                       **" << endl;
    cout << "**********************************************" << endl;

    NuTo::NewmarkDirect myIntegrationScheme(&myStructure);
    myIntegrationScheme.SetTimeStep(Parameters::mTimeStep);
    myIntegrationScheme.SetMinTimeStep(Parameters::mMinTimeStep);
    myIntegrationScheme.SetMaxTimeStep(Parameters::mMaxTimeStep);
    myIntegrationScheme.SetToleranceForce(Parameters::mToleranceForce);
    myIntegrationScheme.SetAutomaticTimeStepping(Parameters::mAutomaticTimeStepping);
    myIntegrationScheme.SetPerformLineSearch(Parameters::mPerformLineSearch);
    myIntegrationScheme.SetResultDirectory(Parameters::mOutputPath.string(), true);


    cout << "**********************************************" << endl;
    cout << "**  section                                 **" << endl;
    cout << "**********************************************" << endl;

    int mySection = myStructure.SectionCreate(NuTo::Section::PLANE_STRESS);
    myStructure.SectionSetThickness(mySection, Parameters::mMatrixThickness);

    cout << "**********************************************" << endl;
    cout << "**  material                                **" << endl;
    cout << "**********************************************" << endl;

    int matrixMaterial = myStructure.ConstitutiveLawCreate(NuTo::Constitutive::eConstitutiveType::PHASE_FIELD);
    myStructure.ConstitutiveLawSetParameterDouble(matrixMaterial, NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, Parameters::mMatrixYoungsModulus);
    myStructure.ConstitutiveLawSetParameterDouble(matrixMaterial, NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO, Parameters::mMatrixPoissonsRatio);
    myStructure.ConstitutiveLawSetParameterDouble(matrixMaterial, NuTo::Constitutive::eConstitutiveParameter::LENGTH_SCALE_PARAMETER, 10);
    myStructure.ConstitutiveLawSetParameterDouble(matrixMaterial, NuTo::Constitutive::eConstitutiveParameter::FRACTURE_ENERGY, 5);

    cout << "**********************************************" << endl;
    cout << "**  geometry                                **" << endl;
    cout << "**********************************************" << endl;

    NuTo::FullMatrix<int, Eigen::Dynamic, Eigen::Dynamic> createdGroupIds = myStructure.ImportFromGmsh(Parameters::mMeshFilePathMatrix.string(), NuTo::ElementData::CONSTITUTIVELAWIP, NuTo::IpData::eIpDataType::STATICDATA);
    int groupId = createdGroupIds.GetValue(0, 0);

    int myInterpolationType = myStructure.InterpolationTypeCreate(NuTo::Interpolation::eShapeType::TRIANGLE2D);
    myStructure.InterpolationTypeAdd(myInterpolationType, NuTo::Node::COORDINATES, NuTo::Interpolation::EQUIDISTANT1);

    switch (dispOrder)
    {
    case 1:
        myStructure.InterpolationTypeAdd(myInterpolationType, NuTo::Node::DISPLACEMENTS, NuTo::Interpolation::EQUIDISTANT1);
        break;
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
        std::cout << "dispOrder either 2,3 or 4." << std::endl;
        return EXIT_FAILURE;
    }

    switch (nlOrder)
    {
    case 1:
        myStructure.InterpolationTypeAdd(myInterpolationType, NuTo::Node::DAMAGE, NuTo::Interpolation::EQUIDISTANT1);
        break;
    default:
        std::cout << "nlOrder either 1,2 or 3." << std::endl;
        return EXIT_FAILURE;
    }

    myStructure.ElementGroupSetInterpolationType(groupId, myInterpolationType);

    switch (ipOrder)
    {
    case 1:
        myStructure.InterpolationTypeSetIntegrationType(myInterpolationType, NuTo::IntegrationType::IntegrationType2D3NGauss1Ip, NuTo::IpData::STATICDATA);
        break;
    case 2:
        myStructure.InterpolationTypeSetIntegrationType(myInterpolationType, NuTo::IntegrationType::IntegrationType2D3NGauss3Ip, NuTo::IpData::STATICDATA);
        break;
    case 3:
        myStructure.InterpolationTypeSetIntegrationType(myInterpolationType, NuTo::IntegrationType::IntegrationType2D3NGauss6Ip, NuTo::IpData::STATICDATA);
        break;
    case 4:
        myStructure.InterpolationTypeSetIntegrationType(myInterpolationType, NuTo::IntegrationType::IntegrationType2D3NGauss12Ip, NuTo::IpData::STATICDATA);
        break;
    default:
        std::cout << "ipOrder either 2, 3 or 4." << std::endl;
        return EXIT_FAILURE;
    }

    myStructure.InterpolationTypeInfo(myInterpolationType);

    myStructure.ElementTotalConvertToInterpolationType(1.e-6, 10);

    myStructure.ElementTotalSetSection(mySection);
    myStructure.ElementTotalSetConstitutiveLaw(matrixMaterial);

    cout << "**********************************************" << endl;
    cout << "**  bc                                      **" << endl;
    cout << "**********************************************" << endl;
    NuTo::FullVector<double, 2> center;
    constexpr double tol = 1.0e-6;
    // bottom left boundary
    center[0] = 25;
    center[1] = 0;
    int grpNodes_bottom_left = myStructure.GroupCreate(NuTo::Groups::Nodes);
    myStructure.GroupAddNodeRadiusRange(grpNodes_bottom_left, center, 0, tol);
    myStructure.ConstraintLinearSetDisplacementNodeGroup(grpNodes_bottom_left, Parameters::mDirectionY, 0);
//    myStructure.ConstraintLinearSetDisplacementNodeGroup(grpNodes_bottom_left, Parameters::mDirectionX, 0);

    int group_nodes_tmp_01 = myStructure.GroupCreate(NuTo::Groups::Nodes);
    myStructure.GroupAddNodeCoordinateRange(group_nodes_tmp_01, 0, 250 - tol, 250 + tol);

    int group_nodes_tmp_02 = myStructure.GroupCreate(NuTo::Groups::Nodes);
    myStructure.GroupAddNodeCoordinateRange(group_nodes_tmp_02, 1, 10 - tol, 100 + tol);


    int group_nodes_symmetry = myStructure.GroupIntersection(group_nodes_tmp_01, group_nodes_tmp_02);






    // remove nodes with only nonlocal eq strains
    int group_only_nonlocal = myStructure.GroupCreate(NuTo::Groups::Nodes);

    auto group_member_ids = myStructure.GroupGetMemberIds(group_nodes_symmetry);
    for (int iNode = 0; iNode < myStructure.GroupGetNumMembers(group_nodes_symmetry); ++iNode)
    {
        int node_id = group_member_ids[iNode];
        auto node_ptr = myStructure.NodeGetNodePtr(node_id);

        if (node_ptr->GetNumDisplacements()==0 and node_ptr->GetNumNonlocalEqStrain()>0)
            myStructure.GroupAddNode(group_only_nonlocal, node_id);

    }


    int group_nodes_test = myStructure.GroupDifference(group_nodes_symmetry, group_only_nonlocal);



    myStructure.ConstraintLinearSetDisplacementNodeGroup(group_nodes_test, Parameters::mDirectionX, 0);

    cout << "**********************************************" << endl;
    cout << "**  load                                    **" << endl;
    cout << "**********************************************" << endl;

    // middle top load
    center[0] = 175;
    center[1] = 100;
    int grpNodes_load = myStructure.GroupCreate(NuTo::Groups::Nodes);

    myStructure.GroupAddNodeRadiusRange(grpNodes_load, center, 0, 1e-6);


    int load = myStructure.ConstraintLinearSetDisplacementNodeGroup(grpNodes_load, Parameters::mDirectionY, 0);

    cout << "**********************************************" << endl;
    cout << "**  visualization                           **" << endl;
    cout << "**********************************************" << endl;

    myStructure.AddVisualizationComponent(groupId, NuTo::VisualizeBase::DISPLACEMENTS);
    myStructure.AddVisualizationComponent(groupId, NuTo::VisualizeBase::ENGINEERING_STRAIN);
    myStructure.AddVisualizationComponent(groupId, NuTo::VisualizeBase::ENGINEERING_STRESS);
    myStructure.AddVisualizationComponent(groupId, NuTo::VisualizeBase::DAMAGE_PHASE_FIELD);

    cout << "**********************************************" << endl;
    cout << "**  solver                                  **" << endl;
    cout << "**********************************************" << endl;


    myStructure.NodeBuildGlobalDofs();
    myStructure.CalculateMaximumIndependentSets();

    myStructure.NodeInfo(10);

//    center[0] = 175;
//    center[1] = 100;
//    int grpNodes_output = myStructure.GroupCreate(NuTo::Groups::Nodes);
//    myStructure.GroupAddNodeRadiusRange(grpNodes_output, center, 0, 5e-1);
//
//    myIntegrationScheme.AddResultGroupNodeForce("myforce", grpNodes_output);
//
//    center[0] = 240;
//    center[1] = 0;
//    int grpNodes_output_disp = myStructure.GroupCreate(NuTo::Groups::Nodes);
//    myStructure.GroupAddNodeRadiusRange(grpNodes_output_disp, center, 0, 7e-1);
//    myIntegrationScheme.AddResultNodeDisplacements("mydisplacements", myStructure.GroupGetMemberIds(grpNodes_output_disp).GetValue(0, 0));

    NuTo::FullMatrix<double, 2, 2> dispRHS;
    dispRHS(0, 0) = 0;
    dispRHS(1, 0) = Parameters::mSimulationTime;
    dispRHS(0, 1) = 0;
    dispRHS(1, 1) = Parameters::mLoad;

    myIntegrationScheme.AddTimeDependentConstraint(load, dispRHS);


    myIntegrationScheme.Solve(Parameters::mSimulationTime);

//    std::string command = "paste " + Parameters::mOutputPath.string() + "myforce.dat " + Parameters::mOutputPath.string() + "mydisplacements.dat > " + Parameters::mOutputPath.string() + "forceDisp.dat";
//            system(command.c_str());
    cout << "**********************************************" << endl;
    cout << "**  end                                     **" << endl;
    cout << "**********************************************" << endl;
}

