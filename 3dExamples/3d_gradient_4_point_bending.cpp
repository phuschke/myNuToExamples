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


constexpr unsigned int dimension = 3;
class Parameters
{
public:
    static const int mDimension = dimension;

    static const bool mPerformLineSearch = true;
    static const bool mAutomaticTimeStepping = true;

    static constexpr double mMatrixYoungsModulus = 3.0e4; // concrete
    static constexpr double mMatrixPoissonsRatio = 0.2;
    static constexpr double mMatrixThickness = 0.1;
    static constexpr double mMatrixNonlocalRadius = 2;
    static constexpr double mMatrixTensileStrength = 3;
    static constexpr double mMatrixCompressiveStrength = 50;
    static constexpr double mMatrixFractureEnergy = 0.1;


    static constexpr double mTimeStep = 1e-2;
    static constexpr double mMinTimeStep = 1e-7;
    static constexpr double mMaxTimeStep = 1e-1;
    static constexpr double mToleranceForce = 1e-6;
    static constexpr double mSimulationTime = 1.0;
    static constexpr double mLoad = 0.2;


    static const boost::filesystem::path mOutputPath;
    static const boost::filesystem::path mMeshFilePath;

    static const NuTo::FullVector<double, dimension> mDirectionX;
    static const NuTo::FullVector<double, dimension> mDirectionY;
    static const NuTo::FullVector<double, dimension> mDirectionZ;
};

const boost::filesystem::path Parameters::mOutputPath("/home/phuschke/develop/nuto/3d_gradient_4_point_bending/");
const boost::filesystem::path
        Parameters::mMeshFilePath("/home/phuschke/develop/nuto/myNutoExamples/MeshFiles/3d/3DfourPointBending.msh");

const NuTo::FullVector<double, dimension> Parameters::mDirectionX({1, 0, 0});
const NuTo::FullVector<double, dimension> Parameters::mDirectionY({0, 1, 0});
const NuTo::FullVector<double, dimension> Parameters::mDirectionZ({0, 0, 1});

int main(int argc, char* argv[])
{

    //    if (argc != 4)
    //    {
    //        std::cout << "input arguments: displacement order, nl eq strain order, integration order" << std::endl;
    //        return 1;
    //    }

    // material
    const int dispOrder = 2; // std::stoi(argv[1]);
    const int nlOrder = 1; // std::stoi(argv[2]);
    const int ipOrder = 2; // std::stoi(argv[3]);


    //**********************************************
    //          Structure
    //**********************************************

    NuTo::Structure myStructure(Parameters::mDimension);

    //**********************************************
    //         Integration Scheme
    //**********************************************

    NuTo::NewmarkDirect myIntegrationScheme(&myStructure);
    myIntegrationScheme.SetTimeStep(Parameters::mTimeStep);
    myIntegrationScheme.SetMinTimeStep(Parameters::mMinTimeStep);
    myIntegrationScheme.SetMaxTimeStep(Parameters::mMaxTimeStep);
    myIntegrationScheme.SetToleranceForce(Parameters::mToleranceForce);
    myIntegrationScheme.SetAutomaticTimeStepping(Parameters::mAutomaticTimeStepping);
    myIntegrationScheme.SetPerformLineSearch(Parameters::mPerformLineSearch);

    //**********************************************
    //          Section
    //**********************************************

    int mySection = myStructure.SectionCreate(NuTo::Section::VOLUME);

    //**********************************************
    //          Material
    //**********************************************

    NuTo::FullVector<double, Eigen::Dynamic> myDamageLaw(1);
    myDamageLaw(0) = NuTo::Constitutive::eDamageLawType::ISOTROPIC_EXPONENTIAL_SOFTENING_RES_LOAD;

    int matrixMaterial = myStructure.ConstitutiveLawCreate(
            NuTo::Constitutive::eConstitutiveType::GRADIENT_DAMAGE_ENGINEERING_STRESS);
    myStructure.ConstitutiveLawSetParameterDouble(matrixMaterial,
                                                  NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS,
                                                  Parameters::mMatrixYoungsModulus);
    myStructure.ConstitutiveLawSetParameterDouble(matrixMaterial,
                                                  NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO,
                                                  Parameters::mMatrixPoissonsRatio);
    myStructure.ConstitutiveLawSetParameterDouble(matrixMaterial,
                                                  NuTo::Constitutive::eConstitutiveParameter::TENSILE_STRENGTH,
                                                  Parameters::mMatrixTensileStrength);
    myStructure.ConstitutiveLawSetParameterDouble(matrixMaterial,
                                                  NuTo::Constitutive::eConstitutiveParameter::COMPRESSIVE_STRENGTH,
                                                  Parameters::mMatrixCompressiveStrength);
    myStructure.ConstitutiveLawSetParameterDouble(matrixMaterial,
                                                  NuTo::Constitutive::eConstitutiveParameter::NONLOCAL_RADIUS,
                                                  Parameters::mMatrixNonlocalRadius);
    myStructure.ConstitutiveLawSetParameterDouble(matrixMaterial,
                                                  NuTo::Constitutive::eConstitutiveParameter::FRACTURE_ENERGY,
                                                  Parameters::mMatrixFractureEnergy);

    //**********************************************
    //          Geometry
    //**********************************************

    auto createdGroupIds =
            myStructure.ImportFromGmsh(Parameters::mMeshFilePath.string(), NuTo::ElementData::CONSTITUTIVELAWIP,
                                       NuTo::IpData::eIpDataType::STATICDATA);
    int groupId = createdGroupIds.GetValue(0, 0);

    int myInterpolationType = myStructure.InterpolationTypeCreate(NuTo::Interpolation::eShapeType::TETRAHEDRON3D);
    myStructure.InterpolationTypeAdd(myInterpolationType, NuTo::Node::COORDINATES, NuTo::Interpolation::EQUIDISTANT1);

    switch (dispOrder)
    {
    case 2:
        myStructure.InterpolationTypeAdd(myInterpolationType, NuTo::Node::DISPLACEMENTS,
                                         NuTo::Interpolation::EQUIDISTANT2);
        break;
    default:
        std::cout << "dispOrder 2." << std::endl;
        return 1;
    }

    switch (nlOrder)
    {
    case 1:
        myStructure.InterpolationTypeAdd(myInterpolationType, NuTo::Node::NONLOCALEQSTRAIN,
                                         NuTo::Interpolation::EQUIDISTANT1);
        break;
    case 2:
        myStructure.InterpolationTypeAdd(myInterpolationType, NuTo::Node::NONLOCALEQSTRAIN,
                                         NuTo::Interpolation::EQUIDISTANT2);
        break;
    default:
        std::cout << "nlOrder either 1 or 2." << std::endl;
        return 1;
    }

    myStructure.ElementGroupSetInterpolationType(groupId, myInterpolationType);

    switch (ipOrder)
    {
    case 1:
        myStructure.InterpolationTypeSetIntegrationType(
                myInterpolationType, NuTo::IntegrationType::IntegrationType3D4NGauss1Ip, NuTo::IpData::STATICDATA);
        break;
    case 2:
        myStructure.InterpolationTypeSetIntegrationType(
                myInterpolationType, NuTo::IntegrationType::IntegrationType3D4NGauss4Ip, NuTo::IpData::STATICDATA);
        break;
    default:
        std::cout << "ipOrder either 1 or 2." << std::endl;
        return 1;
    }

    myStructure.ElementTotalConvertToInterpolationType(1.e-6, 10);
    myStructure.ElementTotalSetSection(mySection);
    myStructure.ElementTotalSetConstitutiveLaw(matrixMaterial);

    //**********************************************
    //          Boundary Conditions
    //**********************************************
    int groupNodesBottom = myStructure.GroupCreate(NuTo::Groups::Nodes);
    myStructure.GroupAddNodeCoordinateRange(groupNodesBottom, 1, 0 - 1e-6, 0 + 1e-6);

    int groupNodesTop = myStructure.GroupCreate(NuTo::Groups::Nodes);
    myStructure.GroupAddNodeCoordinateRange(groupNodesTop, 1, 100 - 1e-6, 100 + 1e-6);


    // bottom right boundary
    int groupNodesTmp = myStructure.GroupCreate(NuTo::Groups::Nodes);
    myStructure.GroupAddNodeCoordinateRange(groupNodesTmp, 0, 475 - 1e-6, 475 + 1e-6);

    int groupNodesBcBottomRight = myStructure.GroupIntersection(groupNodesBottom, groupNodesTmp);
    myStructure.ConstraintLinearSetDisplacementNodeGroup(groupNodesBcBottomRight, Parameters::mDirectionY, 0);
    myStructure.ConstraintLinearSetDisplacementNodeGroup(groupNodesBcBottomRight, Parameters::mDirectionZ, 0);

    myStructure.GroupDelete(groupNodesTmp);


    // bottom left boundary
    groupNodesTmp = myStructure.GroupCreate(NuTo::Groups::Nodes);
    myStructure.GroupAddNodeCoordinateRange(groupNodesTmp, 0, 25 - 1e-6, 25 + 1e-6);

    int groupNodesBcBottomLeft = myStructure.GroupIntersection(groupNodesBottom, groupNodesTmp);
    myStructure.ConstraintLinearSetDisplacementNodeGroup(groupNodesBcBottomLeft, Parameters::mDirectionX, 0);
    myStructure.ConstraintLinearSetDisplacementNodeGroup(groupNodesBcBottomLeft, Parameters::mDirectionY, 0);
    myStructure.ConstraintLinearSetDisplacementNodeGroup(groupNodesBcBottomLeft, Parameters::mDirectionZ, 0);


    myStructure.GroupDelete(groupNodesTmp);


    //**********************************************
    //          Loads
    //**********************************************

    // middle top boundary
    groupNodesTmp = myStructure.GroupCreate(NuTo::Groups::Nodes);
    myStructure.GroupAddNodeCoordinateRange(groupNodesTmp, 0, 175 - 1e-6, 175 + 1e-6);

    int groupNodesLoadLeft = myStructure.GroupIntersection(groupNodesTop, groupNodesTmp);

    myStructure.GroupDelete(groupNodesTmp);

    groupNodesTmp = myStructure.GroupCreate(NuTo::Groups::Nodes);
    myStructure.GroupAddNodeCoordinateRange(groupNodesTmp, 0, 325 - 1e-6, 325 + 1e-6);

    int groupNodesLoadRight = myStructure.GroupIntersection(groupNodesTop, groupNodesTmp);

    myStructure.GroupDelete(groupNodesTmp);

    int groupNodesLoadTotal = myStructure.GroupUnion(groupNodesLoadLeft, groupNodesLoadRight);

    int load = myStructure.ConstraintLinearSetDisplacementNodeGroup(groupNodesLoadTotal, -Parameters::mDirectionY, 0);

    //**********************************************
    //          Visualisation
    //**********************************************

    myStructure.AddVisualizationComponentDisplacements();
    myStructure.AddVisualizationComponentEngineeringStrain();
    myStructure.AddVisualizationComponentEngineeringStress();
    myStructure.AddVisualizationComponentDamage();
    myStructure.AddVisualizationComponentNonlocalEqStrain();

    //**********************************************
    //          Solver
    //**********************************************

    const int numElements = myStructure.GetNumElements();

    char resultDirectoryName[200];
    sprintf(resultDirectoryName, "/home/phuschke/3D_Gradient_Convergence_ele_%04d_dispOrder_%02d_nlOrder_%02d_ip_%02d/",
            numElements, dispOrder, nlOrder, ipOrder);
    boost::filesystem::path resultDirectory(resultDirectoryName);
    boost::filesystem::remove_all(resultDirectory);
    boost::filesystem::create_directory(resultDirectory);

    myStructure.NodeBuildGlobalDofs();
    myStructure.CalculateMaximumIndependentSets();

    NuTo::FullVector<double, Parameters::mDimension> coordinates;
    coordinates[0] = 175;
    coordinates[1] = 100;
    coordinates[2] = 25;

    int groupNodesOutput = myStructure.GroupCreate(NuTo::Groups::Nodes);
    myStructure.GroupAddNodeRadiusRange(groupNodesOutput, coordinates, 0, 5e-1);

    myIntegrationScheme.AddResultGroupNodeForce("myforce", groupNodesOutput);
    myIntegrationScheme.AddResultNodeDisplacements("mydisplacements",
                                                   myStructure.GroupGetMemberIds(groupNodesOutput).GetValue(0, 0));

    NuTo::FullMatrix<double, 2, 2> dispRHS;
    dispRHS(0, 0) = 0;
    dispRHS(1, 0) = Parameters::mSimulationTime;
    dispRHS(0, 1) = 0;
    dispRHS(1, 1) = Parameters::mLoad;

    myIntegrationScheme.SetTimeDependentConstraint(load, dispRHS);

    bool deleteDirectory = false;
    myIntegrationScheme.SetResultDirectory(resultDirectory.string(), deleteDirectory);

    myIntegrationScheme.Solve(Parameters::mSimulationTime);
    myStructure.ExportVtkDataFileElements(resultDirectory.string() + "out.vtk");

    char forceDispName[200];
    sprintf(forceDispName, "/home/phuschke/develop/nuto/3d_gradient_4_point_bending/"
                           "3D_Gradient_forceDisp_ele_%03d_disp_%02d_nl_%02d_ip_%02d.dat",
            numElements, dispOrder, nlOrder, ipOrder);
    std::string command = "paste " + resultDirectory.string() + "myforce.dat " + resultDirectory.string() +
                          "mydisplacements.dat > " + forceDispName;
    system(command.c_str());

    // calculate integral
    boost::filesystem::path forcePath(resultDirectory.string() + "myforce.dat");
    std::ifstream inputFile;

    inputFile.open(forcePath.string(), std::ifstream::in);
    std::vector<double> force;
    while (inputFile.good())
    {
        char temp[256];
        inputFile >> temp;
        inputFile >> temp;
        force.push_back(std::stod(temp));
        inputFile >> temp;
    }

    inputFile.close();

    boost::filesystem::path dispPath(resultDirectory.string() + "mydisplacements.dat");
    inputFile.open(dispPath.string(), std::ifstream::in);
    std::vector<double> disp;
    while (inputFile.good())
    {
        char temp[256];
        inputFile >> temp;
        inputFile >> temp;
        disp.push_back(std::stod(temp));
        inputFile >> temp;
    }
    inputFile.close();

    double integral = 0.;
    assert(force.size() == disp.size());
    for (unsigned int i = 0; i < force.size() - 1; ++i)
    {
        integral += 0.5 * (disp[i + 1] - disp[i]) * (force[i + 1] + force[i]);
    }

    boost::filesystem::path IntegralFile(Parameters::mOutputPath.string() + "3D_Gradient_Integral.dat");
    std::ofstream IntegralOutput;
    IntegralOutput.open(IntegralFile.string(), std::ios_base::app);

    IntegralOutput << numElements << "\t" << dispOrder << "\t" << nlOrder << "\t" << ipOrder << "\t"
                   << myStructure.GetNumDofs() << "\t" << integral << std::endl;
    IntegralOutput.close();

    std::cout << " ===> End Four Point Bending Test <=== " << std::endl;
}
