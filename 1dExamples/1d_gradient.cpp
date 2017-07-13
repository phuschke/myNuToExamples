#include "nuto/mechanics/structures/unstructured/Structure.h"
#include "nuto/mechanics/timeIntegration/NewmarkDirect.h"
#include "nuto/mechanics/sections/SectionTruss.h"
#include <boost/filesystem.hpp>

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <chrono>

constexpr unsigned int dimension = 1;
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
    static constexpr double mMatrixCompressiveStrength = 30;
    static constexpr double mMatrixFractureEnergy = 0.1;


    static constexpr double mCrossSection = 0.1;
    static constexpr double mTimeStep = 1e-4;
    static constexpr double mMinTimeStep = 1e-7;
    static constexpr double mMaxTimeStep = 1e-2;
    static constexpr double mToleranceForce = 1e-6;
    static constexpr double mSimulationTime = 1.0;
    static constexpr double mLoad = 2;

    static const boost::filesystem::path mOutputPath;

    static const NuTo::FullVector<double, dimension> mDirectionX;
};

const boost::filesystem::path Parameters::mOutputPath("/home/phuschke/1d_gradient/");

const NuTo::FullVector<double, dimension> Parameters::mDirectionX = NuTo::FullVector<double, dimension>::UnitX();


int main(int argc, char* argv[])
{

    if (argc != 5)
    {
        std::cout << "input arguments: numElements, displacement order, nl eq strain order, integration order"
                  << std::endl;
        return 1;
    }

    const int numElements = std::stoi(argv[1]);
    const int dispOrder = std::stoi(argv[2]);
    const int nlOrder = std::stoi(argv[3]);
    const int ipOrder = std::stoi(argv[4]);

    const double length = 100;

    //**********************************************
    //          Structure
    //**********************************************

    NuTo::Structure myStructure(dimension);
    myStructure.SetShowTime(false);
    char resultDirectoryName[200];
    sprintf(resultDirectoryName, "/home/phuschke/1D_Gradient_results_ele_%04d_dispOrder_%02d_nlOrder_%02d_ip_%02d/",
            numElements, dispOrder, nlOrder, ipOrder);
    boost::filesystem::path resultDirectory(resultDirectoryName);
    boost::filesystem::remove_all(resultDirectory);
    boost::filesystem::create_directory(resultDirectory);

    //**********************************************
    //         Integration Scheme
    //**********************************************

    NuTo::NewmarkDirect myIntegrationScheme(&myStructure);
    myIntegrationScheme.SetTimeStep(Parameters::mTimeStep);
    myIntegrationScheme.SetMinTimeStep(Parameters::mMinTimeStep);
    myIntegrationScheme.SetMaxTimeStep(Parameters::mMaxTimeStep);
    myIntegrationScheme.SetToleranceForce(Parameters::mToleranceForce);
    myIntegrationScheme.SetAutomaticTimeStepping(Parameters::mAutomaticTimeStepping);
    myIntegrationScheme.SetVerboseLevel(0);
    myIntegrationScheme.SetShowTime(false);
    myIntegrationScheme.SetPerformLineSearch(Parameters::mPerformLineSearch);
    myIntegrationScheme.SetCheckCoefficientMatrix(false);

    //**********************************************
    //          Section
    //**********************************************

    double alpha = 0.1;
    double xWeakSpot = 50.;
    double lWeakSpot = 20.;
    double areaParameters[4];
    areaParameters[0] = xWeakSpot;
    areaParameters[1] = lWeakSpot;
    areaParameters[2] = alpha;
    areaParameters[3] = 1;

    int mySection = myStructure.SectionCreate(NuTo::Section::TRUSS);
    myStructure.SectionSetArea(mySection, Parameters::mCrossSection);

    //    int weakSection = myStructure.SectionCreate(NuTo::Section::TRUSS);
    //    myStructure.SectionSetArea(weakSection, 50);

    myStructure.SectionGetSectionPtr(mySection)->AsSectionTruss()->SetAreaParameters(areaParameters);

    //**********************************************
    //          Material
    //**********************************************

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

    // create nodes
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
        myStructure.InterpolationTypeAdd(myInterpolationType, NuTo::Node::DISPLACEMENTS,
                                         NuTo::Interpolation::EQUIDISTANT2);
        break;
    case 3:
        myStructure.InterpolationTypeAdd(myInterpolationType, NuTo::Node::DISPLACEMENTS,
                                         NuTo::Interpolation::EQUIDISTANT3);
        break;
    case 4:
        myStructure.InterpolationTypeAdd(myInterpolationType, NuTo::Node::DISPLACEMENTS,
                                         NuTo::Interpolation::EQUIDISTANT4);
        break;
    default:
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
    case 3:
        myStructure.InterpolationTypeAdd(myInterpolationType, NuTo::Node::NONLOCALEQSTRAIN,
                                         NuTo::Interpolation::EQUIDISTANT3);
        break;
    default:
        return 1;
    }

    switch (ipOrder)
    {
    case 2:
        myStructure.InterpolationTypeSetIntegrationType(
                myInterpolationType, NuTo::IntegrationType::IntegrationType1D2NGauss2Ip, NuTo::IpData::STATICDATA);
        break;
    case 3:
        myStructure.InterpolationTypeSetIntegrationType(
                myInterpolationType, NuTo::IntegrationType::IntegrationType1D2NGauss3Ip, NuTo::IpData::STATICDATA);
        break;
    case 4:
        myStructure.InterpolationTypeSetIntegrationType(
                myInterpolationType, NuTo::IntegrationType::IntegrationType1D2NGauss4Ip, NuTo::IpData::STATICDATA);
        break;
    case 5:
        myStructure.InterpolationTypeSetIntegrationType(
                myInterpolationType, NuTo::IntegrationType::IntegrationType1D2NGauss5Ip, NuTo::IpData::STATICDATA);
        break;
    default:
        return 1;
    }

    // create elements
    for (int countX = 0; countX < numElements; countX++)
    {
        NuTo::FullVector<int, Eigen::Dynamic> nodes(2);
        nodes(0) = countX;
        nodes(1) = countX + 1;
        myStructure.ElementCreate(myInterpolationType, nodes, NuTo::ElementData::eElementDataType::CONSTITUTIVELAWIP,
                                  NuTo::IpData::eIpDataType::STATICDATA);
    }

    myStructure.ElementTotalConvertToInterpolationType(1.e-6, 10);

    myStructure.ElementTotalSetConstitutiveLaw(matrixMaterial);
    myStructure.ElementTotalSetSection(mySection);

    int weakNodes = myStructure.GroupCreate(NuTo::Groups::Nodes);
    myStructure.GroupAddNodeCoordinateRange(weakNodes, 0, 45 - 1e-6, 55 + 1e-6);
    int weakElemets = myStructure.GroupCreate(NuTo::Groups::Elements);
    myStructure.GroupAddElementsFromNodes(weakElemets, weakNodes, true);

    // myStructure.ElementGroupSetSection(weakElemets, weakSection);

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

    //    myStructure.AddVisualizationComponentDisplacements();
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

    myIntegrationScheme.AddResultGroupNodeForce("myforce", nodesRight);
    myIntegrationScheme.AddResultNodeDisplacements("mydisplacements",
                                                   myStructure.GroupGetMemberIds(nodesRight).GetValue(0, 0));

    myIntegrationScheme.AddTimeDependentConstraint(load, timeDependentLoad);

    myIntegrationScheme.SetResultDirectory(resultDirectory.string(), false);
    myStructure.SetVerboseLevel(10);
    myStructure.Info();

    auto begin = std::chrono::high_resolution_clock::now();
    myIntegrationScheme.Solve(Parameters::mSimulationTime);
    auto end = std::chrono::high_resolution_clock::now();


    //**********************************************
    //          Postprocess
    //**********************************************


    myStructure.ExportVtkDataFileElements(resultDirectory.string() + "out.vtk");

    // output stress at ip
    boost::filesystem::path stressFile(resultDirectory.string() + "stressFile.dat");
    std::ofstream stressOutput;
    stressOutput.open(stressFile.string());

    NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> stress;
    NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> coords;
    double maxStress = 0.;
    myStructure.Info();
    for (int ele = 0; ele < myStructure.GetNumElements(); ++ele)
    {
        stress = myStructure.ElementGetEngineeringStress(ele);
        coords = myStructure.ElementGetIntegrationPointCoordinates(ele);

        int numIP = myStructure.ElementGetElementPtr(ele)->GetNumIntegrationPoints();
        for (int ip = 0; ip < numIP; ++ip)
        {
            stressOutput << stress.GetValue(0, ip) << "\t" << coords.GetValue(0, ip) << std::endl;
            maxStress = std::max(maxStress, stress.GetValue(0, ip));
        }
    }
    stressOutput.close();

    // output damage at ip
    boost::filesystem::path damageFile(resultDirectory.string() + "damageFile.dat");
    std::ofstream damageOutput;
    damageOutput.open(damageFile.string());

    NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> damage;

    myStructure.Info();
    for (int ele = 0; ele < myStructure.GetNumElements(); ++ele)
    {
        coords = myStructure.ElementGetIntegrationPointCoordinates(ele);
        damage = myStructure.ElementGetDamage(ele);

        int numIP = myStructure.ElementGetElementPtr(ele)->GetNumIntegrationPoints();
        for (int ip = 0; ip < numIP; ++ip)
        {
            damageOutput << damage.GetValue(0, ip) << "\t" << coords.GetValue(0, ip) << std::endl;
        }
    }
    damageOutput.close();
    int dofs = myStructure.GetNumDofs(NuTo::Node::eDof::DISPLACEMENTS);
    dofs += myStructure.GetNumDofs(NuTo::Node::eDof::NONLOCALEQSTRAIN);
    boost::filesystem::path maxStressFile(resultDirectory.string() + "maxStress.dat");
    std::ofstream maxStressOutput;
    maxStressOutput.open(maxStressFile.string(), std::ios_base::app);
    maxStressOutput << numElements << "\t" << dispOrder << "\t" << nlOrder << "\t" << ipOrder << "\t" << dofs << "\t"
                    << maxStress << std::endl;
    maxStressOutput.close();

    char forceDispName[200];
    sprintf(forceDispName, "forceDisp_ele_%03d_disp_%02d_nl_%02d_ip_%02d.dat", numElements, dispOrder, nlOrder,
            ipOrder);
    std::string command = "paste " + resultDirectory.string() + "myforce.dat " + resultDirectory.string() +
                          "mydisplacements.dat > " + resultDirectory.string() + forceDispName;
    system(command.c_str());

    char stressName[200];
    sprintf(stressName, "stress_ele_%03d_disp_%02d_nl_%02d_ip_%02d.dat", numElements, dispOrder, nlOrder, ipOrder);
    command = "mv " + resultDirectory.string() + "stressFile.dat " + resultDirectory.string() + stressName;
    system(command.c_str());

    // calculate integral
    boost::filesystem::path forcePath(resultDirectory.string() + "myforce.dat");
    std::ifstream inputFile;

    inputFile.open(forcePath.string(), std::ifstream::in);
    std::vector<double> force;
    while (inputFile.good())
    {
        char temp[256];
        inputFile.getline(temp, 256);
        force.push_back(std::stod(temp));
    }
    inputFile.close();

    boost::filesystem::path dispPath(resultDirectory.string() + "mydisplacements.dat");
    inputFile.open(dispPath.string(), std::ifstream::in);
    std::vector<double> disp;
    while (inputFile.good())
    {
        char temp[256];
        inputFile.getline(temp, 256);
        disp.push_back(std::stod(temp));
    }
    inputFile.close();

    double integral = 0.;
    assert(force.size() == disp.size());
    for (unsigned int i = 0; i < force.size() - 1; ++i)
    {
        integral += 0.5 * (disp[i + 1] - disp[i]) * (force[i + 1] + force[i]);
    }

    // integral
    boost::filesystem::path IntegralFile("/home/phuschke/1d_gradient/Integral.dat");
    std::ofstream IntegralOutput;
    IntegralOutput.open(IntegralFile.string(), std::ios_base::app);
    IntegralOutput << numElements << "\t\t" << dispOrder << "\t\t" << nlOrder << "\t\t" << ipOrder << "\t\t" << dofs
                   << "\t\t" << std::setprecision(10) << integral << "\t\t"
                   << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << std::endl;
    IntegralOutput.close();

    std::cout << " ===> End 1DGradient <=== " << std::endl;
}
