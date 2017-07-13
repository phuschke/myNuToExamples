
#include "mechanics/structures/unstructured/Structure.h"
#include "mechanics/timeIntegration/NewmarkDirect.h"
#include "mechanics/constitutive/laws/PhaseField.h"
#include "mechanics/constitutive/ConstitutiveEnum.h"
#include "mechanics/elements/ElementBase.h"
#include "mechanics/interpolationtypes/InterpolationType.h"
#include "mechanics/integrationtypes/IntegrationTypeBase.h"
#include "mechanics/constitutive/inputoutput/ConstitutiveCalculateStaticData.h"

#include "../../EnumsAndTypedefs.h"
#include <boost/filesystem.hpp>
#include <chrono>

using std::endl;
using NuTo::Constitutive::ePhaseFieldEnergyDecomposition;


int main(int argc, char* argv[])
{

    if (argc != 5)
    {
        cout << "input arguments: displacement order, damage order, integration order, subdirectory" << endl;
        return EXIT_FAILURE;
    }

    const int dispOrder = std::stoi(argv[1]);
    const int phaseFieldOrder = std::stoi(argv[2]);
    const int ipOrder = std::stoi(argv[3]);
    const int subdirectory = std::stoi(argv[4]);
    constexpr unsigned dimension = 2;
    constexpr bool performLineSearch = false;
    constexpr bool automaticTimeStepping = true;
    constexpr double youngsModulus = 2.1e5; // N/mm^2
    constexpr double poissonsRatio = 0.3;
    constexpr double thickness = 1.0; // mm
    constexpr double lengthScaleParameter = 1.5e-2; // mm
    constexpr double fractureEnergy = 2.7; // N/mm
    constexpr double artificialViscosity = 1.0e-3; // Ns/mm^2

    constexpr ePhaseFieldEnergyDecomposition energyDecomposition =
            ePhaseFieldEnergyDecomposition::ANISOTROPIC_SPECTRAL_DECOMPOSITION;


    constexpr double timeStep = 1.0e-3;
    constexpr double minTimeStep = 1.0e-6;
    constexpr double maxTimeStep = 1.0e-3;
    constexpr double simulationTime = 1.6e-2;

    constexpr double toleranceDisp = 1.0e-5;
    constexpr double toleranceCrack = 1.0e-4;
    constexpr double loadFactor = simulationTime;

    constexpr double tol = 1.0e-6;

    boost::filesystem::path resultPath(
            std::string("/home/phuschke/results/2d/2d_miehe_single_edge_notched_shear_test_phase_field/" +
                        std::to_string(subdirectory) + "/"));
    const boost::filesystem::path meshFilePath("/home/phuschke/meshFiles/2d/2d_single_edge_notched_shear_test.msh");

    const Eigen::Vector2d directionX = Eigen::Vector2d::UnitX();
    const Eigen::Vector2d directionY = Eigen::Vector2d::UnitY();


    cout << "**********************************************" << endl;
    cout << "**  strucutre                               **" << endl;
    cout << "**********************************************" << endl;

    NuTo::Structure myStructure(dimension);
    myStructure.SetNumTimeDerivatives(1);
    myStructure.SetShowTime(false);

    cout << subdirectory << "subdirectory" << endl;
    cout << resultPath.string() << endl;
    boost::filesystem::create_directory(resultPath);


    cout << "**********************************************" << endl;
    cout << "**  integration sheme                       **" << endl;
    cout << "**********************************************" << endl;

    NuTo::NewmarkDirect myIntegrationScheme(&myStructure);

    myIntegrationScheme.SetTimeStep(timeStep);
    myIntegrationScheme.SetMinTimeStep(minTimeStep);
    myIntegrationScheme.SetMaxTimeStep(maxTimeStep);
    myIntegrationScheme.SetAutomaticTimeStepping(automaticTimeStepping);
    myIntegrationScheme.SetPerformLineSearch(performLineSearch);
    myIntegrationScheme.SetResultDirectory(resultPath.string(), true);
    myIntegrationScheme.SetToleranceResidual(NuTo::Node::eDof::CRACKPHASEFIELD, toleranceCrack);
    myIntegrationScheme.SetToleranceResidual(NuTo::Node::eDof::DISPLACEMENTS, toleranceDisp);
    myIntegrationScheme.SetExportDataFileNodes(false);

    std::ofstream file;
    file.open(std::string(resultPath.string() + "parameters.txt"));
    cout << "Writing simulation parameters to file: " << std::string(resultPath.string() + "parameters.txt")
         << std::endl;
    file << "dispOrder            " << dispOrder << std::endl;
    file << "phaseFieldOrder      " << phaseFieldOrder << std::endl;
    file << "ipOrder              " << ipOrder << std::endl;
    file << "subdirectory         " << subdirectory << std::endl;
    file << "dimension            " << dimension << std::endl;
    file << "performLineSearch    " << performLineSearch << std::endl;
    file << "automaticTimeStepping" << automaticTimeStepping << std::endl;
    file << "youngsModulus        " << youngsModulus << std::endl;
    file << "poissonsRatio        " << poissonsRatio << std::endl;
    file << "thickness            " << thickness << std::endl;
    file << "lengthScaleParameter " << lengthScaleParameter << std::endl;
    file << "fractureEnergy       " << fractureEnergy << std::endl;
    file << "artificialViscosity  " << artificialViscosity << std::endl;
    file << "timeStep             " << timeStep << std::endl;
    file << "minTimeStep          " << minTimeStep << std::endl;
    file << "maxTimeStep          " << maxTimeStep << std::endl;
    file << "toleranceDisp        " << toleranceDisp << std::endl;
    file << "toleranceCrack       " << toleranceCrack << std::endl;
    file << "simulationTime       " << simulationTime << std::endl;
    file << "loadFactor           " << loadFactor << std::endl;

    file.close();

    cout << "**********************************************" << endl;
    cout << "**  section                                 **" << endl;
    cout << "**********************************************" << endl;

    int mySection = myStructure.SectionCreate(NuTo::eSectionType::PLANE_STRAIN);
    myStructure.SectionSetThickness(mySection, thickness);


    cout << "**********************************************" << endl;
    cout << "**  material                                **" << endl;
    cout << "**********************************************" << endl;

    NuTo::ConstitutiveBase* phaseField = new NuTo::PhaseField(youngsModulus, poissonsRatio, lengthScaleParameter,
                                                              fractureEnergy, artificialViscosity, energyDecomposition);

    int matrixMaterial = myStructure.AddConstitutiveLaw(phaseField);

    cout << "**********************************************" << endl;
    cout << "**  geometry                                **" << endl;
    cout << "**********************************************" << endl;

    auto createdGroupIds = myStructure.ImportFromGmsh(meshFilePath.string());
    int groupId = createdGroupIds[0].first;

    int myInterpolationType = myStructure.InterpolationTypeCreate(eShapeType::TRIANGLE2D);
    myStructure.InterpolationTypeAdd(myInterpolationType, eDof::COORDINATES, eTypeOrder::EQUIDISTANT1);

    switch (dispOrder)
    {
    case 1:
        myStructure.InterpolationTypeAdd(myInterpolationType, eDof::DISPLACEMENTS, eTypeOrder::EQUIDISTANT1);
        break;
    case 2:
        myStructure.InterpolationTypeAdd(myInterpolationType, eDof::DISPLACEMENTS, eTypeOrder::EQUIDISTANT2);
        break;
    case 3:
        myStructure.InterpolationTypeAdd(myInterpolationType, eDof::DISPLACEMENTS, eTypeOrder::EQUIDISTANT3);
        break;
    case 4:
        myStructure.InterpolationTypeAdd(myInterpolationType, eDof::DISPLACEMENTS, eTypeOrder::EQUIDISTANT4);
        break;
    default:
        std::cout << "dispOrder either 2,3 or 4." << std::endl;
        return EXIT_FAILURE;
    }

    switch (phaseFieldOrder)
    {
    case 1:
        myStructure.InterpolationTypeAdd(myInterpolationType, eDof::CRACKPHASEFIELD, eTypeOrder::EQUIDISTANT1);
        break;
    default:
        std::cout << "Crack phase-field order either 1,2 or 3." << std::endl;
        return EXIT_FAILURE;
    }

    myStructure.ElementGroupSetInterpolationType(groupId, myInterpolationType);

    switch (ipOrder)
    {
    case 1:
        myStructure.InterpolationTypeSetIntegrationType(myInterpolationType,
                                                        NuTo::eIntegrationType::IntegrationType2D3NGauss1Ip);
        break;
    case 2:
        myStructure.InterpolationTypeSetIntegrationType(myInterpolationType,
                                                        NuTo::eIntegrationType::IntegrationType2D3NGauss3Ip);
        break;
    case 3:
        myStructure.InterpolationTypeSetIntegrationType(myInterpolationType,
                                                        NuTo::eIntegrationType::IntegrationType2D3NGauss6Ip);
        break;
    case 4:
        myStructure.InterpolationTypeSetIntegrationType(myInterpolationType,
                                                        NuTo::eIntegrationType::IntegrationType2D3NGauss12Ip);
        break;
    default:
        std::cout << "ipOrder either 2, 3 or 4." << std::endl;
        return EXIT_FAILURE;
    }

    myStructure.InterpolationTypeInfo(myInterpolationType);

    myStructure.ElementTotalConvertToInterpolationType(1.e-9, 1);

    myStructure.ElementTotalSetSection(mySection);
    myStructure.ElementTotalSetConstitutiveLaw(matrixMaterial);


    cout << "**********************************************" << endl;
    cout << "**  bc                                      **" << endl;
    cout << "**********************************************" << endl;


    NuTo::FullVector<double, 2> center;

    // constrain top boundary in Y
    int groupNodesTop = myStructure.GroupCreate(NuTo::Groups::Nodes);
    myStructure.GroupAddNodeCoordinateRange(groupNodesTop, 1, 1. - tol, 1. + tol);

    myStructure.ConstraintLinearSetDisplacementNodeGroup(groupNodesTop, directionY, 0);

    // constrain bottom boundary in X and Y
    int groupNodesBottom = myStructure.GroupCreate(NuTo::Groups::Nodes);
    myStructure.GroupAddNodeCoordinateRange(groupNodesBottom, 1, -tol, tol);

    myStructure.ConstraintLinearSetDisplacementNodeGroup(groupNodesBottom, directionX, 0);
    myStructure.ConstraintLinearSetDisplacementNodeGroup(groupNodesBottom, directionY, 0);

    cout << "**********************************************" << endl;
    cout << "**  load                                    **" << endl;
    cout << "**********************************************" << endl;

    int loadId = myStructure.ConstraintLinearSetDisplacementNodeGroup(groupNodesTop, directionX, 0);

    cout << "**********************************************" << endl;
    cout << "**  visualization                           **" << endl;
    cout << "**********************************************" << endl;

    myStructure.AddVisualizationComponent(groupId, NuTo::VisualizeBase::DISPLACEMENTS);
    myStructure.AddVisualizationComponent(groupId, NuTo::VisualizeBase::CRACK_PHASE_FIELD);

    cout << "**********************************************" << endl;
    cout << "**  solver                                  **" << endl;
    cout << "**********************************************" << endl;


    myStructure.NodeBuildGlobalDofs();
    myStructure.CalculateMaximumIndependentSets();

    myIntegrationScheme.AddResultGroupNodeForce("myforce", groupNodesTop);

    center[0] = 0;
    center[1] = 1;
    int grpNodes_output_disp = myStructure.GroupCreate(NuTo::Groups::Nodes);
    myStructure.GroupAddNodeRadiusRange(grpNodes_output_disp, center, 0, tol);
    myIntegrationScheme.AddResultNodeDisplacements("mydisplacements",
                                                   myStructure.GroupGetMemberIds(grpNodes_output_disp).GetValue(0, 0));

    NuTo::FullMatrix<double, 2, 2> dispRHS;
    dispRHS(0, 0) = 0;
    dispRHS(1, 0) = simulationTime;
    dispRHS(0, 1) = 0;
    dispRHS(1, 1) = loadFactor;

    myIntegrationScheme.AddTimeDependentConstraint(loadId, dispRHS);

    try
    {
        myIntegrationScheme.Solve(simulationTime);
    }
    catch (...)
    {
        cout << "!!! SOMETHING WENT WRONG !!!" << endl;
        myStructure.ExportVtkDataFileElements(resultPath.string() + "bla.vtu", true);
    }


    std::string command = "paste " + resultPath.string() + "myforce.dat " + resultPath.string() +
                          "mydisplacements.dat > " + resultPath.string() + "forceDisp.dat";
    system(command.c_str());
    std::string newcommand = "python " + resultPath.parent_path().parent_path().string() + "/paraviewPythonScript.py " +
                             resultPath.string() + "/Group999_ElementsAll.pvd; okular /home/phuschke/test.png";
    system(newcommand.c_str());


    cout << "**********************************************" << endl;
    cout << "**  end                                     **" << endl;
    cout << "**********************************************" << endl;
}
