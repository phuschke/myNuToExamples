#include "mechanics/structures/unstructured/Structure.h"
#include "mechanics/timeIntegration/NewmarkDirect.h"
#include "mechanics/constitutive/laws/PhaseField.h"
#include "mechanics/nodes/NodeEnum.h"
#include "mechanics/groups/GroupEnum.h"
#include "mechanics/sections/SectionPlane.h"

#include "mechanics/constitutive/ConstitutiveEnum.h"
#include "visualize/VisualizeEnum.h"
#include "mechanics/interpolationtypes/InterpolationTypeEnum.h"
#include "mechanics/elements/IpDataEnum.h"
#include "mechanics/integrationtypes/IntegrationTypeEnum.h"
#include "mechanics/nodes/NodeBase.h"
#include "math/SparseMatrixCSR.h"
#include "math/SparseMatrixCSRGeneral.h"
#include "mechanics/structures/StructureOutputBlockMatrix.h"
#include "../../EnumsAndTypedefs.h"
#include <boost/filesystem.hpp>
#include <chrono>

using std::cout;
using std::endl;
using NuTo::Constitutive::ePhaseFieldEnergyDecomposition;

constexpr   int         dimension                   = 2;
constexpr   bool        performLineSearch           = true;
constexpr   bool        automaticTimeStepping       = true;
constexpr   double      youngsModulus               = 2.1e5;
constexpr   double      poissonsRatio               = 0.3;
constexpr   double      thickness                   = 1.0;
constexpr   double      lengthScaleParameter        = 3.0e-2;
constexpr   double      fractureEnergy              = 2.7;
constexpr   double      artificialViscosity         = 0.5;
constexpr   double      timeStep                    = 1e-2;
constexpr   double      minTimeStep                 = 1e-6;
constexpr   double      maxTimeStep                 =  1e-1;
constexpr   double      toleranceForce              = 1e-6;
constexpr   double      simulationTime              = 1.0;
constexpr   double      loadFactor                  = -0.4;
constexpr   ePhaseFieldEnergyDecomposition energyDecomposition = ePhaseFieldEnergyDecomposition::ANISOTROPIC_SPECTRAL_DECOMPOSITION;

void WriteParameters(const int subdirectory, const int dispOrder, const int phaseFieldOrder, const int ipOrder, boost::filesystem::path resultPath)
{
    std::ofstream file;
    file.open(std::string(resultPath.string() + "parameters.txt"));
    cout << "Writing simulation parameters to file: " << std::string(resultPath.string() + "parameters.txt") << std::endl;
    file << "dispOrder            " << dispOrder            << std::endl;
    file << "phaseFieldOrder      " << phaseFieldOrder      << std::endl;
    file << "ipOrder              " << ipOrder              << std::endl;
    file << "subdirectory         " << subdirectory         << std::endl;
    file << "dimension            " << dimension            << std::endl;
    file << "performLineSearch    " << performLineSearch    << std::endl;
    file << "automaticTimeStepping" << automaticTimeStepping<< std::endl;
    file << "youngsModulus        " << youngsModulus        << std::endl;
    file << "poissonsRatio        " << poissonsRatio        << std::endl;
    file << "thickness            " << thickness            << std::endl;
    file << "lengthScaleParameter " << lengthScaleParameter << std::endl;
    file << "fractureEnergy       " << fractureEnergy       << std::endl;
    file << "artificialViscosity  " << artificialViscosity  << std::endl;
    file << "timeStep             " << timeStep             << std::endl;
    file << "minTimeStep          " << minTimeStep          << std::endl;
    file << "maxTimeStep          " << maxTimeStep          << std::endl;
    file << "toleranceForce       " << toleranceForce       << std::endl;
    file << "simulationTime       " << simulationTime       << std::endl;
    file << "loadFactor           " << loadFactor           << std::endl;

    file.close();
}

int main(int argc, char* argv[])
{

    if (argc != 5)
    {
        std::cout << "input arguments: displacement order, damage order, integration order, subdirectory" << std::endl;
        return EXIT_FAILURE;
    }

    const       int         dispOrder                   = std::stoi(argv[1]);
    const       int         phaseFieldOrder             = std::stoi(argv[2]);
    const       int         ipOrder                     = std::stoi(argv[3]);
    const       int         subdirectory                = std::stoi(argv[4]);


    boost::filesystem::path resultPath(std::string("/home/phuschke/results/2d/2d_miehe_three_point_bending_phase_field/" + std::to_string(subdirectory) + "/"));
    const boost::filesystem::path meshFilePath("2d_miehe_symmetric_three_point_bending_no_notch.msh");

    const Eigen::Vector2d directionX    = Eigen::Vector2d::UnitX();
    const Eigen::Vector2d directionY    = Eigen::Vector2d::UnitY();



    cout << "**********************************************" << endl;
    cout << "**  strucutre                               **" << endl;
    cout << "**********************************************" << endl;

    NuTo::Structure myStructure(dimension);
    myStructure.SetNumTimeDerivatives(1);
    myStructure.SetShowTime(false);

    cout << "Writing results to:" << endl;
    cout << resultPath.string() << endl;
    boost::filesystem::create_directory(resultPath);


    cout << "**********************************************" << endl;
    cout << "**  integration sheme                       **" << endl;
    cout << "**********************************************" << endl;

    NuTo::NewmarkDirect myIntegrationScheme(&myStructure);

    myIntegrationScheme.SetTimeStep                 ( timeStep                  );
    myIntegrationScheme.SetMinTimeStep              ( minTimeStep               );
    myIntegrationScheme.SetMaxTimeStep              ( maxTimeStep               );
    myIntegrationScheme.SetToleranceForce           ( toleranceForce            );
    myIntegrationScheme.SetAutomaticTimeStepping    ( automaticTimeStepping     );
    myIntegrationScheme.SetPerformLineSearch        ( performLineSearch         );
    myIntegrationScheme.SetResultDirectory          ( resultPath.string(), true );

    WriteParameters(subdirectory, dispOrder, phaseFieldOrder, ipOrder, resultPath);

    cout << "**********************************************" << endl;
    cout << "**  section                                 **" << endl;
    cout << "**********************************************" << endl;

    auto section = NuTo::SectionPlane::Create(thickness,true);
    myStructure.ElementTotalSetSection(section);

    cout << "**********************************************" << endl;
    cout << "**  material                                **" << endl;
    cout << "**********************************************" << endl;

    NuTo::ConstitutiveBase* phaseField = new NuTo::PhaseField(youngsModulus,
                                                              poissonsRatio,
                                                              lengthScaleParameter,
                                                              fractureEnergy,
                                                              artificialViscosity,
                                                              energyDecomposition
                                                              );

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
            myStructure.InterpolationTypeSetIntegrationType(myInterpolationType, NuTo::eIntegrationType::IntegrationType2D3NGauss1Ip);
            break;
        case 2:
            myStructure.InterpolationTypeSetIntegrationType(myInterpolationType, NuTo::eIntegrationType::IntegrationType2D3NGauss3Ip);
            break;
        case 3:
            myStructure.InterpolationTypeSetIntegrationType(myInterpolationType, NuTo::eIntegrationType::IntegrationType2D3NGauss6Ip);
            break;
        case 4:
            myStructure.InterpolationTypeSetIntegrationType(myInterpolationType, NuTo::eIntegrationType::IntegrationType2D3NGauss12Ip);
            break;
        default:
            std::cout << "ipOrder either 2, 3 or 4." << std::endl;
            return EXIT_FAILURE;
    }

    myStructure.InterpolationTypeInfo(myInterpolationType);

    myStructure.ElementTotalConvertToInterpolationType(1.e-6, 10);

    myStructure.ElementTotalSetSection(section);
    myStructure.ElementTotalSetConstitutiveLaw(matrixMaterial);

    cout << "**********************************************" << endl;
    cout << "**  bc                                      **" << endl;
    cout << "**********************************************" << endl;
    Eigen::Vector2d center;
    constexpr double tol = 1.0e-6;
    // bottom left boundary
    center[0] = 0.0;
    center[1] = 0.0;
    int grpNodes_bottom_left = myStructure.GroupCreate(NuTo::eGroupId::Nodes);
    myStructure.GroupAddNodeRadiusRange(grpNodes_bottom_left, center, 0, tol);
    myStructure.ConstraintLinearSetDisplacementNodeGroup(grpNodes_bottom_left,  directionY, 0);
    myStructure.ConstraintLinearSetDisplacementNodeGroup(grpNodes_bottom_left,  directionX, 0);

    // bottom right boundary
    center[0] = 8.0;
    center[1] = 0.0;

    int grpNodes_bottom_right = myStructure.GroupCreate(NuTo::eGroupId::Nodes);
    myStructure.GroupAddNodeRadiusRange(grpNodes_bottom_right, center, 0, tol);
    myStructure.ConstraintLinearSetDisplacementNodeGroup(grpNodes_bottom_right, directionY, 0);





    cout << "**********************************************" << endl;
    cout << "**  load                                    **" << endl;
    cout << "**********************************************" << endl;

    // middle top load
    center[0] = 4.0;
    center[1] = 2.0;
    int grpNodes_load = myStructure.GroupCreate(NuTo::eGroupId::Nodes);

    myStructure.GroupAddNodeRadiusRange(grpNodes_load, center, 0, 1e-6);


    int loadId = myStructure.ConstraintLinearSetDisplacementNodeGroup(grpNodes_load, directionY, 0);

    cout << "**********************************************" << endl;
    cout << "**  visualization                           **" << endl;
    cout << "**********************************************" << endl;

    myStructure.AddVisualizationComponent(groupId, NuTo::eVisualizeWhat::DISPLACEMENTS);
    myStructure.AddVisualizationComponent(groupId, NuTo::eVisualizeWhat::CRACK_PHASE_FIELD);

    cout << "**********************************************" << endl;
    cout << "**  solver                                  **" << endl;
    cout << "**********************************************" << endl;


    myStructure.NodeBuildGlobalDofs();
    myStructure.CalculateMaximumIndependentSets();


    center[0] = 4;
    center[1] = 2;
    int grpNodes_output = myStructure.GroupCreate(eGroupId::Nodes);
    myStructure.GroupAddNodeRadiusRange(grpNodes_output, center, 0, 5e-1);

    myIntegrationScheme.AddResultGroupNodeForce("myforce", grpNodes_output);

    center[0] = 4;
    center[1] = 2;
    int grpNodes_output_disp = myStructure.GroupCreate(NuTo::eGroupId::Nodes);
    myStructure.GroupAddNodeRadiusRange(grpNodes_output_disp, center, 0, 7e-1);
    myIntegrationScheme.AddResultNodeDisplacements("mydisplacements", myStructure.GroupGetMemberIds(grpNodes_output_disp)[0]);

    Eigen::Matrix2d dispRHS;
    dispRHS(0, 0) = 0;
    dispRHS(1, 0) =  simulationTime;
    dispRHS(0, 1) = 0;
    dispRHS(1, 1) =  loadFactor;

    myIntegrationScheme.AddTimeDependentConstraint(loadId, dispRHS);
    myStructure.ExportVtkDataFileElements(resultPath.string()+"bla.vtu", true);
    myIntegrationScheme.Solve(simulationTime);

//    std::string command = "paste " +  OutputPath.string() + "myforce.dat " +  OutputPath.string() + "mydisplacements.dat > " +  OutputPath.string() + "forceDisp.dat";
//            system(command.c_str());
    cout << "**********************************************" << endl;
    cout << "**  end                                     **" << endl;
    cout << "**********************************************" << endl;
}
