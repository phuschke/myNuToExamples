
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

#include "mechanics/sections/SectionPlane.h"

using std::cout;
using std::endl;
using NuTo::Constitutive::ePhaseFieldEnergyDecomposition;


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

    // geometry
    constexpr   int         dimension                   = 2;
    constexpr   double      thickness                   = 1.0;                  // mm

    // material
    constexpr   double      youngsModulus               = 2.1e5;                // N/mm^2
    constexpr   double      poissonsRatio               = 0.3;
    constexpr   double      lengthScaleParameter        = 7.5e-0;               // mm
    constexpr   double      fractureEnergy              = 2.7;                  // N/mm
    constexpr   double      artificialViscosity         = 2.0e-0;               // Ns/mm^2
    constexpr   ePhaseFieldEnergyDecomposition energyDecomposition = ePhaseFieldEnergyDecomposition::ISOTROPIC;

    constexpr   bool        performLineSearch           = true;
    constexpr   bool        automaticTimeStepping       = true;
    constexpr   double      timeStep                   = 1.e-2;
    constexpr   double      minTimeStep                = 1.e-6;
    constexpr   double      maxTimeStep                = 1.e-2;
    constexpr   double      timeStepPostProcessing     = 5.e-2;

    constexpr   double      toleranceCrack             = 1e-4;
    constexpr   double      toleranceDisp              = 1e-5;
    constexpr   double      simulationTime             = 8.0e-1;
    constexpr   double      loadFactor                 = simulationTime;

    constexpr   double      tol                        = 1.0e-8;
    boost::filesystem::path resultPath(std::string("/home/phuschke/results/2d/2d_tension_test/" + std::to_string(subdirectory) + "/"));
    const boost::filesystem::path meshFilePath("2d_matrix_thought_experiment.msh");

    const Eigen::Vector2d directionX    = Eigen::Vector2d::UnitX();
    const Eigen::Vector2d directionY    = Eigen::Vector2d::UnitY();

    cout << "**********************************************" << endl;
    cout << "**  strucutre                               **" << endl;
    cout << "**********************************************" << endl;

    NuTo::Structure myStructure(dimension);
    myStructure.SetNumTimeDerivatives(1);
    myStructure.SetShowTime(false);

    cout << resultPath.string() << endl;
    boost::filesystem::create_directory(resultPath);

    cout << "**********************************************" << endl;
    cout << "**  integration scheme                       **" << endl;
    cout << "**********************************************" << endl;

    NuTo::NewmarkDirect myIntegrationScheme(&myStructure);

    myIntegrationScheme.SetTimeStep                 ( timeStep                  );
    myIntegrationScheme.SetMinTimeStep              ( minTimeStep               );
    myIntegrationScheme.SetMaxTimeStep              ( maxTimeStep               );
    myIntegrationScheme.SetAutomaticTimeStepping    ( automaticTimeStepping     );
    myIntegrationScheme.SetPerformLineSearch        ( performLineSearch         );
    myIntegrationScheme.SetResultDirectory          ( resultPath.string(), true );
//    myIntegrationScheme.SetToleranceResidual        ( NuTo::Node::eDof::CRACKPHASEFIELD, toleranceCrack  );
//    myIntegrationScheme.SetToleranceResidual        ( NuTo::Node::eDof::DISPLACEMENTS  , toleranceDisp   );
    myIntegrationScheme.SetExportDataFileNodes      ( false );
    myIntegrationScheme.SetMinTimeStepPlot          ( timeStepPostProcessing );

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
    file << "toleranceDisp        " << toleranceDisp        << std::endl;
    file << "toleranceCrack       " << toleranceCrack       << std::endl;
    file << "simulationTime       " << simulationTime       << std::endl;
    file << "loadFactor           " << loadFactor           << std::endl;

    file.close();

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

    myStructure.InterpolationTypeInfo(myInterpolationType);

    myStructure.ElementTotalConvertToInterpolationType(1.e-9, 1);

    myStructure.ElementTotalSetSection(section);
    myStructure.ElementTotalSetConstitutiveLaw(matrixMaterial);


    cout << "**********************************************" << endl;
    cout << "**  bc                                      **" << endl;
    cout << "**********************************************" << endl;

    Eigen::Vector2d center;

    int groupNodesBCLeft = myStructure.GroupCreate(eGroupId::Nodes);
    myStructure.GroupAddNodeCoordinateRange(groupNodesBCLeft, 0, - tol, tol);

    myStructure.ConstraintLinearSetDisplacementNodeGroup(groupNodesBCLeft, directionX, 0);


    cout << "**********************************************" << endl;
    cout << "**  load                                    **" << endl;
    cout << "**********************************************" << endl;

    int groupNodesLoadRight = myStructure.GroupCreate(eGroupId::Nodes);
    myStructure.GroupAddNodeCoordinateRange(groupNodesLoadRight, 0, 120-tol, 120+tol);

    int loadID = myStructure.ConstraintLinearSetDisplacementNodeGroup(groupNodesLoadRight,directionX, 0);

    cout << "**********************************************" << endl;
    cout << "**  visualization                           **" << endl;
    cout << "**********************************************" << endl;

    myStructure.AddVisualizationComponent(groupId, eVisualizeWhat::DISPLACEMENTS);
    myStructure.AddVisualizationComponent(groupId, eVisualizeWhat::CRACK_PHASE_FIELD);

    cout << "**********************************************" << endl;
    cout << "**  solver                                  **" << endl;
    cout << "**********************************************" << endl;

    myStructure.NodeBuildGlobalDofs();
    myStructure.CalculateMaximumIndependentSets();

    Eigen::MatrixXd dispRHS(dimension, dimension);
    dispRHS(0, 0) = 0;
    dispRHS(1, 0) =  simulationTime;
    dispRHS(0, 1) = 0;
    dispRHS(1, 1) =  loadFactor;

    myIntegrationScheme.AddTimeDependentConstraint(loadID, dispRHS);

    try
    {
        myIntegrationScheme.Solve(simulationTime);
    }
    catch(...)
    {
        cout << "!!! SOMETHING WENT WRONG !!!" << endl;
        myStructure.ExportVtkDataFileElements(resultPath.string()+"bla.vtu", true);
    }



    cout << "**********************************************" << endl;
    cout << "**  end                                     **" << endl;
    cout << "**********************************************" << endl;
}
















