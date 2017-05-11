
#include "mechanics/structures/unstructured/Structure.h"
#include "mechanics/constitutive/laws/PhaseField.h"
#include "mechanics/timeIntegration/NewmarkDirect.h"
#include "mechanics/sections/SectionPlane.h"
#include "../../EnumsAndTypedefs.h"

#include <boost/filesystem.hpp>
#include <chrono>

using std::cout;
using std::endl;

constexpr int dimension = 2;
constexpr double thickness = 1.0;

// material
constexpr double youngsModulus = 4.0e4;
constexpr double poissonsRatio = 0.2;
constexpr double tensileStrength = 3;
constexpr double compressiveStrength = 30;
constexpr double lengthScaleParameter = 3.0e-2;
constexpr double fractureEnergy = 2.7;
constexpr double artificialViscosity = 0.5;

// integration
constexpr bool performLineSearch = true;
constexpr bool automaticTimeStepping = true;
constexpr double timeStep = 1e-3;
constexpr double minTimeStep = 1e-5;
constexpr double maxTimeStep = 1e-1;

constexpr double toleranceDisp = 1e-8;
constexpr double toleranceNlEqStrain = 1e-8;
constexpr double tolerance = 1e-5;

constexpr double simulationTime = 1.0;
constexpr double loadFactor = 0.5;
constexpr double maxIterations = 10;


const Eigen::Vector2d directionX = Eigen::Vector2d::UnitX();
const Eigen::Vector2d directionY = Eigen::Vector2d::UnitY();

constexpr ePhaseFieldEnergyDecomposition energyDecomposition = ePhaseFieldEnergyDecomposition::ISOTROPIC;


int main(int argc, char* argv[])
{

    cout << "**********************************************" << endl;
    cout << "**  strucutre                               **" << endl;
    cout << "**********************************************" << endl;

    NuTo::Structure structure(dimension);
    structure.SetNumTimeDerivatives(1);
    structure.SetShowTime(false);

    cout << "**********************************************" << endl;
    cout << "**  integration sheme                       **" << endl;
    cout << "**********************************************" << endl;

    boost::filesystem::path resultPath(boost::filesystem::path(getenv("HOME")).string() +
                                       std::string("/results/L_shaped_panel/"));

    NuTo::NewmarkDirect integrationScheme(&structure);
    integrationScheme.SetTimeStep(timeStep);
    integrationScheme.SetMinTimeStep(minTimeStep);
    integrationScheme.SetMaxTimeStep(maxTimeStep);
    integrationScheme.SetToleranceResidual(eDof::DISPLACEMENTS, toleranceDisp);
    integrationScheme.SetToleranceResidual(eDof::NONLOCALEQSTRAIN, toleranceNlEqStrain);
    integrationScheme.SetAutomaticTimeStepping(automaticTimeStepping);
    integrationScheme.SetPerformLineSearch(performLineSearch);
    integrationScheme.SetResultDirectory(resultPath.string(), true);

    structure.GetLogger() << "*********************************** \n"
                          << "**      geometry                 ** \n"
                          << "*********************************** \n\n";

    auto importContainer = structure.ImportFromGmsh(argv[1]);

    const int interpolationTypeId = importContainer[0].second;
    structure.InterpolationTypeAdd(interpolationTypeId, eDof::DISPLACEMENTS, eTypeOrder::EQUIDISTANT1);
    structure.InterpolationTypeAdd(interpolationTypeId, eDof::CRACKPHASEFIELD, eTypeOrder::EQUIDISTANT1);
    structure.ElementTotalConvertToInterpolationType();

    structure.ElementTotalConvertToInterpolationType(1.e-6, 10);

    structure.GetLogger() << "*********************************** \n"
                          << "**      section                  ** \n"
                          << "*********************************** \n\n";

    auto section = NuTo::SectionPlane::Create(thickness, true);
    structure.ElementTotalSetSection(section);

    structure.GetLogger() << "*********************************** \n"
                          << "**      material                 ** \n"
                          << "*********************************** \n\n";

    NuTo::ConstitutiveBase* phaseField = new NuTo::PhaseField(youngsModulus, poissonsRatio, lengthScaleParameter,
                                                              fractureEnergy, artificialViscosity, energyDecomposition);

    const int materialId = structure.AddConstitutiveLaw(phaseField);

    structure.ElementTotalSetConstitutiveLaw(materialId);

    cout << "**********************************************" << endl;
    cout << "**  bc                                      **" << endl;
    cout << "**********************************************" << endl;

    int groupNodesBoundary = structure.GroupCreate(eGroupId::Nodes);
    structure.GroupAddNodeCoordinateRange(groupNodesBoundary, 1, -1.e-6, +1.e-6);
    structure.ConstraintLinearSetDisplacementNodeGroup(groupNodesBoundary, Vector2d::UnitX(), 0.0);
    structure.ConstraintLinearSetDisplacementNodeGroup(groupNodesBoundary, Vector2d::UnitY(), 0.0);

    cout << "**********************************************" << endl;
    cout << "**  load                                    **" << endl;
    cout << "**********************************************" << endl;

    Eigen::VectorXd coordinates(dimension);
    coordinates << 470., 250.;

    int groupNodesLoad = structure.GroupCreate(NuTo::eGroupId::Nodes);
    structure.GroupAddNodeRadiusRange(groupNodesLoad, coordinates, 0, 1e-6);

    int load = structure.ConstraintLinearSetDisplacementNodeGroup(groupNodesLoad, Vector2d::UnitY(), 0);

    cout << "**********************************************" << endl;
    cout << "**  visualization                           **" << endl;
    cout << "**********************************************" << endl;

    int groupAllElements = 9999;
    structure.GroupCreate(groupAllElements, eGroupId::Elements);
    structure.GroupAddElementsTotal(groupAllElements);

    structure.AddVisualizationComponent(groupAllElements, NuTo::eVisualizeWhat::DISPLACEMENTS);
    structure.AddVisualizationComponent(groupAllElements, NuTo::eVisualizeWhat::ENGINEERING_STRAIN);
    structure.AddVisualizationComponent(groupAllElements, NuTo::eVisualizeWhat::ENGINEERING_STRESS);
    structure.AddVisualizationComponent(groupAllElements, NuTo::eVisualizeWhat::CRACK_PHASE_FIELD);

    cout << "**********************************************" << endl;
    cout << "**  solver                                  **" << endl;
    cout << "**********************************************" << endl;

    structure.NodeBuildGlobalDofs();

    Eigen::Matrix2d dispRHS;
    dispRHS(0, 0) = 0;
    dispRHS(1, 0) = simulationTime;
    dispRHS(0, 1) = 0;
    dispRHS(1, 1) = loadFactor;

    integrationScheme.AddTimeDependentConstraint(load, dispRHS);

    integrationScheme.Solve(simulationTime);

    cout << "**********************************************" << endl;
    cout << "**  end                                     **" << endl;
    cout << "**********************************************" << endl;
}
