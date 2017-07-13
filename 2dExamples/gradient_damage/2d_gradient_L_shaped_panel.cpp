
#include "mechanics/structures/unstructured/Structure.h"
#include "mechanics/constitutive/damageLaws/DamageLawExponential.h"
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
constexpr double nonlocalRadius = 4; // mm
constexpr double youngsModulus = 4.0e4;
constexpr double poissonsRatio = 0.2;
constexpr double tensileStrength = 3;
constexpr double compressiveStrength = 30;
constexpr double fractureEnergy = 0.01;
constexpr double alpha = 0.99;

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


int main(int argc, char* argv[])
{

    cout << "**********************************************" << endl;
    cout << "**  strucutre                               **" << endl;
    cout << "**********************************************" << endl;

    NuTo::Structure structure(dimension);
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
    structure.InterpolationTypeAdd(interpolationTypeId, eDof::NONLOCALEQSTRAIN, eTypeOrder::EQUIDISTANT1);
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

    int materialId = structure.ConstitutiveLawCreate(eConstitutiveType::GRADIENT_DAMAGE_ENGINEERING_STRESS);
    structure.ConstitutiveLawSetParameterDouble(materialId, eConstitutiveParameter::YOUNGS_MODULUS, youngsModulus);
    structure.ConstitutiveLawSetParameterDouble(materialId, eConstitutiveParameter::POISSONS_RATIO, poissonsRatio);
    structure.ConstitutiveLawSetParameterDouble(materialId, eConstitutiveParameter::TENSILE_STRENGTH, tensileStrength);
    structure.ConstitutiveLawSetParameterDouble(materialId, eConstitutiveParameter::COMPRESSIVE_STRENGTH,
                                                compressiveStrength);
    structure.ConstitutiveLawSetParameterDouble(materialId, eConstitutiveParameter::NONLOCAL_RADIUS, nonlocalRadius);

    structure.ConstitutiveLawSetDamageLaw(
            materialId, NuTo::Constitutive::DamageLawExponential::Create(tensileStrength / youngsModulus,
                                                                         tensileStrength / fractureEnergy, alpha));

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
    structure.AddVisualizationComponent(groupAllElements, NuTo::eVisualizeWhat::DAMAGE);

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
