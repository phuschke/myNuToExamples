
#include "mechanics/structures/unstructured/Structure.h"
#include "mechanics/timeIntegration/NewmarkDirect.h"
#include "mechanics/constitutive/laws/PhaseField.h"
#include "mechanics/constraints/Constraints.h"
#include "mechanics/elements/ElementBase.h"
#include "mechanics/interpolationtypes/InterpolationType.h"
#include "mechanics/integrationtypes/IntegrationTypeBase.h"
#include "mechanics/constitutive/inputoutput/ConstitutiveCalculateStaticData.h"

#include "../../EnumsAndTypedefs.h"
#include <boost/filesystem.hpp>
#include <chrono>
#include <mechanics/groups/Group.h>

#include "mechanics/sections/SectionPlane.h"
#include "mechanics/constraints/ConstraintCompanion.h"
#include "mechanics/timeIntegration/postProcessing/PostProcessor.h"
using std::cout;
using std::endl;
using NuTo::Constitutive::ePhaseFieldEnergyDecomposition;
using NuTo::Constraint::Component;
using NuTo::Constraint::RhsRamp;

int main(int argc, char* argv[])
{

    // geometry
    constexpr int dimension = 2;
    constexpr double thickness = 1.0; // mm

    // material
    constexpr double youngsModulus = 2.1e5; // N/mm^2
    constexpr double poissonsRatio = 0.3;
    constexpr double lengthScaleParameter = 0.02; // mm
    constexpr double fractureEnergy = 2.7; // N/mm
    constexpr double artificialViscosity = 0.01; // Ns/mm^2
    constexpr ePhaseFieldEnergyDecomposition energyDecomposition = ePhaseFieldEnergyDecomposition::ISOTROPIC;

    constexpr bool performLineSearch = false;
    constexpr bool automaticTimeStepping = true;
    constexpr double timeStep = 1.e-3;
    constexpr double minTimeStep = 1.e-10;
    constexpr double maxTimeStep = 1.e-2;
    constexpr double timeStepPostProcessing = 5.e-5;

    constexpr double toleranceCrack = 1e-6;
    constexpr double toleranceDisp = 1e-6;
    constexpr double simulationTime = 20.0e-3;
    constexpr double loadFactor = simulationTime;

    constexpr double tol = 1.0e-8;
    boost::filesystem::path resultPath("results_wedge_split_test/");
    const boost::filesystem::path meshFilePath("wedge_split_test.msh");

    cout << "**********************************************" << endl;
    cout << "**  strucutre                               **" << endl;
    cout << "**********************************************" << endl;

    NuTo::Structure structure(dimension);
    structure.SetNumTimeDerivatives(1);
    structure.SetShowTime(false);

    cout << resultPath.string() << endl;
    boost::filesystem::create_directory(resultPath);


    cout << "**********************************************" << endl;
    cout << "**  integration scheme                       **" << endl;
    cout << "**********************************************" << endl;

    NuTo::NewmarkDirect integrationScheme(&structure);

    integrationScheme.SetTimeStep(timeStep);
    integrationScheme.SetMinTimeStep(minTimeStep);
    integrationScheme.SetMaxTimeStep(maxTimeStep);
    integrationScheme.SetAutomaticTimeStepping(automaticTimeStepping);
    integrationScheme.SetPerformLineSearch(performLineSearch);
    integrationScheme.PostProcessing().SetResultDirectory(resultPath.string(), true);
    integrationScheme.SetToleranceResidual(NuTo::Node::eDof::CRACKPHASEFIELD, toleranceCrack);
    integrationScheme.SetToleranceResidual(NuTo::Node::eDof::DISPLACEMENTS, toleranceDisp);
    // integrationScheme.SetMinTimeStepPlot(timeStepPostProcessing);

    std::ofstream file;
    file.open(std::string(resultPath.string() + "parameters.txt"));
    cout << "Writing simulation parameters to file: " << std::string(resultPath.string() + "parameters.txt")
         << std::endl;
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

    auto section = NuTo::SectionPlane::Create(thickness, true);
    structure.ElementTotalSetSection(section);

    cout << "**********************************************" << endl;
    cout << "**  material                                **" << endl;
    cout << "**********************************************" << endl;

    NuTo::ConstitutiveBase* phaseField = new NuTo::PhaseField(youngsModulus, poissonsRatio, lengthScaleParameter,
                                                              fractureEnergy, artificialViscosity, energyDecomposition);

    int matrixMaterial = structure.AddConstitutiveLaw(phaseField);

    cout << "**********************************************" << endl;
    cout << "**  geometry                                **" << endl;
    cout << "**********************************************" << endl;

    auto createdGroupIds = structure.ImportFromGmsh(meshFilePath.string());
    int groupId = createdGroupIds[0].first;

    const int interpolationTypeId = structure.InterpolationTypeCreate(eShapeType::QUAD2D);
    structure.InterpolationTypeAdd(interpolationTypeId, eDof::COORDINATES, eTypeOrder::EQUIDISTANT1);
    structure.InterpolationTypeAdd(interpolationTypeId, eDof::DISPLACEMENTS, eTypeOrder::EQUIDISTANT1);
    structure.InterpolationTypeAdd(interpolationTypeId, eDof::CRACKPHASEFIELD, eTypeOrder::EQUIDISTANT1);
    structure.ElementGroupSetInterpolationType(groupId, interpolationTypeId);

    structure.InterpolationTypeInfo(interpolationTypeId);
    structure.ElementTotalConvertToInterpolationType(1.e-9, 1);
    structure.ElementTotalSetSection(section);
    structure.ElementTotalSetConstitutiveLaw(matrixMaterial);

    cout << "**********************************************" << endl;
    cout << "**  bc                                      **" << endl;
    cout << "**********************************************" << endl;

    cout << "**********************************************" << endl;
    cout << "**  load                                    **" << endl;
    cout << "**********************************************" << endl;

    Vector2d coordinateTop({0.1, 0.3});
    auto& groupNodesLoad01 = structure.GroupGetNodeRadiusRange(coordinateTop, 0, 0.1);

    Vector2d coordinateBottom({0.1, -0.3});
    auto& groupNodesLoad02 = structure.GroupGetNodeRadiusRange(coordinateBottom, 0, 0.1);

    structure.Constraints().Add(eDof::DISPLACEMENTS, Component(groupNodesLoad01, {NuTo::eDirection::Y},
                                                               RhsRamp(simulationTime, loadFactor)));

    structure.Constraints().Add(eDof::DISPLACEMENTS, Component(groupNodesLoad02, {NuTo::eDirection::Y},
                                                               RhsRamp(simulationTime, -loadFactor)));

    cout << "**********************************************" << endl;
    cout << "**  visualization                           **" << endl;
    cout << "**********************************************" << endl;

    structure.AddVisualizationComponent(groupId, eVisualizeWhat::DISPLACEMENTS);
    structure.AddVisualizationComponent(groupId, eVisualizeWhat::CRACK_PHASE_FIELD);

    cout << "**********************************************" << endl;
    cout << "**  solver                                  **" << endl;
    cout << "**********************************************" << endl;

    structure.NodeBuildGlobalDofs();
    structure.CalculateMaximumIndependentSets();

    Eigen::MatrixXd dispRHS(dimension, dimension);
    dispRHS(0, 0) = 0;
    dispRHS(1, 0) = simulationTime;
    dispRHS(0, 1) = 0;
    dispRHS(1, 1) = loadFactor;

    int groupNodesLoadIdTop = structure.GroupCreateNodeGroup();
    structure.GroupAddNodeRadiusRange(groupNodesLoadIdTop,coordinateTop,0.,0.1);

    integrationScheme.PostProcessing().AddResultGroupNodeForce("forceTop", groupNodesLoadIdTop);
    integrationScheme.PostProcessing().AddResultNodeDisplacements("displacementsTop", structure.GroupGetMemberIds(groupNodesLoadIdTop)[0]);

    int groupNodesLoadIdBottom = structure.GroupCreateNodeGroup();
    structure.GroupAddNodeRadiusRange(groupNodesLoadIdBottom,coordinateBottom,0.,0.1);

    integrationScheme.PostProcessing().AddResultGroupNodeForce("forceBottom", groupNodesLoadIdBottom);
    integrationScheme.PostProcessing().AddResultNodeDisplacements("displacementsBottom", structure.GroupGetMemberIds(groupNodesLoadIdBottom)[0]);


    std::cout << "Number of elements: \t" << structure.GetNumElements() << std::endl;
    integrationScheme.Solve(simulationTime);


    cout << "**********************************************" << endl;
    cout << "**  end                                     **" << endl;
    cout << "**********************************************" << endl;
}
