#include "mechanics/structures/unstructured/Structure.h"
#include "mechanics/timeIntegration/NewmarkDirect.h"
#include "mechanics/constitutive/laws/PhaseField.h"
#include "mechanics/nodes/NodeEnum.h"
#include "mechanics/groups/GroupEnum.h"
#include "mechanics/constitutive/ConstitutiveEnum.h"
#include "visualize/VisualizeEnum.h"
#include "mechanics/interpolationtypes/InterpolationTypeEnum.h"
#include "mechanics/elements/IpDataEnum.h"
#include "mechanics/integrationtypes/IntegrationTypeEnum.h"
#include "mechanics/nodes/NodeBase.h"
#include "math/SparseMatrixCSR.h"
#include "mechanics/sections/SectionPlane.h"
#include "math/SparseMatrixCSRGeneral.h"
#include "mechanics/structures/StructureOutputBlockMatrix.h"
#include "../../../EnumsAndTypedefs.h"
#include <boost/filesystem.hpp>
#include <chrono>

using std::cout;
using std::endl;
using NuTo::Constitutive::ePhaseFieldEnergyDecomposition;

constexpr int dim = 2;
constexpr bool performLineSearch = true;
constexpr bool automaticTimeStepping = true;
constexpr double youngsModulus = 2.1e5;
constexpr double poissonsRatio = 0.3;
constexpr double thickness = 1.0;
constexpr double lengthScaleParameter = 3.0e-0;
constexpr double fractureEnergy = 2.7;
constexpr double artificialViscosity = 0.5;
constexpr double timeStep = 1e-2;
constexpr double minTimeStep = 1e-6;
constexpr double maxTimeStep = 1e-1;
constexpr double toleranceForce = 1e-6;
constexpr double simulationTime = 1.0;
constexpr double loadFactor = -0.4;
constexpr ePhaseFieldEnergyDecomposition energyDecomposition =
        ePhaseFieldEnergyDecomposition::ANISOTROPIC_SPECTRAL_DECOMPOSITION;


void AssignSection(NuTo::Structure& structure);
void AssignMaterial(NuTo::Structure& structure);


int main(int argc, char* argv[])
{


    NuTo::Structure structure(dim);
    structure.SetNumTimeDerivatives(1);
    structure.SetVerboseLevel(10);
    structure.SetShowTime(false);
    structure.GetLogger().OpenFile("output_compare");

    auto importContainer = structure.ImportFromGmsh(argv[1]);

    const int interpolationTypeId = importContainer[0].second;
    structure.InterpolationTypeAdd(interpolationTypeId, eDof::DISPLACEMENTS, eTypeOrder::EQUIDISTANT1);
    structure.InterpolationTypeAdd(interpolationTypeId, eDof::CRACKPHASEFIELD, eTypeOrder::EQUIDISTANT1);
    structure.ElementTotalConvertToInterpolationType();

    AssignMaterial(structure);
    AssignSection(structure);

    structure.GetLogger() << "*********************************** \n"
                          << "**      real constraints         ** \n"
                          << "*********************************** \n\n";

    Eigen::VectorXd nodeCoords(2);


    nodeCoords[0] = 0;
    nodeCoords[1] = 0;
    int groupNodesLeftBoundary = structure.GroupCreate(eGroupId::Nodes);
    structure.GroupAddNodeRadiusRange(groupNodesLeftBoundary, nodeCoords, 0, 1.e-6);
    structure.ConstraintLinearSetDisplacementNodeGroup(groupNodesLeftBoundary, Vector2d::UnitX(), 0.0);
    structure.ConstraintLinearSetDisplacementNodeGroup(groupNodesLeftBoundary, Vector2d::UnitY(), 0.0);

    nodeCoords[0] = 60;
    nodeCoords[1] = 0;
    int groupNodesRightBoundary = structure.GroupCreate(eGroupId::Nodes);
    structure.GroupAddNodeRadiusRange(groupNodesRightBoundary, nodeCoords, 0, 1.e-6);
    structure.ConstraintLinearSetDisplacementNodeGroup(groupNodesRightBoundary, Vector2d::UnitY(), 0.0);


    structure.GetLogger() << "**********************************************"
                          << "\n";
    structure.GetLogger() << "**  load                                    **"
                          << "\n";
    structure.GetLogger() << "**********************************************"
                          << "\n\n";


    int loadNodeGroup = structure.GroupCreate(eGroupId::Nodes);

    nodeCoords[0] = 30;
    nodeCoords[1] = 10;
    structure.GroupAddNodeRadiusRange(loadNodeGroup, nodeCoords, 0, 1.e-6);
    int loadId = structure.ConstraintLinearSetDisplacementNodeGroup(loadNodeGroup, Vector2d::UnitY(), 1);


    structure.GetLogger() << "***********************************"
                          << "\n";
    structure.GetLogger() << "**      Visualization            **"
                          << "\n";
    structure.GetLogger() << "***********************************"
                          << "\n\n";

    int groupAllElements = 9999;
    structure.GroupCreate(groupAllElements, eGroupId::Elements);
    structure.GroupAddElementsTotal(groupAllElements);
    structure.AddVisualizationComponent(groupAllElements, eVisualizeWhat::DISPLACEMENTS);
    structure.AddVisualizationComponent(groupAllElements, eVisualizeWhat::CRACK_PHASE_FIELD);
    //    structure.AddVisualizationComponent(groupAllElements, eVisualizeWhat::ENGINEERING_STRAIN);


    structure.GetLogger() << "**********************************************"
                          << "\n";
    structure.GetLogger() << "**  integration sheme                       **"
                          << "\n";
    structure.GetLogger() << "**********************************************"
                          << "\n\n";


    NuTo::NewmarkDirect myIntegrationScheme(&structure);
    boost::filesystem::path resultPath(std::string("/home/phuschke/results/feti/compare"));

    myIntegrationScheme.SetTimeStep(timeStep);
    myIntegrationScheme.SetMinTimeStep(minTimeStep);
    myIntegrationScheme.SetMaxTimeStep(maxTimeStep);
    myIntegrationScheme.SetAutomaticTimeStepping(automaticTimeStepping);
    myIntegrationScheme.SetResultDirectory(resultPath.string(), true);
    myIntegrationScheme.SetPerformLineSearch(performLineSearch);

    Eigen::Matrix2d dispRHS;
    dispRHS(0, 0) = 0;
    dispRHS(1, 0) = simulationTime;
    dispRHS(0, 1) = 0;
    dispRHS(1, 1) = loadFactor;

    myIntegrationScheme.AddResultGroupNodeForce("myforcebla", loadNodeGroup);

    nodeCoords[0] = 30;
    nodeCoords[1] = 10;
    int grpNodes_output_disp = structure.GroupCreate(eGroupId::Nodes);
    structure.GroupAddNodeRadiusRange(grpNodes_output_disp, nodeCoords, 0, 1.e-6);
    myIntegrationScheme.AddResultNodeDisplacements("mydisplacements",
                                                   structure.GroupGetMemberIds(grpNodes_output_disp)[0]);

    myIntegrationScheme.AddTimeDependentConstraint(loadId, dispRHS);
    //    myIntegrationScheme.SetTimeDependentLoadCase(loadId, dispRHS);


    structure.GetLogger() << "*********************************** \n"
                          << "**      solve                    ** \n"
                          << "*********************************** \n\n";

    myIntegrationScheme.Solve(simulationTime);

    structure.GetLogger() << "Total number of Dofs: \t" << structure.GetNumTotalDofs() << "\n\n";
}


void AssignSection(NuTo::Structure& structure)
{
    auto section = NuTo::SectionPlane::Create(thickness, true);
    structure.ElementTotalSetSection(section);
}

void AssignMaterial(NuTo::Structure& structure)
{

    NuTo::ConstitutiveBase* phaseField = new NuTo::PhaseField(youngsModulus, poissonsRatio, lengthScaleParameter,
                                                              fractureEnergy, artificialViscosity, energyDecomposition);

    int material00 = structure.AddConstitutiveLaw(phaseField);

    structure.ElementTotalSetConstitutiveLaw(material00);
}
