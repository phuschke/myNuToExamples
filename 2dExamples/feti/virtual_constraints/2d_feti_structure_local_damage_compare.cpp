

#include <iostream>
#include <vector>
#include "mechanics/feti/StructureFeti.h"
#include <ctime>
#include <chrono>
#include "mechanics/timeIntegration/NewmarkDirect.h"
#include "../../../EnumsAndTypedefs.h"

#include "mechanics/nodes/NodeBase.h"

#include "math/SparseMatrixCSR.h"
#include "math/SparseMatrixCSRGeneral.h"
#include "boost/filesystem.hpp"

constexpr int dim = 2;
using Eigen::VectorXd;
using Eigen::MatrixXd;


constexpr int dimension = 2;
constexpr double thickness = 1.0;

// material
constexpr double youngsModulus = 4.0e4;
constexpr double poissonsRatio = 0.2;
constexpr double tensileStrength = 3;
constexpr double compressiveStrength = 30;
constexpr double fractureEnergy = 0.1;

// integration
constexpr bool performLineSearch = true;
constexpr bool automaticTimeStepping = true;
constexpr double timeStep = 1e-0;
constexpr double minTimeStep = 1e-3;
constexpr double maxTimeStep = 1e-0;
constexpr double toleranceDisp = 1e-6;
constexpr double simulationTime = 1.0;
constexpr double loadFactor = -2e-2;
constexpr double maxInterations = 10;

const Eigen::Vector2d directionX = Eigen::Vector2d::UnitX();
const Eigen::Vector2d directionY = Eigen::Vector2d::UnitY();

void AssignSection(NuTo::Structure& structure);
void AssignMaterial(NuTo::Structure& structure);


int main(int argc, char* argv[])
{

    NuTo::Structure structure(dim);
    structure.SetNumTimeDerivatives(0);
    structure.SetVerboseLevel(10);
    structure.SetShowTime(false);
    structure.GetLogger().OpenFile("output_compare");

    std::string meshFile = argv[1];
    structure.GetLogger() << meshFile << "\n";

    auto bla = structure.ImportFromGmsh(meshFile);


    const int interpolationTypeId = bla[0].second;

    //    structure.InterpolationTypeAdd(interpolationTypeId, eDof::COORDINATES,     eTypeOrder::EQUIDISTANT1);
    structure.InterpolationTypeAdd(interpolationTypeId, eDof::DISPLACEMENTS, eTypeOrder::EQUIDISTANT1);
    structure.InterpolationTypeSetIntegrationType(interpolationTypeId,
                                                  NuTo::eIntegrationType::IntegrationType2D4NGauss4Ip);


    std::cout << " set interpolation type" << std::endl;
    //    structure.ElementTotalSetInterpolationType(interpolationTypeId);
    std::cout << " convert to interpolation type" << std::endl;
    structure.ElementTotalConvertToInterpolationType();

    structure.InterpolationTypeInfo(interpolationTypeId);
    structure.NodeInfo(10);

    AssignMaterial(structure);
    AssignSection(structure);


    structure.NodeBuildGlobalDofs(__PRETTY_FUNCTION__);

    structure.GetLogger() << "**********************************************"
                          << "\n";
    structure.GetLogger() << "**  real constraints                        **"
                          << "\n";
    structure.GetLogger() << "**********************************************"
                          << "\n\n";

    Eigen::VectorXd nodeCoords(2);
    nodeCoords[0] = 0;
    nodeCoords[1] = 0;

    int groupNodesLeftBoundary = structure.GroupCreate(eGroupId::Nodes);


    structure.GroupAddNodeRadiusRange(groupNodesLeftBoundary, nodeCoords, 0, 1.e-6);

    structure.ConstraintLinearSetDisplacementNodeGroup(groupNodesLeftBoundary, directionX, 0.0);
    structure.ConstraintLinearSetDisplacementNodeGroup(groupNodesLeftBoundary, directionY, 0.0);

    nodeCoords[0] = 60;
    nodeCoords[1] = 0;

    int groupNodesRightBoundary = structure.GroupCreate(eGroupId::Nodes);
    structure.GroupAddNodeRadiusRange(groupNodesRightBoundary, nodeCoords, 0, 1.e-6);

    structure.ConstraintLinearSetDisplacementNodeGroup(groupNodesRightBoundary, directionX, 0.0);
    structure.ConstraintLinearSetDisplacementNodeGroup(groupNodesRightBoundary, directionY, 0.0);


    structure.GetLogger() << "**********************************************"
                          << "\n";
    structure.GetLogger() << "**  load                                    **"
                          << "\n";
    structure.GetLogger() << "**********************************************"
                          << "\n\n";


    int loadNodeGroup = structure.GroupCreate(eGroupId::Nodes);

    nodeCoords[0] = 20;
    nodeCoords[1] = 10;
    structure.GroupAddNodeRadiusRange(loadNodeGroup, nodeCoords, 0, 1.e-6);
    int loadId = structure.ConstraintLinearSetDisplacementNodeGroup(loadNodeGroup, directionY, 1);


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
    //    structure.AddVisualizationComponent(groupAllElements, eVisualizeWhat::ENGINEERING_STRAIN);
    //    structure.AddVisualizationComponent(groupAllElements, eVisualizeWhat::ENGINEERING_STRESS);
    structure.AddVisualizationComponent(groupAllElements, eVisualizeWhat::DAMAGE);


    structure.GetLogger() << "**********************************************"
                          << "\n";
    structure.GetLogger() << "**  integration sheme                       **"
                          << "\n";
    structure.GetLogger() << "**********************************************"
                          << "\n\n";


    NuTo::NewmarkDirect myIntegrationScheme(&structure);
    boost::filesystem::path resultPath(std::string("/home/phuschke/results/feti/compare"));

    myIntegrationScheme.SetTimeStep(timeStep);
    myIntegrationScheme.SetMaxNumIterations(maxInterations);
    myIntegrationScheme.SetMinTimeStep(minTimeStep);
    myIntegrationScheme.SetMaxTimeStep(maxTimeStep);
    myIntegrationScheme.SetAutomaticTimeStepping(automaticTimeStepping);
    myIntegrationScheme.SetResultDirectory(resultPath.string(), true);
    myIntegrationScheme.SetPerformLineSearch(performLineSearch);
    myIntegrationScheme.SetToleranceResidual(eDof::DISPLACEMENTS, toleranceDisp);

    Eigen::Matrix2d dispRHS;
    dispRHS(0, 0) = 0;
    dispRHS(1, 0) = simulationTime;
    dispRHS(0, 1) = 0;
    dispRHS(1, 1) = loadFactor;

    myIntegrationScheme.AddResultGroupNodeForce("myforcebla", loadNodeGroup);

    nodeCoords[0] = 20;
    nodeCoords[1] = 10;
    int grpNodes_output_disp = structure.GroupCreate(eGroupId::Nodes);
    structure.GroupAddNodeRadiusRange(grpNodes_output_disp, nodeCoords, 0, 1.e-6);
    myIntegrationScheme.AddResultNodeDisplacements("mydisplacements",
                                                   structure.GroupGetMemberIds(grpNodes_output_disp)[0]);

    myIntegrationScheme.AddTimeDependentConstraint(loadId, dispRHS);
    //    myIntegrationScheme.SetTimeDependentLoadCase(loadId, dispRHS);

    structure.GetLogger() << "***********************************"
                          << "\n";
    structure.GetLogger() << "**      Solve                    **"
                          << "\n";
    structure.GetLogger() << "***********************************"
                          << "\n\n";

    structure.GetLogger() << "Total number of Dofs: \t" << structure.GetNumTotalDofs() << "\n\n";

    myIntegrationScheme.Solve(simulationTime);
}


void AssignSection(NuTo::Structure& structure)
{
    structure.GetLogger() << "***********************************"
                          << "\n";
    structure.GetLogger() << "**      Section                  **"
                          << "\n";
    structure.GetLogger() << "***********************************"
                          << "\n\n";

    int section00 = structure.SectionCreate(NuTo::eSectionType::PLANE_STRESS);
    structure.SectionSetThickness(section00, thickness);

    structure.ElementTotalSetSection(section00);
}

void AssignMaterial(NuTo::Structure& structure)
{
    structure.GetLogger() << "***********************************"
                          << "\n";
    structure.GetLogger() << "**      Material                 **"
                          << "\n";
    structure.GetLogger() << "***********************************"
                          << "\n\n";

    int material00 = structure.ConstitutiveLawCreate(eConstitutiveType::LOCAL_DAMAGE_MODEL);

    structure.ConstitutiveLawSetParameterDouble(material00, eConstitutiveParameter::YOUNGS_MODULUS, youngsModulus);
    structure.ConstitutiveLawSetParameterDouble(material00, eConstitutiveParameter::POISSONS_RATIO, poissonsRatio);
    structure.ConstitutiveLawSetParameterDouble(material00, eConstitutiveParameter::TENSILE_STRENGTH, tensileStrength);
    structure.ConstitutiveLawSetParameterDouble(material00, eConstitutiveParameter::FRACTURE_ENERGY, fractureEnergy);
    structure.ConstitutiveLawSetParameterDouble(material00, eConstitutiveParameter::COMPRESSIVE_STRENGTH,
                                                compressiveStrength);

    structure.ElementTotalSetConstitutiveLaw(material00);
}
