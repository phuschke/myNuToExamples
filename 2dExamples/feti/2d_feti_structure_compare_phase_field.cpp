
#include <iostream>
#include <vector>
#include "mechanics/structures/unstructured/Structure.h"
#include "../ConjugateProjectedGradient.h"
#include <ctime>
#include <chrono>
#include "mechanics/timeIntegration/NewmarkDirect.h"
#include "../../EnumsAndTypedefs.h"
#include "mechanics/elements/ElementBase.h"
#include "mechanics/constitutive/laws/PhaseField.h"

#include "mechanics/nodes/NodeBase.h"
#include "mechanics/interpolationtypes/InterpolationType.h"
#include "mechanics/integrationtypes/IntegrationTypeBase.h"

#include "math/SparseMatrixCSR.h"
#include "math/SparseMatrixCSRGeneral.h"
#include "mechanics/structures/StructureOutputBlockMatrix.h"
#include "boost/filesystem.hpp"
#include "mechanics/sections/SectionPlane.h"


constexpr int dim = 2;
using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::cout;
using std::endl;


constexpr int dimension = 2;
constexpr double thickness = 1.0;

// material
constexpr double youngsModulus = 2.0; // N/mm^2
constexpr double poissonsRatio = 0.3;
constexpr double lengthScaleParameter = 1.0e-1; // mm
constexpr double fractureEnergy = 2.0e-4; // N/mm
constexpr double artificialViscosity = 0; // Ns/mm^2
constexpr ePhaseFieldEnergyDecomposition energyDecomposition = ePhaseFieldEnergyDecomposition::ISOTROPIC;

// integration
constexpr bool performLineSearch = false;
constexpr bool automaticTimeStepping = true;
constexpr double timeStep = 5e-2;
constexpr double minTimeStep = 1e-5;
constexpr double maxTimeStep = 5e-2;
constexpr double toleranceDisp = 1e-8;
constexpr double toleranceCrack = 1e-8;
constexpr double simulationTime = 1.0;
constexpr double loadFactor = 1;
constexpr double maxIterations = 10;
constexpr double tol = 1.e-6;
constexpr double lengthX = 60.0;


const Eigen::Vector2d directionX = Eigen::Vector2d::UnitX();
const Eigen::Vector2d directionY = Eigen::Vector2d::UnitY();

void AssignSection(NuTo::Structure& rStructure)
{
    auto section = NuTo::SectionPlane::Create(thickness, true);
    rStructure.ElementTotalSetSection(section);
}

void AssignMaterial(NuTo::Structure& rStructure)
{

    NuTo::ConstitutiveBase* phaseField = new NuTo::PhaseField(youngsModulus, poissonsRatio, lengthScaleParameter,
                                                              fractureEnergy, artificialViscosity, energyDecomposition);

    int material00 = rStructure.AddConstitutiveLaw(phaseField);

    rStructure.ElementTotalSetConstitutiveLaw(material00);
}


int main(int argc, char* argv[])
{


    std::string meshFile = argv[1];
    std::cout << meshFile << std::endl;

    NuTo::Structure structure(dim);
    structure.SetNumTimeDerivatives(1);
    structure.SetVerboseLevel(10);
    structure.SetShowTime(false);


    structure.ImportFromGmsh(meshFile);


    int interpolationType = structure.InterpolationTypeCreate(NuTo::Interpolation::eShapeType::QUAD2D);
    structure.InterpolationTypeAdd(interpolationType, NuTo::Node::eDof::COORDINATES,
                                   NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
    structure.InterpolationTypeAdd(interpolationType, NuTo::Node::eDof::DISPLACEMENTS,
                                   NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
    structure.InterpolationTypeAdd(interpolationType, NuTo::Node::eDof::CRACKPHASEFIELD,
                                   NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
    structure.ElementTotalSetInterpolationType(interpolationType);

    structure.ElementTotalConvertToInterpolationType();
    structure.NodeBuildGlobalDofs();
    AssignMaterial(structure);
    AssignSection(structure);

    cout << "**********************************************" << endl;
    cout << "**  constraints                             **" << endl;
    cout << "**********************************************" << endl;
    Eigen::VectorXd nodeCoords(2);

    int groupNodesLeftBoundary = structure.GroupCreate(eGroupId::Nodes);

    structure.GroupAddNodeCoordinateRange(groupNodesLeftBoundary, 0, -1.e-6, +1.e-6);

    //        nodeCoords[0] = 0;
    //        nodeCoords[1] = 0;
    //        structure.GroupAddNodeRadiusRange(groupNodesLeftBoundary, nodeCoords, 0, 1.e-6);

    structure.ConstraintLinearSetDisplacementNodeGroup(groupNodesLeftBoundary, directionX, 0);
    structure.ConstraintLinearSetDisplacementNodeGroup(groupNodesLeftBoundary, directionY, 0);

    //        int groupNodesRightBoundary = structure.GroupCreate(eGroupId::Nodes);
    ////            structure.GroupAddNodeCoordinateRange(groupNodesRightBoundary,0,60-1.e-6,60+1.e-6);

    //        nodeCoords[0] = 60;
    //        nodeCoords[1] = 0;
    //        structure.GroupAddNodeRadiusRange(groupNodesRightBoundary, nodeCoords, 0, 1.e-6);

    //        structure.ConstraintLinearSetDisplacementNodeGroup(groupNodesRightBoundary, directionX, 0);
    //        structure.ConstraintLinearSetDisplacementNodeGroup(groupNodesRightBoundary, directionY, 0);

    cout << "**********************************************" << endl;
    cout << "**  load                                    **" << endl;
    cout << "**********************************************" << endl;


    int loadNodeGroup = structure.GroupCreate(eGroupId::Nodes);

    //        nodeCoords[0] = 30;
    //        nodeCoords[1] = 10;
    //        structure.GroupAddNodeRadiusRange(loadNodeGroup, nodeCoords, 0, 1.e-6);

    structure.GroupAddNodeCoordinateRange(loadNodeGroup, 0, 60 - 1.e-6, 60 + 1.e-6);

    //    nodeCoords[0] = 29;
    //    nodeCoords[1] = 10;
    //    structure.GroupAddNodeRadiusRange(loadNodeGroup, nodeCoords, 0, 1.e-6);


    int loadId = structure.ConstraintLinearSetDisplacementNodeGroup(loadNodeGroup, directionX, 0);
    //    int loadId = structure.LoadCreateNodeGroupForce(0, loadNodeGroup, directionY, 1);


    cout << "**********************************************" << endl;
    cout << "** section weak                       **" << endl;
    cout << "**********************************************" << endl;

    auto sectionWeak = NuTo::SectionPlane::Create(thickness, true);
    structure.ElementTotalSetSection(sectionWeak);


    int groupNodeSectionWeak = structure.GroupCreate(eGroupId::Nodes);
    structure.GroupAddNodeCoordinateRange(groupNodeSectionWeak, 0, 30 - 2.1, 30 + 2.1);

    int groupEleSectionWeak = structure.GroupCreate(eGroupId::Elements);
    structure.GroupAddElementsFromNodes(groupEleSectionWeak, groupNodeSectionWeak, true);

    structure.ElementGroupSetSection(groupEleSectionWeak, sectionWeak);


    std::cout << "***********************************" << std::endl;
    std::cout << "**      Visualization rank = " << std::endl;
    std::cout << "***********************************" << std::endl;
    int groupAllElements = 9999;
    structure.GroupCreate(groupAllElements, eGroupId::Elements);
    structure.GroupAddElementsTotal(groupAllElements);
    structure.AddVisualizationComponent(groupAllElements, eVisualizeWhat::DISPLACEMENTS);
    structure.AddVisualizationComponent(groupAllElements, eVisualizeWhat::CRACK_PHASE_FIELD);

    cout << "**********************************************" << endl;
    cout << "**  integration sheme                       **" << endl;
    cout << "**********************************************" << endl;


    NuTo::NewmarkDirect myIntegrationScheme(&structure);
    boost::filesystem::path resultPath(std::string("/home/phuschke/results/feti/compare"));

    myIntegrationScheme.SetTimeStep(timeStep);
    myIntegrationScheme.SetMaxNumIterations(maxIterations);
    myIntegrationScheme.SetMinTimeStep(minTimeStep);
    myIntegrationScheme.SetMaxTimeStep(maxTimeStep);
    myIntegrationScheme.SetAutomaticTimeStepping(automaticTimeStepping);
    myIntegrationScheme.SetResultDirectory(resultPath.string(), true);
    myIntegrationScheme.SetPerformLineSearch(performLineSearch);
    myIntegrationScheme.SetToleranceResidual(eDof::DISPLACEMENTS, toleranceDisp);
    myIntegrationScheme.SetToleranceResidual(eDof::CRACKPHASEFIELD, toleranceCrack);


    //    myIntegrationScheme.AddResultGroupNodeForce("myforce", loadNodeGroup);
    //    nodeCoords[0]            = 60;
    //    nodeCoords[1]            = 10;
    //    int grpNodes_output_disp = structure.GroupCreate(eGroupId::Nodes);
    //    structure.GroupAddNodeRadiusRange(grpNodes_output_disp, nodeCoords, 0, 1.e-6);
    //    myIntegrationScheme.AddResultNodeDisplacements("mydisplacements",
    //                                                   structure.GroupGetMemberIds(grpNodes_output_disp)[0]);


    Eigen::Matrix2d dispRHS;
    dispRHS(0, 0) = 0;
    dispRHS(1, 0) = simulationTime;
    dispRHS(0, 1) = 0;
    dispRHS(1, 1) = loadFactor;

    myIntegrationScheme.AddTimeDependentConstraint(loadId, dispRHS);
    //    myIntegrationScheme.SetTimeDependentLoadCase(loadId, dispRHS);

    cout << "***********************************" << std::endl;
    cout << "**      Solve                    **" << std::endl;
    cout << "***********************************" << std::endl;

    myIntegrationScheme.Solve(simulationTime);


    std::cout << "END" << std::endl;
}
