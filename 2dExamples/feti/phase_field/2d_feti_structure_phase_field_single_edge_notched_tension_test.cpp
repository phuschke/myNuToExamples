
#include <mpi.h>
#include <boost/mpi.hpp>
#include "mechanics/constitutive/laws/PhaseField.h"
#include "mechanics/feti/StructureFeti.h"
#include <ctime>
#include <chrono>
#include "mechanics/groups/Group.h"
#include "mechanics/feti/NewmarkFeti.h"
#include "../../../EnumsAndTypedefs.h"

#include "boost/filesystem.hpp"

#include "mechanics/sections/SectionPlane.h"
#include "mechanics/timeIntegration/postProcessing/PostProcessor.h"

using namespace NuTo;
using namespace Constitutive;
using namespace Interpolation;
using Node::eDof;

using Eigen::VectorXd;
using Eigen::Vector2d;
using Eigen::Matrix2d;
using EigenSolver = Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>>;
using FetiIterativeSolver = NewmarkFeti<EigenSolver>::eIterativeSolver;
using FetiScaling = NewmarkFeti<EigenSolver>::eFetiScaling;

// geometry
constexpr int dimension = 2;
constexpr double thickness = 1.0;
constexpr double lengthY = 1.0;
const Vector2d coordinateAtBottomLeft(0., 0.);

constexpr double tol = 1e-6;
constexpr double toleranceDisp = 1e-6;
constexpr double toleranceCrack = 1e-6;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//          phase-field model parameters
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// material
constexpr double youngsModulus = 2.1e5; // N/mm^2
constexpr double poissonsRatio = 0.3;
constexpr double lengthScaleParameter = 0.04; // mm
constexpr double fractureEnergy = 2.7; // N/mm
constexpr double artificialViscosity = 0.01; // Ns/mm^2
constexpr ePhaseFieldEnergyDecomposition energyDecomposition = ePhaseFieldEnergyDecomposition::ISOTROPIC;


constexpr bool performLineSearch = false;
constexpr bool automaticTimeStepping = true;
constexpr double timeStep = 1.e-4;
constexpr double minTimeStep = 1.e-8;
constexpr double maxTimeStep = 1.e-2;
constexpr double timeStepPostProcessing = 1.e-5;

constexpr double simulationTime = 10.0e-3;
constexpr double loadFactor = simulationTime;

void AssignMaterial(NuTo::StructureFeti& structure);


int main(int argc, char* argv[])
{
    boost::mpi::environment env(argc, argv);
    boost::mpi::communicator world;

    const int rank = world.rank();

    NuTo::StructureFeti structure(dimension);
    structure.SetNumTimeDerivatives(1);
    structure.SetVerboseLevel(10);
    structure.SetShowTime(true);
    structure.GetLogger().OpenFile("output" + std::to_string(rank));
    structure.GetLogger().SetQuiet(true);


    std::string meshFile = argv[1] + std::to_string(rank);
    structure.GetLogger() << meshFile << "\n";

    const int interpolationTypeId = structure.InterpolationTypeCreate(eShapeType::QUAD2D);
    structure.InterpolationTypeAdd(interpolationTypeId, eDof::COORDINATES, eTypeOrder::EQUIDISTANT1);
    structure.InterpolationTypeAdd(interpolationTypeId, eDof::DISPLACEMENTS, eTypeOrder::EQUIDISTANT1);
    structure.InterpolationTypeAdd(interpolationTypeId, eDof::CRACKPHASEFIELD, eTypeOrder::EQUIDISTANT1);

    structure.ImportMeshJson(meshFile, interpolationTypeId);


    AssignMaterial(structure);

    auto section = NuTo::SectionPlane::Create(thickness, true);
    structure.ElementTotalSetSection(section);

    structure.GetLogger() << "*********************************** \n"
                          << "**      create node groups       ** \n"
                          << "*********************************** \n\n";

    const auto& groupNodesBottomBoundary = structure.GroupGetNodesAtCoordinate(eDirection::Y, 0.);
    const auto& groupNodesLoad = structure.GroupGetNodesAtCoordinate(eDirection::Y, lengthY);

    structure.GetLogger() << "*********************************** \n"
                          << "**      virtual constraints      ** \n"
                          << "*********************************** \n\n";

    structure.ApplyVirtualConstraints(groupNodesBottomBoundary.GetMemberIds(), groupNodesLoad.GetMemberIds());

    structure.GetLogger() << "*********************************** \n"
                          << "**      real constraints         ** \n"
                          << "*********************************** \n\n";

    structure.NodeBuildGlobalDofs(__PRETTY_FUNCTION__);

    std::vector<int> boundaryDofIds;
    for (const int nodeId : groupNodesBottomBoundary.GetMemberIds())
    {
        std::vector<int> dofIds = structure.NodeGetDofIds(nodeId, eDof::DISPLACEMENTS);
        boundaryDofIds.push_back(dofIds[1]);
    }

    const auto& groupNodeAtBottomLeftCorner = structure.GroupGetNodeRadiusRange(coordinateAtBottomLeft);
    for (const int nodeId : groupNodeAtBottomLeftCorner.GetMemberIds())
    {
        std::vector<int> dofIds = structure.NodeGetDofIds(nodeId, eDof::DISPLACEMENTS);
        boundaryDofIds.push_back(dofIds[0]);
    }

    structure.ApplyConstraintsTotalFeti(boundaryDofIds);

    structure.GetLogger() << "*********************************** \n"
                          << "**      load                     ** \n"
                          << "*********************************** \n\n";

    std::map<int, double> dofIdAndPrescribedDisplacementMap;
    for (auto const& nodeId : groupNodesLoad.GetMemberIds())
    {
        std::vector<int> dofIds = structure.NodeGetDofIds(nodeId, eDof::DISPLACEMENTS);
        dofIdAndPrescribedDisplacementMap.emplace(dofIds[1], 1.);
    }

    structure.ApplyPrescribedDisplacements(dofIdAndPrescribedDisplacementMap);

    int loadId = structure.LoadCreateNodeGroupForce(&groupNodesLoad, Vector2d::UnitY(), 0.);

    structure.GetLogger() << "*********************************** \n"
                          << "**      visualization            ** \n"
                          << "*********************************** \n\n";

    int groupAllElements = 9999;
    structure.GroupCreate(groupAllElements, eGroupId::Elements);
    structure.GroupAddElementsTotal(groupAllElements);
    structure.AddVisualizationComponent(groupAllElements, eVisualizeWhat::DISPLACEMENTS);
    structure.AddVisualizationComponent(groupAllElements, eVisualizeWhat::CRACK_PHASE_FIELD);

    structure.GetLogger() << "*********************************** \n"
                          << "**      integration scheme       ** \n"
                          << "*********************************** \n\n";

    NuTo::NewmarkFeti<EigenSolver> newmarkFeti(&structure);
    boost::filesystem::path resultPath(boost::filesystem::initial_path().string() + "/results_" + std::to_string(rank));

    newmarkFeti.SetTimeStep(timeStep);
    newmarkFeti.SetMinTimeStep(minTimeStep);
    newmarkFeti.SetMaxTimeStep(maxTimeStep);
    newmarkFeti.SetMaxNumIterations(20);
    newmarkFeti.SetAutomaticTimeStepping(automaticTimeStepping);
    newmarkFeti.PostProcessing().SetResultDirectory(resultPath.string(), true);
    newmarkFeti.SetPerformLineSearch(performLineSearch);
    newmarkFeti.SetToleranceResidual(eDof::DISPLACEMENTS, toleranceDisp);
    newmarkFeti.SetToleranceResidual(eDof::CRACKPHASEFIELD, toleranceCrack);
    newmarkFeti.SetIterativeSolver(NuTo::NewmarkFeti<EigenSolver>::eIterativeSolver::ProjectedGmres);
    newmarkFeti.SetFetiPreconditioner(std::make_unique<NuTo::FetiLumpedPreconditioner>());
    newmarkFeti.SetMaxNumberOfFetiIterations(50);
    newmarkFeti.SetToleranceIterativeSolver(1.e-10);
    newmarkFeti.SetFetiScaling(FetiScaling::Multiplicity);
    newmarkFeti.PostProcessing().SetMinTimeStepPlot(timeStepPostProcessing);

    Matrix2d dispRHS;
    dispRHS(0, 0) = 0;
    dispRHS(1, 0) = simulationTime;
    dispRHS(0, 1) = 0;
    dispRHS(1, 1) = loadFactor;

    newmarkFeti.SetTimeDependentLoadCase(loadId, dispRHS);

    //    if (rank == 3)
    //    {
    //        int groupNodesLoadId = structure.GroupCreateNodeGroup();
    //        structure.GroupAddNodeCoordinateRange(groupNodesLoadId, eDirection::Y, 1.-tol, 1.+tol);
    //
    //        newmarkFeti.PostProcessing().AddResultGroupNodeForce("force", groupNodesLoadId);
    //        Vector2d coordinateAtTopLeft(1., 1.);
    //        int grpNodes_output_disp = structure.GroupCreate(eGroupId::Nodes);
    //        structure.GroupAddNodeRadiusRange(grpNodes_output_disp, coordinateAtTopLeft, 0, tol);
    //        newmarkFeti.PostProcessing().AddResultNodeDisplacements("displacements",
    //        structure.GroupGetMemberIds(grpNodes_output_disp)[0]);
    //    }
    //    else if (rank == 5)
    //    {
    //        int groupNodesLoadId = structure.GroupCreateNodeGroup();
    //        structure.GroupAddNodeCoordinateRange(groupNodesLoadId, eDirection::Y, 1.-tol, 1.+tol);
    //
    //        newmarkFeti.PostProcessing().AddResultGroupNodeForce("force", groupNodesLoadId);
    //        Vector2d coordinateAtTopRight(0., 1.);
    //        int grpNodes_output_disp = structure.GroupCreate(eGroupId::Nodes);
    //        structure.GroupAddNodeRadiusRange(grpNodes_output_disp, coordinateAtTopRight, 0, tol);
    //        newmarkFeti.PostProcessing().AddResultNodeDisplacements("displacements",
    //        structure.GroupGetMemberIds(grpNodes_output_disp)[0]);
    //    }

    structure.GetLogger() << "*********************************** \n"
                          << "**      solve                    ** \n"
                          << "*********************************** \n\n";

    newmarkFeti.Solve(simulationTime);
    structure.GetLogger() << "Total number of Dofs: \t" << structure.GetNumTotalDofs() << "\n\n";
}


void AssignMaterial(NuTo::StructureFeti& structure)
{

    structure.GetLogger() << "*********************************** \n"
                          << "**      material                 ** \n"
                          << "*********************************** \n\n";

    NuTo::ConstitutiveBase* phaseField = new NuTo::PhaseField(youngsModulus, poissonsRatio, lengthScaleParameter,
                                                              fractureEnergy, artificialViscosity, energyDecomposition);

    const int materialId = structure.AddConstitutiveLaw(phaseField);

    structure.ElementTotalSetConstitutiveLaw(materialId);
}
