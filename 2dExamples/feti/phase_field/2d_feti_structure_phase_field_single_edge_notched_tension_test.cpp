
#include <mpi.h>
#include <boost/mpi.hpp>
#include "mechanics/constitutive/laws/PhaseField.h"
#include "mechanics/feti/StructureFeti.h"
#include <ctime>
#include <chrono>
#include "mechanics/feti/NewmarkFeti.h"
#include "../../../EnumsAndTypedefs.h"

#include "boost/filesystem.hpp"

#include "mechanics/sections/SectionPlane.h"

constexpr int dim = 2;
using Eigen::VectorXd;
using Eigen::MatrixXd;
using EigenSolver = Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>>;


constexpr int dimension = 2;
constexpr double thickness = 1.0;

constexpr double tol = 1e-6;
constexpr double toleranceDisp = 1e-8;
constexpr double toleranceCrack = 1e-8;

const Eigen::Vector2d directionX = Eigen::Vector2d::UnitX();
const Eigen::Vector2d directionY = Eigen::Vector2d::UnitY();

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

constexpr bool performLineSearch = true;
constexpr bool automaticTimeStepping = true;
constexpr double timeStep = 1.e-3;
constexpr double minTimeStep = 1.e-8;
constexpr double maxTimeStep = 1.e-2;
constexpr double timeStepPostProcessing = 1.e-5;

constexpr double simulationTime = 10.0e-3;
constexpr double loadFactor = simulationTime;


void AssignSection(NuTo::StructureFeti& structure);
void AssignMaterial(NuTo::StructureFeti& structure);


int main(int argc, char* argv[])
{
    boost::mpi::environment env(argc, argv);
    boost::mpi::communicator world;

    const int rank = world.rank();

    NuTo::StructureFeti structure(dim);
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
    AssignSection(structure);


    structure.GetLogger() << "*********************************** \n"
                          << "**      create node groups       ** \n"
                          << "*********************************** \n\n";

    Eigen::VectorXd nodeCoords(dim);

    int groupNodesBoundary = structure.GroupCreate(eGroupId::Nodes);
    structure.GroupAddNodeCoordinateRange(groupNodesBoundary, 1, -tol, tol);

    int groupNodesLoad = structure.GroupCreate(eGroupId::Nodes);
    structure.GroupAddNodeCoordinateRange(groupNodesLoad, 1, 1. - tol, 1. + tol);


    structure.GetLogger() << "*********************************** \n"
                          << "**      virtual constraints      ** \n"
                          << "*********************************** \n\n";


    std::vector<int> nodeIdsBoundaries = structure.GroupGetMemberIds(groupNodesBoundary);
    std::vector<int> nodeIdsLoads = structure.GroupGetMemberIds(groupNodesLoad);


    structure.ApplyVirtualConstraints(nodeIdsBoundaries, nodeIdsLoads);

    structure.GetLogger() << "*********************************** \n"
                          << "**      real constraints         ** \n"
                          << "*********************************** \n\n";


    structure.NodeBuildGlobalDofs(__PRETTY_FUNCTION__);

    structure.NodeInfo(10);

    int groupNodesBoundaryBottom = structure.GroupCreate(eGroupId::Nodes);
    structure.GroupAddNodeCoordinateRange(groupNodesBoundaryBottom, 1, -tol, tol);

    std::vector<int> boundaryDofIds;
    std::vector<int> nodeIdsBoundaryBottom = structure.GroupGetMemberIds(groupNodesBoundaryBottom);
    for (const int nodeId : nodeIdsBoundaryBottom)
    {
        std::vector<int> dofIds = structure.NodeGetDofIds(nodeId, eDof::DISPLACEMENTS);
        boundaryDofIds.push_back(dofIds[1]);
    }

    nodeCoords[0] = 0;
    nodeCoords[1] = 0;
    int groupNodesBoundaryBottomX = structure.GroupCreate(eGroupId::Nodes);
    structure.GroupAddNodeRadiusRange(groupNodesBoundaryBottomX, nodeCoords, 0, tol);


    std::vector<int> nodeIdsBoundaryBottomX = structure.GroupGetMemberIds(groupNodesBoundaryBottomX);
    for (const int nodeId : nodeIdsBoundaryBottomX)
    {
        std::vector<int> dofIds = structure.NodeGetDofIds(nodeId, eDof::DISPLACEMENTS);
        boundaryDofIds.push_back(dofIds[0]);
    }


    structure.ApplyConstraintsTotalFeti(boundaryDofIds);


    structure.GetLogger() << "*********************************** \n"
                          << "**      load                     ** \n"
                          << "*********************************** \n\n";


    std::map<int, double> dofIdAndPrescribedDisplacementMap;
    std::vector<int> nodeIds = structure.GroupGetMemberIds(groupNodesLoad);
    for (auto const& nodeId : nodeIds)
    {
        std::vector<int> dofIds = structure.NodeGetDofIds(nodeId, eDof::DISPLACEMENTS);
        dofIdAndPrescribedDisplacementMap.emplace(dofIds[1], 1.);
    }

    structure.ApplyPrescribedDisplacements(dofIdAndPrescribedDisplacementMap);

    int loadId = structure.LoadCreateNodeGroupForce(0, groupNodesLoad, directionY, 0);


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
    boost::filesystem::path resultPath(boost::filesystem::path(getenv("HOME")).string() +
                                       std::string("/results/feti/phase_field/") + std::to_string(structure.mRank));

    newmarkFeti.SetTimeStep(timeStep);

    newmarkFeti.SetMinTimeStep(minTimeStep);
    newmarkFeti.SetMaxTimeStep(maxTimeStep);
    newmarkFeti.SetMaxNumIterations(5);
    newmarkFeti.SetAutomaticTimeStepping(automaticTimeStepping);
    newmarkFeti.SetResultDirectory(resultPath.string(), true);
    newmarkFeti.SetPerformLineSearch(performLineSearch);
    newmarkFeti.SetToleranceResidual(eDof::DISPLACEMENTS, toleranceDisp);
    newmarkFeti.SetToleranceResidual(eDof::CRACKPHASEFIELD, toleranceCrack);
    newmarkFeti.SetIterativeSolver(NuTo::NewmarkFeti<EigenSolver>::eIterativeSolver::ProjectedGmres);
    newmarkFeti.SetFetiPreconditioner(NuTo::NewmarkFeti<EigenSolver>::eFetiPreconditioner::Lumped);
    newmarkFeti.SetMinTimeStepPlot(timeStepPostProcessing);


    Eigen::Matrix2d dispRHS;
    dispRHS(0, 0) = 0;
    dispRHS(1, 0) = simulationTime;
    dispRHS(0, 1) = 0;
    dispRHS(1, 1) = loadFactor;

    //    newmarkFeti.AddTimeDependentConstraint(loadId, dispRHS);
    newmarkFeti.SetTimeDependentLoadCase(loadId, dispRHS);

    structure.GetLogger() << "*********************************** \n"
                          << "**      solve                    ** \n"
                          << "*********************************** \n\n";


    //    if (rank == 0)
    //    {
    //        newmarkFeti.AddResultGroupNodeForce("force", groupNodesLoad);
    //        nodeCoords[0] = 0;
    //        nodeCoords[1] = 1;
    //        int grpNodes_output_disp = structure.GroupCreate(eGroupId::Nodes);
    //        structure.GroupAddNodeRadiusRange(grpNodes_output_disp, nodeCoords, 0, tol);
    //        newmarkFeti.AddResultNodeDisplacements("displacements",
    //                                                       structure.GroupGetMemberIds(grpNodes_output_disp)[0]);
    //    }
    //    else if (rank == 1)
    //    {
    //        newmarkFeti.AddResultGroupNodeForce("force", groupNodesLoad);
    //        nodeCoords[0] = 1;
    //        nodeCoords[1] = 1;
    //        int grpNodes_output_disp = structure.GroupCreate(eGroupId::Nodes);
    //        structure.GroupAddNodeRadiusRange(grpNodes_output_disp, nodeCoords, 0, tol);
    //        newmarkFeti.AddResultNodeDisplacements("displacements",
    //                                                       structure.GroupGetMemberIds(grpNodes_output_disp)[0]);
    //    }


    newmarkFeti.Solve(simulationTime);
    structure.GetLogger() << "Total number of Dofs: \t" << structure.GetNumTotalDofs() << "\n\n";
}


void AssignSection(NuTo::StructureFeti& structure)
{

    auto section = NuTo::SectionPlane::Create(thickness, true);
    structure.ElementTotalSetSection(section);
}

void AssignMaterial(NuTo::StructureFeti& structure)
{

    structure.GetLogger() << "*********************************** \n"
                          << "**      material                 ** \n"
                          << "*********************************** \n\n";

    NuTo::ConstitutiveBase* phaseField = new NuTo::PhaseField(youngsModulus, poissonsRatio, lengthScaleParameter,
                                                              fractureEnergy, artificialViscosity, energyDecomposition);

    int material00 = structure.AddConstitutiveLaw(phaseField);

    structure.ElementTotalSetConstitutiveLaw(material00);
}
