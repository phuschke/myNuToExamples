
#include <mpi.h>
#include <boost/mpi.hpp>

#include "mechanics/feti/StructureFeti.h"
#include <ctime>
#include <chrono>
#include "mechanics/feti/NewmarkFeti.h"
#include "../../../EnumsAndTypedefs.h"

#include "mechanics/nodes/NodeBase.h"

#include "boost/filesystem.hpp"
#include "mechanics/sections/SectionPlane.h"

#include "mechanics/mesh/MeshGenerator.h"
constexpr int dim = 2;
using Eigen::VectorXd;
using Eigen::MatrixXd;

using EigenSolver = Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>>;

// geometry
constexpr double domainLength = 20.;
constexpr double thickness = 1.0;

// material
double youngsModulus = 4094;
constexpr double poissonsRatio = 0.3;


// integration
constexpr bool performLineSearch = false;
constexpr bool automaticTimeStepping = false;
constexpr double timeStep = 1e-0;
constexpr double minTimeStep = 1e-1;
constexpr double maxTimeStep = 1e-0;
constexpr double toleranceDisp = 1e-5;
constexpr double simulationTime = 1.0;
constexpr double loadFactor = 1.0;
constexpr double maxIterations = 1;
constexpr double tol = 1.e-6;

const Eigen::Vector2d directionX = Eigen::Vector2d::UnitX();
const Eigen::Vector2d directionY = Eigen::Vector2d::UnitY();

int main(int argc, char* argv[])
{
    boost::mpi::environment env(argc, argv);
    boost::mpi::communicator world;

    const int rank = world.rank();

    NuTo::StructureFeti structure(dim);
    structure.SetNumTimeDerivatives(0);
    structure.SetVerboseLevel(5);
    structure.SetShowTime(false);
    structure.GetLogger().OpenFile("output" + std::to_string(rank));
    structure.GetLogger().SetQuiet(true);

    std::string meshFile = argv[1] + std::to_string(rank);
    structure.GetLogger() << meshFile << "\n";

    const int interpolationTypeId = structure.InterpolationTypeCreate(eShapeType::QUAD2D);
    structure.InterpolationTypeAdd(interpolationTypeId, eDof::COORDINATES, eTypeOrder::EQUIDISTANT1);
    structure.InterpolationTypeAdd(interpolationTypeId, eDof::DISPLACEMENTS, eTypeOrder::EQUIDISTANT1);

    structure.ImportMeshJson(meshFile, interpolationTypeId);

    if (rank == 8 or rank == 9 or rank == 10 or rank == 11)
    {
        youngsModulus = 1;
    }


    const int materialId = structure.ConstitutiveLawCreate(eConstitutiveType::LINEAR_ELASTIC_ENGINEERING_STRESS);
    structure.ConstitutiveLawSetParameterDouble(materialId, eConstitutiveParameter::YOUNGS_MODULUS, youngsModulus);
    structure.ConstitutiveLawSetParameterDouble(materialId, eConstitutiveParameter::POISSONS_RATIO, poissonsRatio);
    structure.ElementTotalSetConstitutiveLaw(materialId);


    auto section = NuTo::SectionPlane::Create(thickness, true);
    structure.ElementTotalSetSection(section);

    structure.ExportVtkDataFileElements("elements_" + std::to_string(rank) + ".vtk");

    structure.GetLogger() << "*********************************** \n"
                          << "**      create node groups       ** \n"
                          << "*********************************** \n\n";

    Eigen::VectorXd nodeCoords(dim);

    int groupNodesBoundary = structure.GroupCreate(eGroupId::Nodes);
    structure.GroupAddNodeCoordinateRange(groupNodesBoundary, 0, -tol, tol);

    int groupNodesLoad = structure.GroupCreate(eGroupId::Nodes);
    structure.GroupAddNodeCoordinateRange(groupNodesLoad, 0, domainLength - tol, domainLength + tol);


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

    std::vector<int> boundaryDofIds;
    for (const int nodeId : nodeIdsBoundaries)
    {
        std::vector<int> dofIds = structure.NodeGetDofIds(nodeId, eDof::DISPLACEMENTS);
        boundaryDofIds.push_back(dofIds[0]);
        boundaryDofIds.push_back(dofIds[1]);
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
        dofIdAndPrescribedDisplacementMap.emplace(dofIds[0], 0.);
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
    structure.AddVisualizationComponent(groupAllElements, eVisualizeWhat::ENGINEERING_STRESS);


    structure.GetLogger() << "*********************************** \n"
                          << "**      integration scheme       ** \n"
                          << "*********************************** \n\n";

    NuTo::NewmarkFeti<EigenSolver> newmarkFeti(&structure);
    boost::filesystem::path resultPath(boost::filesystem::path(getenv("HOME")).string() +
                                       std::string("/results/feti/linear_elastic/rixen_heterogeneous_example/") +
                                       std::to_string(structure.mRank));

    newmarkFeti.SetTimeStep(timeStep);

    newmarkFeti.SetMinTimeStep(minTimeStep);
    newmarkFeti.SetMaxTimeStep(maxTimeStep);
    newmarkFeti.SetMaxNumIterations(maxIterations);
    newmarkFeti.SetAutomaticTimeStepping(automaticTimeStepping);
    newmarkFeti.SetResultDirectory(resultPath.string(), true);
    newmarkFeti.SetPerformLineSearch(performLineSearch);
    newmarkFeti.SetToleranceResidual(eDof::DISPLACEMENTS, toleranceDisp);
    newmarkFeti.SetIterativeSolver(NuTo::NewmarkFeti<EigenSolver>::eIterativeSolver::ConjugateGradient);
    newmarkFeti.SetFetiPreconditioner(NuTo::NewmarkFeti<EigenSolver>::eFetiPreconditioner::Lumped);
    newmarkFeti.SetMaxNumberOfFetiIterations(1000);
    newmarkFeti.SetToleranceIterativeSolver(1.e-6);

    Eigen::Matrix2d dispRHS;
    dispRHS(0, 0) = 0;
    dispRHS(1, 0) = simulationTime;
    dispRHS(0, 1) = 0;
    dispRHS(1, 1) = loadFactor;

    newmarkFeti.SetTimeDependentLoadCase(loadId, dispRHS);

    structure.GetLogger() << "*********************************** \n"
                          << "**      solve                    ** \n"
                          << "*********************************** \n\n";


    newmarkFeti.Solve(simulationTime);
    structure.GetLogger() << "Total number of Dofs: \t" << structure.GetNumTotalDofs() << "\n\n";
}
