
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
// using EigenSolver = Eigen::SparseQR<Eigen::SparseMatrix<double>,Eigen::COLAMDOrdering<int>>;
// using EigenSolver = Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>>;
// using EigenSolver = Eigen::PardisoLU<Eigen::SparseMatrix<double>>;
using EigenSolver = Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>>;

constexpr double thickness = 1.0;

// material
constexpr double youngsModulus = 2.1e5;
constexpr double poissonsRatio = 0.3;


// integration
constexpr bool performLineSearch = false;
constexpr bool automaticTimeStepping = false;
constexpr double timeStep = 1e-0;
constexpr double minTimeStep = 1e-1;
constexpr double maxTimeStep = 1e-0;
constexpr double toleranceDisp = 1e-6;
constexpr double simulationTime = 1.0;
constexpr double loadFactor = -10.0;
constexpr double maxIterations = 1;

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
    structure.SetShowTime(true);
    structure.GetLogger().OpenFile("output" + std::to_string(rank));
    structure.GetLogger().SetQuiet(true);

    std::vector<double> meshDimensions;
    meshDimensions.push_back(60.);
    meshDimensions.push_back(10.);

    std::vector<int> numElements;
    numElements.push_back(100);
    numElements.push_back(50);

    auto importContainer = structure.CreateRectangularMesh2D(meshDimensions, numElements);
    const int interpolationTypeId = importContainer.second;
    structure.InterpolationTypeAdd(interpolationTypeId, eDof::DISPLACEMENTS,   eTypeOrder::EQUIDISTANT1);
    structure.ElementTotalConvertToInterpolationType();

    const int materialId = structure.ConstitutiveLawCreate(eConstitutiveType::LINEAR_ELASTIC_ENGINEERING_STRESS);
    structure.ConstitutiveLawSetParameterDouble(materialId, eConstitutiveParameter::YOUNGS_MODULUS, youngsModulus);
    structure.ConstitutiveLawSetParameterDouble(materialId, eConstitutiveParameter::POISSONS_RATIO, poissonsRatio);
    structure.ElementTotalSetConstitutiveLaw(materialId);


    auto section = NuTo::SectionPlane::Create(thickness, true);
    structure.ElementTotalSetSection(section);


    structure.GetLogger() << "*********************************** \n"
                          << "**      node groups              ** \n"
                          << "*********************************** \n\n";

    Eigen::VectorXd coordinates(dim);

    int groupNodesLeftBoundary = structure.GroupCreate(eGroupId::Nodes);
    structure.GroupAddNodeCoordinateRange(groupNodesLeftBoundary, 0, -1.e-6, +1.e-6);

    int loadNodeGroup = structure.GroupCreate(eGroupId::Nodes);
    coordinates[0] = 60;
    coordinates[1] = 0;
    structure.GroupAddNodeRadiusRange(loadNodeGroup, coordinates, 0, 1.e-6);


    structure.GetLogger() << "*********************************** \n"
                          << "**      virtual constraints      ** \n"
                          << "*********************************** \n\n";

    std::vector<int> nodeIdsBoundaries = structure.GroupGetMemberIds(groupNodesLeftBoundary);
    std::vector<int> nodeIdsLoads = structure.GroupGetMemberIds(loadNodeGroup);

//    std::cout << nodeIdsBoundaries << "\n";

    structure.ApplyVirtualConstraints(nodeIdsBoundaries, nodeIdsLoads);

    structure.GetLogger() << "*********************************** \n"
                          << "**      real constraints         ** \n"
                          << "*********************************** \n\n";

    structure.NodeBuildGlobalDofs(__PRETTY_FUNCTION__);
    structure.ApplyConstraintsTotalFeti(groupNodesLeftBoundary);

    structure.GetLogger() << "*********************************** \n"
                          << "**      load                     ** \n"
                          << "*********************************** \n\n";


    // prescribe displacement of loadNodeGroup in Y direction
    std::map<int, double> dofIdAndPrescribedDisplacementMap;
    std::vector<int> nodeIds = structure.GroupGetMemberIds(loadNodeGroup);
    for (auto const& nodeId : nodeIds)
    {
        std::vector<int> dofIds = structure.NodeGetDofIds(nodeId, eDof::DISPLACEMENTS);
        dofIdAndPrescribedDisplacementMap.emplace(dofIds[1], 1.);
    }

    structure.ApplyPrescribedDisplacements(dofIdAndPrescribedDisplacementMap);

    int loadId = structure.LoadCreateNodeGroupForce(0, loadNodeGroup, directionY, 0);

    structure.GetLogger() << "*********************************** \n"
                          << "**      visualization            ** \n"
                          << "*********************************** \n\n";

    int groupAllElements = 9999;
    structure.GroupCreate(groupAllElements, eGroupId::Elements);
    structure.GroupAddElementsTotal(groupAllElements);
    structure.AddVisualizationComponent(groupAllElements, eVisualizeWhat::DISPLACEMENTS);
    structure.AddVisualizationComponent(groupAllElements, eVisualizeWhat::ENGINEERING_STRAIN);
    structure.AddVisualizationComponent(groupAllElements, eVisualizeWhat::ENGINEERING_STRESS);

    structure.GetLogger() << "*********************************** \n"
                          << "**      integration scheme       ** \n"
                          << "*********************************** \n\n";


    NuTo::NewmarkFeti<EigenSolver> myIntegrationScheme(&structure);

    boost::filesystem::path resultPath(boost::filesystem::path(getenv("HOME")).string() +
                                       std::string("/results/feti/") + std::to_string(structure.mRank));


    myIntegrationScheme.SetTimeStep(timeStep);
    myIntegrationScheme.SetMaxNumIterations(maxIterations);
    myIntegrationScheme.SetMinTimeStep(minTimeStep);
    myIntegrationScheme.SetMaxTimeStep(maxTimeStep);
    myIntegrationScheme.SetAutomaticTimeStepping(automaticTimeStepping);
    myIntegrationScheme.SetResultDirectory(resultPath.string(), true);
    myIntegrationScheme.SetPerformLineSearch(performLineSearch);
    myIntegrationScheme.SetToleranceResidual(eDof::DISPLACEMENTS, toleranceDisp);
    myIntegrationScheme.SetToleranceIterativeSolver(1.e-4);
    myIntegrationScheme.SetIterativeSolver(NuTo::NewmarkFeti<EigenSolver>::eIterativeSolver::ProjectedGmres);
    myIntegrationScheme.SetFetiPreconditioner(NuTo::NewmarkFeti<EigenSolver>::eFetiPreconditioner::Lumped);

    Eigen::Matrix2d dispRHS;
    dispRHS(0, 0) = 0;
    dispRHS(1, 0) = simulationTime;
    dispRHS(0, 1) = 0;
    dispRHS(1, 1) = loadFactor;

    myIntegrationScheme.SetTimeDependentLoadCase(loadId, dispRHS);

    structure.GetLogger() << "*********************************** \n"
                          << "**      solve                    ** \n"
                          << "*********************************** \n\n";


    myIntegrationScheme.Solve(simulationTime);
    structure.GetLogger() << "Total number of Dofs: \t" << structure.GetNumTotalDofs() << "\n\n";
}
