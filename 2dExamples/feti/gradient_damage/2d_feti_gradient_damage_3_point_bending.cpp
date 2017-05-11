
#include <mpi.h>
#include <boost/mpi.hpp>

#include "mechanics/feti/StructureFeti.h"
#include <ctime>
#include <chrono>
#include "mechanics/feti/NewmarkFeti.h"
#include "../../../EnumsAndTypedefs.h"

#include "mechanics/nodes/NodeBase.h"
#include "mechanics/constitutive/damageLaws/DamageLawExponential.h"
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
constexpr double nonlocalRadius = 0.5; // mm
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
constexpr double loadFactor = -0.1;
constexpr double maxIterations = 10;


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
    numElements.push_back(20);
    numElements.push_back(10);

    const auto importContainer = structure.CreateRectangularMesh2D(meshDimensions, numElements);

    const int interpolationTypeIdDamage = importContainer.second;
    structure.InterpolationTypeAdd(interpolationTypeIdDamage, eDof::DISPLACEMENTS, eTypeOrder::EQUIDISTANT1);
    structure.InterpolationTypeAdd(interpolationTypeIdDamage, eDof::NONLOCALEQSTRAIN, eTypeOrder::EQUIDISTANT1);
    structure.ElementTotalConvertToInterpolationType();

    structure.GetLogger() << "*********************************** \n"
                          << "**      material                 ** \n"
                          << "*********************************** \n\n";

    const int materialId = structure.ConstitutiveLawCreate(eConstitutiveType::GRADIENT_DAMAGE_ENGINEERING_STRESS);
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

    auto section = NuTo::SectionPlane::Create(thickness, true);
    structure.ElementTotalSetSection(section);


    structure.GetLogger() << "*********************************** \n"
                          << "**      node groups              ** \n"
                          << "*********************************** \n\n";

    Eigen::VectorXd coordinates(dim);

    coordinates << 0., 0.;
    const int groupNodesBottomLeft = structure.GroupCreate(eGroupId::Nodes);
    structure.GroupAddNodeRadiusRange(groupNodesBottomLeft, coordinates, 0, tolerance);

    coordinates << 60., 0.;
    const int groupNodesBottomRight = structure.GroupCreate(eGroupId::Nodes);
    structure.GroupAddNodeRadiusRange(groupNodesBottomRight, coordinates, 0, tolerance);

    const int groupNodesBoundary = structure.GroupUnion(groupNodesBottomLeft, groupNodesBottomRight);

    int loadNodeGroup = structure.GroupCreate(eGroupId::Nodes);
    coordinates[0] = 30;
    coordinates[1] = 10;
    structure.GroupAddNodeRadiusRange(loadNodeGroup, coordinates, 0, 1.e-6);


    structure.GetLogger() << "*********************************** \n"
                          << "**      virtual constraints      ** \n"
                          << "*********************************** \n\n";

    std::vector<int> nodeIdsBoundaries = structure.GroupGetMemberIds(groupNodesBoundary);
    std::vector<int> nodeIdsLoads = structure.GroupGetMemberIds(loadNodeGroup);

    //    std::cout << nodeIdsBoundaries << "\n";

    structure.ApplyVirtualConstraints(nodeIdsBoundaries, nodeIdsLoads);

    structure.GetLogger() << "*********************************** \n"
                          << "**      real constraints         ** \n"
                          << "*********************************** \n\n";

    structure.NodeBuildGlobalDofs(__PRETTY_FUNCTION__);

    std::vector<int> boundaryDofIds;
    std::vector<int> nodeIdsBoundaryLeft = structure.GroupGetMemberIds(groupNodesBottomLeft);
    for (const int nodeId : nodeIdsBoundaryLeft)
    {
        std::vector<int> dofIds = structure.NodeGetDofIds(nodeId, eDof::DISPLACEMENTS);
        boundaryDofIds.push_back(dofIds[0]);
        boundaryDofIds.push_back(dofIds[1]);
    }

    std::vector<int> nodeIdsBoundaryRight = structure.GroupGetMemberIds(groupNodesBottomRight);
    for (const int nodeId : nodeIdsBoundaryRight)
    {
        std::vector<int> dofIds = structure.NodeGetDofIds(nodeId, eDof::DISPLACEMENTS);
        boundaryDofIds.push_back(dofIds[1]);
    }

    structure.ApplyConstraintsTotalFeti(boundaryDofIds);


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
    structure.AddVisualizationComponent(groupAllElements, eVisualizeWhat::DAMAGE);

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
    myIntegrationScheme.SetToleranceResidual(eDof::NONLOCALEQSTRAIN, toleranceNlEqStrain);
    myIntegrationScheme.SetToleranceIterativeSolver(1.e-6);
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
