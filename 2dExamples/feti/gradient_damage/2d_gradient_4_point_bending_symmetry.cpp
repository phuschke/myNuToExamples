
#include <mpi.h>


#include "mechanics/feti/StructureFeti.h"

#include <ctime>
#include <chrono>
#include <mechanics/constitutive/damageLaws/DamageLawExponential.h>
#include "mechanics/feti/NewmarkFeti.h"
#include "../../../EnumsAndTypedefs.h"

#include "mechanics/nodes/NodeBase.h"

#include "boost/filesystem.hpp"

#include "mechanics/sections/SectionPlane.h"

constexpr int dim = 2;
using Eigen::VectorXd;
using Eigen::MatrixXd;

// using EigenSolver = Eigen::SparseQR<Eigen::SparseMatrix<double>,Eigen::COLAMDOrdering<int>>;
// using EigenSolver = Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>>;
// using EigenSolver = Eigen::PardisoLU<Eigen::SparseMatrix<double>>;
using EigenSolver = Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>>;


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
constexpr double loadFactor = -0.5;
constexpr double maxIterations = 10;


const Eigen::Vector2d directionX = Eigen::Vector2d::UnitX();
const Eigen::Vector2d directionY = Eigen::Vector2d::UnitY();


int main(int argc, char* argv[])
{
    MPI_Init(&argc, &argv);

    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    NuTo::StructureFeti structure(dim);
    structure.SetNumTimeDerivatives(0);
    structure.SetVerboseLevel(10);
    structure.SetShowTime(false);
    structure.GetLogger().OpenFile("output" + std::to_string(rank));
    structure.GetLogger().SetQuiet(true);


    std::string meshFile = argv[1] + std::to_string(rank);
    structure.GetLogger() << meshFile << "\n";

    const int interpolationTypeId = structure.InterpolationTypeCreate(eShapeType::QUAD2D);
    structure.InterpolationTypeAdd(interpolationTypeId, eDof::COORDINATES, eTypeOrder::EQUIDISTANT1);
    structure.InterpolationTypeAdd(interpolationTypeId, eDof::DISPLACEMENTS, eTypeOrder::EQUIDISTANT1);
    structure.InterpolationTypeAdd(interpolationTypeId, eDof::NONLOCALEQSTRAIN, eTypeOrder::EQUIDISTANT1);

    structure.ImportMeshJson(meshFile, interpolationTypeId);

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


    structure.GetLogger() << "*********************************** \n"
                          << "**      section                  ** \n"
                          << "*********************************** \n\n";

    auto section = NuTo::SectionPlane::Create(thickness, true);
    structure.ElementTotalSetSection(section);

    structure.GetLogger() << "*********************************** \n"
                          << "**      node groups              ** \n"
                          << "*********************************** \n\n";

    Eigen::VectorXd coordinates(dim);
    coordinates << 25., 0.;

    const int groupNodesBottomLeft = structure.GroupCreate(eGroupId::Nodes);
    structure.GroupAddNodeRadiusRange(groupNodesBottomLeft, coordinates, 0, tolerance);

    const int groupNodesTmp01 = structure.GroupCreate(NuTo::eGroupId::Nodes);
    structure.GroupAddNodeCoordinateRange(groupNodesTmp01, 0, 250 - tolerance, 250 + tolerance);

    const int groupNodesTmp02 = structure.GroupCreate(NuTo::eGroupId::Nodes);
    structure.GroupAddNodeCoordinateRange(groupNodesTmp02, 1, 10 - tolerance, 100 + tolerance);

    const int groupNodesSymmetry = structure.GroupIntersection(groupNodesTmp01, groupNodesTmp02);

    const int groupNodesBoundary = structure.GroupUnion(groupNodesBottomLeft, groupNodesSymmetry);

    coordinates << 175., 100.;
    const int groupNodesLoad = structure.GroupCreate(NuTo::eGroupId::Nodes);
    structure.GroupAddNodeRadiusRange(groupNodesLoad, coordinates, 0, tolerance);

    structure.GetLogger() << "*********************************** \n"
                          << "**      virtual constraints      ** \n"
                          << "*********************************** \n\n";

    structure.Info();
    std::vector<int> nodeIdsBoundaries = structure.GroupGetMemberIds(groupNodesBoundary);
    std::vector<int> nodeIdsLoads = structure.GroupGetMemberIds(groupNodesLoad);


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
        boundaryDofIds.push_back(dofIds[1]);
    }

    std::vector<int> nodeIdsBoundaryRight = structure.GroupGetMemberIds(groupNodesSymmetry);
    for (const int nodeId : nodeIdsBoundaryRight)
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
    structure.AddVisualizationComponent(groupAllElements, eVisualizeWhat::NONLOCAL_EQ_STRAIN);
    structure.AddVisualizationComponent(groupAllElements, eVisualizeWhat::ENGINEERING_STRAIN);
    structure.AddVisualizationComponent(groupAllElements, eVisualizeWhat::ENGINEERING_STRESS);
    structure.AddVisualizationComponent(groupAllElements, eVisualizeWhat::DAMAGE);


    structure.GetLogger() << "*********************************** \n"
                          << "**      integration scheme       ** \n"
                          << "*********************************** \n\n";

    NuTo::NewmarkFeti<EigenSolver> myIntegrationScheme(&structure);
    boost::filesystem::path resultPath(boost::filesystem::path(getenv("HOME")).string() +
                                       std::string("/results/feti/gradient_damage/") + std::to_string(structure.mRank));

    myIntegrationScheme.SetTimeStep(timeStep);
    myIntegrationScheme.SetMaxNumIterations(maxIterations);
    myIntegrationScheme.SetMinTimeStep(minTimeStep);
    myIntegrationScheme.SetMaxTimeStep(maxTimeStep);
    myIntegrationScheme.SetAutomaticTimeStepping(automaticTimeStepping);
    myIntegrationScheme.SetResultDirectory(resultPath.string(), true);
    myIntegrationScheme.SetPerformLineSearch(performLineSearch);
    myIntegrationScheme.SetToleranceResidual(eDof::DISPLACEMENTS, toleranceDisp);
    myIntegrationScheme.SetToleranceResidual(eDof::NONLOCALEQSTRAIN, toleranceNlEqStrain);
    myIntegrationScheme.SetIterativeSolver(NuTo::NewmarkFeti<EigenSolver>::eIterativeSolver::ProjectedGmres);
    myIntegrationScheme.SetToleranceIterativeSolver(tolerance);

    Eigen::Matrix2d dispRHS;
    dispRHS(0, 0) = 0;
    dispRHS(1, 0) = simulationTime;
    dispRHS(0, 1) = 0;
    dispRHS(1, 1) = loadFactor;

    myIntegrationScheme.SetTimeDependentLoadCase(loadId, dispRHS);

    structure.GetLogger() << "*********************************** \n"
                          << "**      solve                    ** \n"
                          << "*********************************** \n\n";

    structure.GetLogger() << "Total number of Dofs: \t" << structure.GetNumTotalDofs() << "\n\n";

    myIntegrationScheme.Solve(simulationTime);

    MPI_Finalize();
}
