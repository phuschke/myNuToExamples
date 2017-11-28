
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
#include <mechanics/feti/FetiDirichletPreconditioner.h>
#include "mechanics/mesh/MeshGenerator.h"
#include "mechanics/groups/Group.h"
using namespace NuTo;
constexpr int dim = 2;
using Eigen::VectorXd;
using Eigen::MatrixXd;
using EigenSolver = Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>>;
using FetiScaling = NewmarkFeti<EigenSolver>::eFetiScaling;
constexpr double thickness = 1.0;

// material
constexpr double nonlocalRadius = 0.2; // mm
constexpr double youngsModulus = 4.0e4;
constexpr double poissonsRatio = 0.2;
constexpr double tensileStrength = 3;
constexpr double compressiveStrength = 30;
constexpr double fractureEnergy = 0.01;
constexpr double alpha = 0.99;

// integration
constexpr bool performLineSearch = false;
constexpr bool automaticTimeStepping = false;
constexpr double timeStep = 1e-1;
constexpr double minTimeStep = 1e-5;
constexpr double maxTimeStep = 1e-1;

constexpr double toleranceDisp = 1e-8;
constexpr double toleranceNlEqStrain = 1e-8;
constexpr double tolerance = 1e-5;

constexpr double simulationTime = 1.0;
constexpr double loadFactor = -0.03;
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
    meshDimensions.push_back(10.);
    meshDimensions.push_back(60.);

    std::vector<int> numElements;
    numElements.push_back(5);
    numElements.push_back(120);

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

//    const int materialId = structure.ConstitutiveLawCreate(eConstitutiveType::LINEAR_ELASTIC_ENGINEERING_STRESS);
//    structure.ConstitutiveLawSetParameterDouble(materialId, eConstitutiveParameter::YOUNGS_MODULUS, youngsModulus);
//    structure.ConstitutiveLawSetParameterDouble(materialId, eConstitutiveParameter::POISSONS_RATIO, poissonsRatio);


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

    coordinates << 0., 60.;
    const int groupNodesBottomRight = structure.GroupCreate(eGroupId::Nodes);
    structure.GroupAddNodeRadiusRange(groupNodesBottomRight, coordinates, 0, tolerance);

    const int groupNodesBoundary = structure.GroupUnion(groupNodesBottomLeft, groupNodesBottomRight);

    coordinates[0] = 10.;
    coordinates[1] = 30;
    const auto& groupNodeLoad = structure.GroupGetNodeRadiusRange(coordinates);


    structure.GetLogger() << "*********************************** \n"
                          << "**      virtual constraints      ** \n"
                          << "*********************************** \n\n";

    std::vector<int> nodeIdsBoundaries = structure.GroupGetMemberIds(groupNodesBoundary);
    std::vector<int> nodeIdsLoads = groupNodeLoad.GetMemberIds();

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
        boundaryDofIds.push_back(dofIds[0]);
    }

    structure.ApplyConstraintsTotalFeti(boundaryDofIds);


    structure.GetLogger() << "*********************************** \n"
                          << "**      load                     ** \n"
                          << "*********************************** \n\n";


    // prescribe displacement of loadNodeGroup in Y direction
    std::map<int, double> dofIdAndPrescribedDisplacementMap;
    std::vector<int> nodeIds = groupNodeLoad.GetMemberIds();
    for (auto const& nodeId : nodeIds)
    {
        std::vector<int> dofIds = structure.NodeGetDofIds(nodeId, eDof::DISPLACEMENTS);
        dofIdAndPrescribedDisplacementMap.emplace(dofIds[0], 1.);
    }

    structure.ApplyPrescribedDisplacements(dofIdAndPrescribedDisplacementMap);

    int loadId = structure.LoadCreateNodeGroupForce(&groupNodeLoad, Vector2d::UnitX(), 0.);

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

    NuTo::NewmarkFeti<EigenSolver> newmarkFeti(&structure);

    boost::filesystem::path resultPath("results_three_point_bending_" + std::to_string(structure.mRank));

    newmarkFeti.SetTimeStep(timeStep);
    newmarkFeti.SetMaxNumIterations(maxIterations);
    newmarkFeti.SetMinTimeStep(minTimeStep);
    newmarkFeti.SetMaxTimeStep(maxTimeStep);
    newmarkFeti.SetAutomaticTimeStepping(automaticTimeStepping);
    newmarkFeti.PostProcessing().SetResultDirectory(resultPath.string(), true);
    newmarkFeti.SetPerformLineSearch(performLineSearch);
    newmarkFeti.SetToleranceResidual(eDof::DISPLACEMENTS, toleranceDisp);
    newmarkFeti.SetToleranceResidual(eDof::NONLOCALEQSTRAIN, toleranceNlEqStrain);
    newmarkFeti.SetToleranceIterativeSolver(1.e-6);
    newmarkFeti.SetIterativeSolver(NuTo::NewmarkFeti<EigenSolver>::eIterativeSolver::ProjectedGmres);
    newmarkFeti.SetFetiPreconditioner(std::make_unique<NuTo::FetiDirichletPreconditioner>());
    newmarkFeti.SetFetiScaling(FetiScaling::Multiplicity);

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
