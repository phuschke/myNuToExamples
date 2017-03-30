
#include <mpi.h>
#include <boost/mpi.hpp>

#include "mechanics/feti/StructureFeti.h"
#include <ctime>
#include <chrono>
#include <mechanics/constitutive/damageLaws/DamageLawExponential.h>
#include "mechanics/feti/NewmarkFeti.h"
#include "../../../EnumsAndTypedefs.h"

#include "mechanics/nodes/NodeBase.h"

#include "boost/filesystem.hpp"


using Eigen::VectorXd;
using Eigen::MatrixXd;
// using EigenSolver = Eigen::SparseQR<Eigen::SparseMatrix<double>,Eigen::COLAMDOrdering<int>>;
// using EigenSolver = Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>>;
// using EigenSolver = Eigen::PardisoLU<Eigen::SparseMatrix<double>>;
using EigenSolver = Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>>;

constexpr int dim = 3;

// material
// constexpr   double      youngsModulus               = 4.0e4 * 1.e-4;
// constexpr   double      poissonsRatio               = 0.2;
// constexpr   double      tensileStrength             = 3 * 1.e-4;
// constexpr   double      compressiveStrength         = 30 * 1.e-4;
// constexpr   double      fractureEnergy              = 1.e-4;
// constexpr   double      nonlocalRadius              = 0.05;                   // mm

constexpr double youngsModulus = 4.0e4;
constexpr double poissonsRatio = 0.2;
constexpr double tensileStrength = 3;
constexpr double compressiveStrength = 30;
constexpr double fractureEnergy = 0.01;
constexpr double nonlocalRadius = 2; // mm
constexpr double alpha = 0.95; // mm

// integration
constexpr double timeStep = 1.e-0;
constexpr double minTimeStep = 1.e-2 * timeStep;
constexpr double toleranceDisp = 1.e-6;
constexpr double toleranceNlEqStrain = 1.e-6;
constexpr double simulationTime = 30.0;
constexpr double loadFactor = -0.2;
constexpr bool performLineSearch = true;
constexpr bool automaticTimeStepping = true;

const Eigen::Vector3d directionX = Eigen::Vector3d::UnitX();
const Eigen::Vector3d directionY = Eigen::Vector3d::UnitY();
const Eigen::Vector3d directionZ = Eigen::Vector3d::UnitZ();


void AssignMaterial(NuTo::StructureFeti& structure);


int main(int argc, char* argv[])
{
    boost::mpi::environment env(argc, argv);
    boost::mpi::communicator world;

    const int rank = world.rank();

    NuTo::StructureFeti structure(dim);
    structure.SetNumTimeDerivatives(0);
    structure.SetVerboseLevel(10);
    structure.SetShowTime(true);
    structure.GetLogger().OpenFile("output" + std::to_string(rank));
    structure.GetLogger().SetQuiet(true);


    std::string meshFile = argv[1] + std::to_string(rank);
    structure.GetLogger() << meshFile << "\n";

    const int interpolationTypeId = structure.InterpolationTypeCreate(eShapeType::BRICK3D);
    structure.InterpolationTypeAdd(interpolationTypeId, eDof::COORDINATES, eTypeOrder::EQUIDISTANT1);
    structure.InterpolationTypeAdd(interpolationTypeId, eDof::DISPLACEMENTS, eTypeOrder::EQUIDISTANT1);
    structure.InterpolationTypeAdd(interpolationTypeId, eDof::NONLOCALEQSTRAIN, eTypeOrder::EQUIDISTANT1);

    structure.ImportMeshJson(meshFile, interpolationTypeId);


    AssignMaterial(structure);

    structure.GetLogger() << "**********************************************"
                          << "\n";
    structure.GetLogger() << "**  create node groups                      **"
                          << "\n";
    structure.GetLogger() << "**********************************************"
                          << "\n\n";

    Eigen::VectorXd nodeCoords(dim);

    int groupNodesBoundary = structure.GroupCreate(eGroupId::Nodes);
    nodeCoords << 0., 0., 0.;
    structure.GroupAddNodeCylinderRadiusRange(groupNodesBoundary, nodeCoords, directionZ, 0, 1.e-3);

    nodeCoords << 100., 0., 0.;
    structure.GroupAddNodeCylinderRadiusRange(groupNodesBoundary, nodeCoords, directionZ, 0, 1.e-3);

    int loadNodeGroup = structure.GroupCreate(eGroupId::Nodes);

    nodeCoords << 50., 10., 0.;
    structure.GroupAddNodeCylinderRadiusRange(loadNodeGroup, nodeCoords, directionZ, 0, 1.e-3);


    structure.GetLogger() << "**********************************************"
                          << "\n";
    structure.GetLogger() << "**  virtual constraints                     **"
                          << "\n";
    structure.GetLogger() << "**********************************************"
                          << "\n\n";

    std::vector<int> nodeIdsBoundaries = structure.GroupGetMemberIds(groupNodesBoundary);

    std::vector<int> nodeIdsLoads = structure.GroupGetMemberIds(loadNodeGroup);


    structure.ApplyVirtualConstraints(nodeIdsBoundaries, nodeIdsLoads);

    structure.GetLogger() << "**********************************************"
                          << "\n";
    structure.GetLogger() << "**  real constraints                        **"
                          << "\n";
    structure.GetLogger() << "**********************************************"
                          << "\n\n";

    structure.NodeBuildGlobalDofs(__PRETTY_FUNCTION__);
    //        structure.ApplyConstraintsTotalFeti(groupNodesBoundary);

    //////////////////////////////////////////
    //// Testing a new way to apply constraints
    int groupNodesBoundaryLeft = structure.GroupCreate(eGroupId::Nodes);
    nodeCoords << 0., 0., 0.;
    structure.GroupAddNodeCylinderRadiusRange(groupNodesBoundaryLeft, nodeCoords, directionZ, 0, 1.e-3);

    std::vector<int> boundaryDofIds;
    std::vector<int> nodeIdsBoundaryLeft = structure.GroupGetMemberIds(groupNodesBoundaryLeft);
    for (const int nodeId : nodeIdsBoundaryLeft)
    {
        std::vector<int> dofIds = structure.NodeGetDofIds(nodeId, eDof::DISPLACEMENTS);

        for (const auto& id : dofIds)
            boundaryDofIds.push_back(id);
    }

    int groupNodesBoundaryRight = structure.GroupCreate(eGroupId::Nodes);
    nodeCoords << 100., 0., 0.;
    structure.GroupAddNodeCylinderRadiusRange(groupNodesBoundaryRight, nodeCoords, directionZ, 0, 1.e-3);
    std::vector<int> nodeIdsBoundaryRight = structure.GroupGetMemberIds(groupNodesBoundaryRight);
    for (const int nodeId : nodeIdsBoundaryRight)
    {
        std::vector<int> dofIds = structure.NodeGetDofIds(nodeId, eDof::DISPLACEMENTS);


        boundaryDofIds.push_back(dofIds[1]);
        boundaryDofIds.push_back(dofIds[2]);
    }

    //////////////////////////////////////////

    structure.ApplyConstraintsTotalFeti(boundaryDofIds);

    structure.GetLogger() << "**********************************************"
                          << "\n";
    structure.GetLogger() << "**  load                                    **"
                          << "\n";
    structure.GetLogger() << "**********************************************"
                          << "\n\n";


    //    structure.NodeInfo(10);
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
    structure.AddVisualizationComponent(groupAllElements, eVisualizeWhat::ENGINEERING_STRAIN);
    structure.AddVisualizationComponent(groupAllElements, eVisualizeWhat::ENGINEERING_STRESS);
    structure.AddVisualizationComponent(groupAllElements, eVisualizeWhat::DAMAGE);


    structure.GetLogger() << "**********************************************"
                          << "\n";
    structure.GetLogger() << "**  integration sheme                       **"
                          << "\n";
    structure.GetLogger() << "**********************************************"
                          << "\n\n";


    NuTo::NewmarkFeti<EigenSolver> myIntegrationScheme(&structure);
    boost::filesystem::path resultPath(std::string("/home/phuschke/results/feti/3d/" + std::string(argv[2]) + "/" +
                                                   std::to_string(structure.mRank)));

    myIntegrationScheme.SetTimeStep(timeStep);
    myIntegrationScheme.SetMinTimeStep(minTimeStep);
    myIntegrationScheme.SetAutomaticTimeStepping(automaticTimeStepping);
    myIntegrationScheme.SetResultDirectory(resultPath.string(), true);
    myIntegrationScheme.SetPerformLineSearch(performLineSearch);

    myIntegrationScheme.SetToleranceResidual(eDof::DISPLACEMENTS, toleranceDisp);
    myIntegrationScheme.SetToleranceResidual(eDof::NONLOCALEQSTRAIN, toleranceNlEqStrain);
    myIntegrationScheme.SetIterativeSolver(
            NuTo::NewmarkFeti<EigenSolver>::eIterativeSolver::BiconjugateGradientStabilized);

    Eigen::Matrix2d dispRHS;
    dispRHS(0, 0) = 0;
    dispRHS(1, 0) = simulationTime;
    dispRHS(0, 1) = 0;
    dispRHS(1, 1) = loadFactor;

    //    myIntegrationScheme.AddTimeDependentConstraint(loadId, dispRHS);
    myIntegrationScheme.SetTimeDependentLoadCase(loadId, dispRHS);

    structure.GetLogger() << "***********************************"
                          << "\n";
    structure.GetLogger() << "**      Solve                    **"
                          << "\n";
    structure.GetLogger() << "***********************************"
                          << "\n\n";


    myIntegrationScheme.Solve(simulationTime);
    structure.GetLogger() << "Total number of Dofs: \t" << structure.GetNumTotalDofs() << "\n\n";
}


void AssignMaterial(NuTo::StructureFeti& structure)
{
    structure.GetLogger() << "***********************************"
                          << "\n";
    structure.GetLogger() << "**      Material                 **"
                          << "\n";
    structure.GetLogger() << "***********************************"
                          << "\n\n";

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
}
