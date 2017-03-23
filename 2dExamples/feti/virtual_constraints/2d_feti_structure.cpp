
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
constexpr double youngsModulus = 1.0e4;
constexpr double poissonsRatio = 0.2;


// integration
constexpr bool performLineSearch = true;
constexpr bool automaticTimeStepping = true;
constexpr double timeStep = 1e-0;
constexpr double minTimeStep = 1e-1;
constexpr double maxTimeStep = 1e-0;
constexpr double toleranceDisp = 1e-6;
constexpr double simulationTime = 1.0;
constexpr double loadFactor = 10.0;
constexpr double maxIterations = 1;


const Eigen::Vector2d directionX = Eigen::Vector2d::UnitX();
const Eigen::Vector2d directionY = Eigen::Vector2d::UnitY();

void AssignSection(NuTo::StructureFeti& structure);
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

    const int interpolationTypeId = structure.InterpolationTypeCreate(eShapeType::QUAD2D);
    structure.InterpolationTypeAdd(interpolationTypeId, eDof::COORDINATES, eTypeOrder::EQUIDISTANT1);
    structure.InterpolationTypeAdd(interpolationTypeId, eDof::DISPLACEMENTS, eTypeOrder::EQUIDISTANT1);

    structure.ImportMeshJson(meshFile, interpolationTypeId);


    AssignMaterial(structure);
    AssignSection(structure);


    structure.GetLogger() << "**********************************************"
                          << "\n";
    structure.GetLogger() << "**  create node groups                      **"
                          << "\n";
    structure.GetLogger() << "**********************************************"
                          << "\n\n";

    Eigen::VectorXd nodeCoords(2);

    int groupNodesLeftBoundary = structure.GroupCreate(eGroupId::Nodes);
    structure.GroupAddNodeCoordinateRange(groupNodesLeftBoundary, 0, -1.e-6, +1.e-6);

    int loadNodeGroup = structure.GroupCreate(eGroupId::Nodes);
    //    nodeCoords[0] = 60;
    //    nodeCoords[1] = 0;
    //    structure.GroupAddNodeRadiusRange(loadNodeGroup, nodeCoords, 0, 1.e-6);
    structure.GroupAddNodeCoordinateRange(loadNodeGroup, 0, 60 - 1.e-6, 60 + 1.e-6);

    structure.GetLogger() << "**********************************************"
                          << "\n";
    structure.GetLogger() << "**  virtual constraints                     **"
                          << "\n";
    structure.GetLogger() << "**********************************************"
                          << "\n\n";

    std::vector<int> nodeIdsBoundaries = structure.GroupGetMemberIds(groupNodesLeftBoundary);
    std::vector<int> nodeIdsLoads = structure.GroupGetMemberIds(loadNodeGroup);


    structure.ApplyVirtualConstraints(nodeIdsBoundaries, nodeIdsLoads);

    structure.GetLogger() << "**********************************************"
                          << "\n";
    structure.GetLogger() << "**  real constraints                        **"
                          << "\n";
    structure.GetLogger() << "**********************************************"
                          << "\n\n";

    structure.NodeBuildGlobalDofs(__PRETTY_FUNCTION__);
    structure.ApplyConstraintsTotalFeti(groupNodesLeftBoundary);

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
        dofIdAndPrescribedDisplacementMap.emplace(dofIds[0], 1.);
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
    //    structure.AddVisualizationComponent(groupAllElements, eVisualizeWhat::DAMAGE);


    structure.GetLogger() << "**********************************************"
                          << "\n";
    structure.GetLogger() << "**  integration sheme                       **"
                          << "\n";
    structure.GetLogger() << "**********************************************"
                          << "\n\n";


    NuTo::NewmarkFeti<EigenSolver> myIntegrationScheme(&structure);
    boost::filesystem::path resultPath(std::string("/home/phuschke/results/feti/" + std::to_string(structure.mRank)));

    myIntegrationScheme.SetTimeStep(timeStep);
    myIntegrationScheme.SetMaxNumIterations(maxIterations);
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

    //
    //    structure.GetLogger().CloseFile();
    //    structure.GetLogger().OpenFile("PreconditionerTest" + std::to_string(rank));
    //
    //    structure.NodeInfo(10);
    //    structure.GetLogger() << "B=\n" << structure.GetConnectivityMatrix() << "\n\n";
    //    structure.GetLogger() << "RBM=\n" << structure.GetRigidBodyModes() << "\n\n";
    //    structure.GetLogger() << std::setprecision(5) << "Projection=\n" << structure.GetProjectionMatrix() << "\n\n";
    //
    //    std::map<NuTo::eStructureOutput, NuTo::StructureOutputBase*> evalHessian0;
    //
    //    NuTo::StructureOutputBlockMatrix hessian0(structure.GetDofStatus(), true);
    //    evalHessian0[NuTo::eStructureOutput::HESSIAN0] = &hessian0;
    //
    //    NuTo::ConstitutiveInputMap inputMap;
    //    inputMap[NuTo::Constitutive::eInput::CALCULATE_STATIC_DATA] =
    //            std::make_unique<NuTo::ConstitutiveCalculateStaticData>(NuTo::eCalculateStaticData::EULER_BACKWARD);
    //
    //    structure.Evaluate(inputMap, evalHessian0);
    //
    //    structure.GetLogger() << "Kjj=\n" << hessian0.JJ.ExportToEigenSparseMatrix() << "\n\n";
    //
    //    Eigen::Matrix2d bla = 3*Eigen::Matrix2d::Ones();
    //    Eigen::Matrix2d bla2;
    //    structure.GetLogger() << "F_I=\n" << bla << "\n\n";
    //
    //    MPI_Allreduce(bla.data(), bla2.data(), 4, MPI_DOUBLE, MPI_PROD, MPI_COMM_WORLD);
    //    structure.GetLogger() << "F_I=\n" << bla2 << "\n\n";
}


void AssignSection(NuTo::StructureFeti& structure)
{

    auto section = NuTo::SectionPlane::Create(thickness, true);
    structure.ElementTotalSetSection(section);
}

void AssignMaterial(NuTo::StructureFeti& structure)
{
    structure.GetLogger() << "***********************************"
                          << "\n";
    structure.GetLogger() << "**      Material                 **"
                          << "\n";
    structure.GetLogger() << "***********************************"
                          << "\n\n";

    int material00 = structure.ConstitutiveLawCreate(eConstitutiveType::LINEAR_ELASTIC_ENGINEERING_STRESS);

    structure.ConstitutiveLawSetParameterDouble(material00, eConstitutiveParameter::YOUNGS_MODULUS, youngsModulus);
    structure.ConstitutiveLawSetParameterDouble(material00, eConstitutiveParameter::POISSONS_RATIO, poissonsRatio);

    structure.ElementTotalSetConstitutiveLaw(material00);
}
