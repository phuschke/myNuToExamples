
#include <mpi.h>
#include <boost/mpi.hpp>
#include "mechanics/constitutive/laws/PhaseField.h"
#include "mechanics/feti/StructureFeti.h"
#include <ctime>
#include <chrono>
#include "mechanics/feti/NewmarkFeti.h"
#include "../../../EnumsAndTypedefs.h"
#include "mechanics/sections/SectionPlane.h"
#include "mechanics/nodes/NodeBase.h"

#include "boost/filesystem.hpp"

constexpr int dim = 2;
using Eigen::VectorXd;
using Eigen::MatrixXd;
using EigenSolver = Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>>;
// using EigenSolver = Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>>;
// using EigenSolver = Eigen::PardisoLU<Eigen::SparseMatrix<double>>;
// using EigenSolver = Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>>;


constexpr int dimension = 2;
constexpr double thickness = 1.0;


constexpr double toleranceDisp = 1e-4;
constexpr double toleranceCrack = 1e-4;
constexpr double toleranceNlEqStrain = 1e-6;

const Eigen::Vector2d directionX = Eigen::Vector2d::UnitX();
const Eigen::Vector2d directionY = Eigen::Vector2d::UnitY();

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//          gradient enhanced damage model parameters
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//// material
// constexpr double nonlocalRadius = 5; // mm
//
// constexpr double youngsModulus       = 4.0; // N/mm^2
// constexpr double poissonsRatio       = 0.2;
// constexpr double fractureEnergy      = 1e-6; // N/mm
// constexpr double compressiveStrength = 30.e-4; // N/mm
// constexpr double tensileStrength     = 3.e-4; // N/mm
//
//
//// integration
// constexpr bool performLineSearch     = true;
// constexpr bool automaticTimeStepping = true;
// constexpr double timeStep            = 1e-1;
// constexpr double minTimeStep         = 1e-3;
// constexpr double maxTimeStep         = 1e-1;
//
//
// constexpr double simulationTime = 1.0;
// constexpr double loadFactor     = -2e-2;
// constexpr double maxIterations  = 10;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//          phase-field model parameters
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// material
constexpr double youngsModulus = 2.0; // N/mm^2
constexpr double poissonsRatio = 0.3;
constexpr double lengthScaleParameter = 1.0e-1; // mm
constexpr double fractureEnergy = 2.0e-4; // N/mm
constexpr double artificialViscosity = 0; // Ns/mm^2
constexpr ePhaseFieldEnergyDecomposition energyDecomposition = ePhaseFieldEnergyDecomposition::ISOTROPIC;


// integration
constexpr double timeStep = 1e-6;
constexpr double minTimeStep = 1e-6;
constexpr double maxTimeStep = 1e-1;
constexpr double simulationTime = 1.0;
constexpr double loadFactor = -0.5;
constexpr bool performLineSearch = true;
constexpr bool automaticTimeStepping = true;


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
    //    structure.InterpolationTypeAdd(interpolationTypeId, eDof::NONLOCALEQSTRAIN, eTypeOrder::EQUIDISTANT1);
    structure.InterpolationTypeAdd(interpolationTypeId, eDof::CRACKPHASEFIELD, eTypeOrder::EQUIDISTANT1);

    structure.ImportMeshJson(meshFile, interpolationTypeId);


    AssignMaterial(structure);
    AssignSection(structure);


    structure.GetLogger() << "*********************************** \n"
                          << "**      create node groups       ** \n"
                          << "*********************************** \n\n";

    Eigen::VectorXd nodeCoords(dim);

    int groupNodesBoundary = structure.GroupCreate(eGroupId::Nodes);
    nodeCoords << 0., 0.;
    structure.GroupAddNodeRadiusRange(groupNodesBoundary, nodeCoords, 0, 1.e-6);

    nodeCoords << 60., 0.;
    structure.GroupAddNodeRadiusRange(groupNodesBoundary, nodeCoords, 0, 1.e-6);

    int groupNodesLoad = structure.GroupCreate(eGroupId::Nodes);

    nodeCoords << 30., 10.;
    structure.GroupAddNodeRadiusRange(groupNodesLoad, nodeCoords, 0, 1.e-6);

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

    int groupNodesBoundaryLeft = structure.GroupCreate(eGroupId::Nodes);
    nodeCoords << 0., 0.;
    structure.GroupAddNodeRadiusRange(groupNodesBoundaryLeft, nodeCoords, 0, 1.e-6);

    std::vector<int> boundaryDofIds;
    std::vector<int> nodeIdsBoundaryLeft = structure.GroupGetMemberIds(groupNodesBoundaryLeft);
    for (const int nodeId : nodeIdsBoundaryLeft)
    {
        std::vector<int> dofIds = structure.NodeGetDofIds(nodeId, eDof::DISPLACEMENTS);

        for (const auto& id : dofIds)
            boundaryDofIds.push_back(id);
    }

    int groupNodesBoundaryRight = structure.GroupCreate(eGroupId::Nodes);
    nodeCoords << 60., 0.;
    structure.GroupAddNodeRadiusRange(groupNodesBoundaryRight, nodeCoords, 0, 1.e-6);
    std::vector<int> nodeIdsBoundaryRight = structure.GroupGetMemberIds(groupNodesBoundaryRight);
    for (const int nodeId : nodeIdsBoundaryRight)
    {
        std::vector<int> dofIds = structure.NodeGetDofIds(nodeId, eDof::DISPLACEMENTS);
        boundaryDofIds.push_back(dofIds[1]);
    }

    //////////////////////////////////////////

    structure.ApplyConstraintsTotalFeti(boundaryDofIds);


    structure.GetLogger() << "*********************************** \n"
                          << "**      load                     ** \n"
                          << "*********************************** \n\n";


    //    structure.NodeInfo(10);
    // prescribe displacement of loadNodeGroup in Y direction
    std::map<int, double> dofIdAndPrescribedDisplacementMap;
    std::vector<int> nodeIds = structure.GroupGetMemberIds(groupNodesLoad);
    for (auto const& nodeId : nodeIds)
    {
        std::vector<int> dofIds = structure.NodeGetDofIds(nodeId, eDof::DISPLACEMENTS);
        dofIdAndPrescribedDisplacementMap.emplace(dofIds[1], 1.);
    }

    structure.ApplyPrescribedDisplacements(dofIdAndPrescribedDisplacementMap);

    int loadId = structure.LoadCreateNodeGroupForce(0, groupNodesLoad, directionY, 0);


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
    //        structure.AddVisualizationComponent(groupAllElements, eVisualizeWhat::DAMAGE);
    structure.AddVisualizationComponent(groupAllElements, eVisualizeWhat::CRACK_PHASE_FIELD);

    structure.GetLogger() << "*********************************** \n"
                          << "**      integration scheme       ** \n"
                          << "*********************************** \n\n";

    NuTo::NewmarkFeti<EigenSolver> myIntegrationScheme(&structure);
    boost::filesystem::path resultPath(std::string("/home/phuschke/results/feti/" + std::to_string(structure.mRank)));

    myIntegrationScheme.SetTimeStep(timeStep);

    myIntegrationScheme.SetMinTimeStep(minTimeStep);
    myIntegrationScheme.SetMaxTimeStep(maxTimeStep);
    myIntegrationScheme.SetAutomaticTimeStepping(automaticTimeStepping);
    myIntegrationScheme.SetResultDirectory(resultPath.string(), true);
    myIntegrationScheme.SetPerformLineSearch(performLineSearch);
    myIntegrationScheme.SetToleranceResidual(eDof::DISPLACEMENTS, toleranceDisp);
    myIntegrationScheme.SetToleranceResidual(eDof::CRACKPHASEFIELD, toleranceCrack);
    myIntegrationScheme.SetToleranceResidual(eDof::NONLOCALEQSTRAIN, toleranceNlEqStrain);


    Eigen::Matrix2d dispRHS;
    dispRHS(0, 0) = 0;
    dispRHS(1, 0) = simulationTime;
    dispRHS(0, 1) = 0;
    dispRHS(1, 1) = loadFactor;

    //    myIntegrationScheme.AddTimeDependentConstraint(loadId, dispRHS);
    myIntegrationScheme.SetTimeDependentLoadCase(loadId, dispRHS);

    structure.GetLogger() << "*********************************** \n"
                          << "**      solve                    ** \n"
                          << "*********************************** \n\n";


    myIntegrationScheme.Solve(simulationTime);
    structure.GetLogger() << "Total number of Dofs: \t" << structure.GetNumTotalDofs() << "\n\n";
}


void AssignSection(NuTo::StructureFeti& structure)
{
    auto section = NuTo::SectionPlane::Create(thickness, true);
    structure.ElementTotalSetSection(section);
}

void AssignMaterial(NuTo::StructureFeti& structure)
{

    NuTo::ConstitutiveBase* phaseField = new NuTo::PhaseField(youngsModulus, poissonsRatio, lengthScaleParameter,
                                                              fractureEnergy, artificialViscosity, energyDecomposition);

    int material00 = structure.AddConstitutiveLaw(phaseField);

    structure.ElementTotalSetConstitutiveLaw(material00);

    //        int material00 = structure.ConstitutiveLawCreate(eConstitutiveType::LINEAR_ELASTIC_ENGINEERING_STRESS);
    //
    //        structure.ConstitutiveLawSetParameterDouble(material00, eConstitutiveParameter::YOUNGS_MODULUS,
    //        youngsModulus);
    //        structure.ConstitutiveLawSetParameterDouble(material00, eConstitutiveParameter::POISSONS_RATIO,
    //        poissonsRatio);
    //
    //        structure.ElementTotalSetConstitutiveLaw(material00);

    //
    //    int material00 = structure.ConstitutiveLawCreate(eConstitutiveType::GRADIENT_DAMAGE_ENGINEERING_STRESS);
    //    structure.ConstitutiveLawSetParameterDouble(material00, eConstitutiveParameter::YOUNGS_MODULUS,
    //    youngsModulus);
    //    structure.ConstitutiveLawSetParameterDouble(material00, eConstitutiveParameter::POISSONS_RATIO,
    //    poissonsRatio);
    //    structure.ConstitutiveLawSetParameterDouble(material00, eConstitutiveParameter::TENSILE_STRENGTH,
    //    tensileStrength);
    //    structure.ConstitutiveLawSetParameterDouble(material00, eConstitutiveParameter::COMPRESSIVE_STRENGTH,
    //                                                compressiveStrength);
    //    structure.ConstitutiveLawSetParameterDouble(material00, eConstitutiveParameter::NONLOCAL_RADIUS,
    //    nonlocalRadius);
    //    structure.ConstitutiveLawSetParameterDouble(material00, eConstitutiveParameter::FRACTURE_ENERGY,
    //    fractureEnergy);
    //
    //    structure.ElementTotalSetConstitutiveLaw(material00);
}
