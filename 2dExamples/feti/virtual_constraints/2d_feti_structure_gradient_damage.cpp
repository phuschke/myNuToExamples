
#include <mpi.h>


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

//using EigenSolver = Eigen::SparseQR<Eigen::SparseMatrix<double>,Eigen::COLAMDOrdering<int>>;
//using EigenSolver = Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>>;
//using EigenSolver = Eigen::PardisoLU<Eigen::SparseMatrix<double>>;
using EigenSolver = Eigen::SparseLU<Eigen::SparseMatrix<double>,Eigen::COLAMDOrdering<int>>;


constexpr   int         dimension                   = 2;
constexpr   double      thickness                   = 1.0;

// material
constexpr   double      nonlocalRadius              = 5;                   // mm

constexpr   double      youngsModulus               = 4.0;                // N/mm^2
constexpr   double      poissonsRatio               = 0.2;
constexpr   double      fractureEnergy              = 1e-6;                   // N/mm
constexpr   double      compressiveStrength         = 30.e-4;                  // N/mm
constexpr   double      tensileStrength             = 3.e-4;                  // N/mm

//constexpr   double      youngsModulus               = 4.0e4;
//constexpr   double      poissonsRatio               = 0.2;
//constexpr   double      tensileStrength             = 3;
//constexpr   double      compressiveStrength         = 30;
//constexpr   double      fractureEnergy              = 0.1;

// integration
constexpr   bool        performLineSearch           = true;
constexpr   bool        automaticTimeStepping       = true;
constexpr   double      timeStep                    = 1e-1;
constexpr   double      minTimeStep                 = 1e-3;
constexpr   double      maxTimeStep                 =  1e-1;

constexpr   double      toleranceDisp              = 1e-6;
constexpr   double      toleranceNlEqStrain        = 1e-6;

constexpr   double      simulationTime              = 1.0;
constexpr   double      loadFactor                  = -2e-2;
constexpr   double      maxIterations              = 10;


const Eigen::Vector2d directionX    = Eigen::Vector2d::UnitX();
const Eigen::Vector2d directionY    = Eigen::Vector2d::UnitY();

void AssignSection(NuTo::StructureFeti& structure);
void AssignMaterial(NuTo::StructureFeti& structure);



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
    structure.InterpolationTypeAdd(interpolationTypeId, eDof::COORDINATES,      eTypeOrder::EQUIDISTANT1);
    structure.InterpolationTypeAdd(interpolationTypeId, eDof::DISPLACEMENTS,    eTypeOrder::EQUIDISTANT1);
    structure.InterpolationTypeAdd(interpolationTypeId, eDof::NONLOCALEQSTRAIN, eTypeOrder::EQUIDISTANT1);

    structure.ImportMeshJson(meshFile,interpolationTypeId);


    AssignMaterial(structure);
    AssignSection(structure);

    structure.GetLogger() << "**********************************************" << "\n";
    structure.GetLogger() << "**  virtual constraints                     **" << "\n";
    structure.GetLogger() << "**********************************************" << "\n\n";

    Eigen::VectorXd nodeCoords(2);


    if (structure.mRank == 1)
    {
        int groupNodesFakeConstraints00 = structure.GroupCreate(eGroupId::Nodes);
        int groupNodesFakeConstraints01 = structure.GroupCreate(eGroupId::Nodes);

        nodeCoords[0] = 10;
        nodeCoords[1] = 0;
        structure.GroupAddNodeRadiusRange(groupNodesFakeConstraints00, nodeCoords, 0, 1.e-6);
        structure.ConstraintLinearSetDisplacementNodeGroup(groupNodesFakeConstraints00, directionX, 0);


        nodeCoords[0] = 10;
        nodeCoords[1] = 10;
        structure.GroupAddNodeRadiusRange(groupNodesFakeConstraints01, nodeCoords, 0, 1.e-6);
        structure.ConstraintLinearSetDisplacementNodeGroup(groupNodesFakeConstraints01, directionX, 0);
        structure.ConstraintLinearSetDisplacementNodeGroup(groupNodesFakeConstraints01, directionY, 0);

        structure.GetLogger() << "Number of nodes that are constraint in 1st group: \t"
                              << structure.GroupGetNumMembers(groupNodesFakeConstraints00) << "\n";

        structure.GetLogger() << "Number of nodes that are constraint in 2nd group: \t"
                              << structure.GroupGetNumMembers(groupNodesFakeConstraints01) << "\n";


    }

    if (structure.mRank == 0)
    {
        int groupNodesFakeConstraints00 = structure.GroupCreate(eGroupId::Nodes);
        int groupNodesFakeConstraints01 = structure.GroupCreate(eGroupId::Nodes);

        nodeCoords[0] = 40;
        nodeCoords[1] = 0;
        structure.GroupAddNodeRadiusRange(groupNodesFakeConstraints00, nodeCoords, 0, 1.e-6);
        structure.ConstraintLinearSetDisplacementNodeGroup(groupNodesFakeConstraints00, directionX, 0);


        nodeCoords[0] = 50;
        nodeCoords[1] = 10;
        structure.GroupAddNodeRadiusRange(groupNodesFakeConstraints01, nodeCoords, 0, 1.e-6);
        structure.ConstraintLinearSetDisplacementNodeGroup(groupNodesFakeConstraints01, directionX, 0);
        structure.ConstraintLinearSetDisplacementNodeGroup(groupNodesFakeConstraints01, directionY, 0);

        structure.GetLogger() << "Number of nodes that are constraint in 1st group: \t"
                              << structure.GroupGetNumMembers(groupNodesFakeConstraints00) << "\n";

        structure.GetLogger() << "Number of nodes that are constraint in 2nd group: \t"
                              << structure.GroupGetNumMembers(groupNodesFakeConstraints01) << "\n";


    }

    if (structure.mRank == 2)
    {
        int groupNodesFakeConstraints00 = structure.GroupCreate(eGroupId::Nodes);
        int groupNodesFakeConstraints01 = structure.GroupCreate(eGroupId::Nodes);

        nodeCoords[0] = 40;
        nodeCoords[1] = 0;
        structure.GroupAddNodeRadiusRange(groupNodesFakeConstraints00, nodeCoords, 0, 1.e-6);
        structure.ConstraintLinearSetDisplacementNodeGroup(groupNodesFakeConstraints00, directionX, 0);


        nodeCoords[0] = 40;
        nodeCoords[1] = 10;
        structure.GroupAddNodeRadiusRange(groupNodesFakeConstraints01, nodeCoords, 0, 1.e-6);
        structure.ConstraintLinearSetDisplacementNodeGroup(groupNodesFakeConstraints01, directionX, 0);
        structure.ConstraintLinearSetDisplacementNodeGroup(groupNodesFakeConstraints01, directionY, 0);

        structure.GetLogger() << "Number of nodes that are constraint in 1st group: \t"
                              << structure.GroupGetNumMembers(groupNodesFakeConstraints00) << "\n";

        structure.GetLogger() << "Number of nodes that are constraint in 2nd group: \t"
                              << structure.GroupGetNumMembers(groupNodesFakeConstraints01) << "\n";


    }


    if (structure.mRank == 3)
    {
        int groupNodesFakeConstraints00 = structure.GroupCreate(eGroupId::Nodes);
        int groupNodesFakeConstraints01 = structure.GroupCreate(eGroupId::Nodes);

        nodeCoords[0] = 50;
        nodeCoords[1] = 0;
        structure.GroupAddNodeRadiusRange(groupNodesFakeConstraints00, nodeCoords, 0, 1.e-6);
        structure.ConstraintLinearSetDisplacementNodeGroup(groupNodesFakeConstraints00, directionX, 0);


        nodeCoords[0] = 50;
        nodeCoords[1] = 10;
        structure.GroupAddNodeRadiusRange(groupNodesFakeConstraints01, nodeCoords, 0, 1.e-6);
        structure.ConstraintLinearSetDisplacementNodeGroup(groupNodesFakeConstraints01, directionX, 0);
        structure.ConstraintLinearSetDisplacementNodeGroup(groupNodesFakeConstraints01, directionY, 0);

        structure.GetLogger() << "Number of nodes that are constraint in 1st group: \t"
                              << structure.GroupGetNumMembers(groupNodesFakeConstraints00) << "\n";

        structure.GetLogger() << "Number of nodes that are constraint in 2nd group: \t"
                              << structure.GroupGetNumMembers(groupNodesFakeConstraints01) << "\n";

    }

    structure.NodeBuildGlobalDofs(__PRETTY_FUNCTION__);

    structure.GetLogger() << "**********************************************" << "\n";
    structure.GetLogger() << "**  real constraints                        **" << "\n";
    structure.GetLogger() << "**********************************************" << "\n\n";


    if (rank == 1)
    {
        nodeCoords[0] = 0;
        nodeCoords[1] = 0;

        int groupNodesLeftBoundary = structure.GroupCreate(eGroupId::Nodes);

        structure.GroupAddNodeRadiusRange(groupNodesLeftBoundary, nodeCoords, 0, 1.e-6);
        structure.ApplyConstraintsTotalFeti(groupNodesLeftBoundary);
    }


    if (rank == 0)
    {

        nodeCoords[0] = 60;
        nodeCoords[1] = 0;

        int groupNodesRightBoundary = structure.GroupCreate(eGroupId::Nodes);

        structure.GroupAddNodeRadiusRange(groupNodesRightBoundary, nodeCoords, 0, 1.e-6);
        structure.ApplyConstraintsTotalFeti(groupNodesRightBoundary);


    }

    structure.GetLogger() << "**********************************************" << "\n";
    structure.GetLogger() << "**  load                                    **" << "\n";
    structure.GetLogger() << "**********************************************" << "\n\n";



    int loadNodeGroup = structure.GroupCreate(eGroupId::Nodes);

    nodeCoords[0] = 14;
    nodeCoords[1] = 10;
    structure.GroupAddNodeRadiusRange(loadNodeGroup, nodeCoords, 0, 1.e-6);

    structure.NodeInfo(10);
    // prescribe displacement of loadNodeGroup in Y direction
    std::map<int, double> dofIdAndPrescribedDisplacementMap;
    std::vector<int> nodeIds = structure.GroupGetMemberIds(loadNodeGroup);
    for (auto const& nodeId : nodeIds)
    {
        std::vector<int> dofIds = structure.NodeGetDofIds(nodeId, eDof::DISPLACEMENTS);
        dofIdAndPrescribedDisplacementMap.emplace(dofIds[1], -1.);
    }

    structure.ApplyPrescribedDisplacements(dofIdAndPrescribedDisplacementMap);

//    int loadId = structure.ConstraintLinearSetDisplacementNodeGroup(loadNodeGroup, directionY, 0);


    int loadId = structure.LoadCreateNodeGroupForce(0, loadNodeGroup, directionY, 0);


    structure.GetLogger() << "***********************************" << "\n";
    structure.GetLogger() << "**      Visualization            **" << "\n";
    structure.GetLogger() << "***********************************" << "\n\n";

    int groupAllElements = 9999;
    structure.GroupCreate(groupAllElements, eGroupId::Elements);
    structure.GroupAddElementsTotal(groupAllElements);
    structure.AddVisualizationComponent(groupAllElements, eVisualizeWhat::DISPLACEMENTS);
    structure.AddVisualizationComponent(groupAllElements, eVisualizeWhat::NONLOCAL_EQ_STRAIN);
    structure.AddVisualizationComponent(groupAllElements, eVisualizeWhat::ENGINEERING_STRAIN);
    structure.AddVisualizationComponent(groupAllElements, eVisualizeWhat::ENGINEERING_STRESS);
    structure.AddVisualizationComponent(groupAllElements, eVisualizeWhat::DAMAGE);


    structure.GetLogger() << "**********************************************" << "\n";
    structure.GetLogger() << "**  integration sheme                       **" << "\n";
    structure.GetLogger() << "**********************************************" << "\n\n";



    NuTo::NewmarkFeti<EigenSolver> myIntegrationScheme(&structure);
    boost::filesystem::path resultPath(std::string("/home/phuschke/results/feti/" + std::to_string(structure.mRank)));

    myIntegrationScheme.SetTimeStep                 ( timeStep                  );
    myIntegrationScheme.SetMaxNumIterations         ( maxIterations            );
    myIntegrationScheme.SetMinTimeStep              ( minTimeStep               );
    myIntegrationScheme.SetMaxTimeStep              ( maxTimeStep               );
    myIntegrationScheme.SetAutomaticTimeStepping    ( automaticTimeStepping     );
    myIntegrationScheme.SetResultDirectory          ( resultPath.string(), true );
    myIntegrationScheme.SetPerformLineSearch        ( performLineSearch         );
    myIntegrationScheme.SetToleranceResidual        ( eDof::DISPLACEMENTS,      toleranceDisp );
    myIntegrationScheme.SetToleranceResidual        ( eDof::NONLOCALEQSTRAIN,   toleranceNlEqStrain );

    Eigen::Matrix2d dispRHS;
    dispRHS(0, 0) = 0;
    dispRHS(1, 0) = simulationTime;
    dispRHS(0, 1) = 0;
    dispRHS(1, 1) = loadFactor;


    if (rank == 1)
    {
        myIntegrationScheme.AddResultGroupNodeForce("myforcebla", loadNodeGroup);

        nodeCoords[0] = 20;
        nodeCoords[1] = 10;
        int grpNodes_output_disp = structure.GroupCreate(eGroupId::Nodes);
        structure.GroupAddNodeRadiusRange(grpNodes_output_disp, nodeCoords, 0, 1.e-6);
        myIntegrationScheme.AddResultNodeDisplacements("mydisplacements", structure.GroupGetMemberIds(grpNodes_output_disp)[0]);
    }


//    myIntegrationScheme.AddTimeDependentConstraint(loadId, dispRHS);
    myIntegrationScheme.SetTimeDependentLoadCase(loadId, dispRHS);

    structure.GetLogger() << "***********************************" << "\n";
    structure.GetLogger() << "**      Solve                    **" << "\n";
    structure.GetLogger() << "***********************************" << "\n\n";

    structure.GetLogger()   << "Total number of Dofs: \t"
                            << structure.GetNumTotalDofs() << "\n\n";

    myIntegrationScheme.Solve(simulationTime);


    MPI_Finalize();
}


void AssignSection(NuTo::StructureFeti& structure)
{
    auto section = NuTo::SectionPlane::Create(thickness,true);
    structure.ElementTotalSetSection(section);
}

void AssignMaterial(NuTo::StructureFeti& structure)
{
    structure.GetLogger() << "***********************************" << "\n";
    structure.GetLogger() << "**      Material                 **" << "\n";
    structure.GetLogger() << "***********************************" << "\n\n";

    int material00 = structure.ConstitutiveLawCreate(eConstitutiveType::GRADIENT_DAMAGE_ENGINEERING_STRESS);
    structure.ConstitutiveLawSetParameterDouble(material00, eConstitutiveParameter::YOUNGS_MODULUS, youngsModulus);
    structure.ConstitutiveLawSetParameterDouble(material00, eConstitutiveParameter::POISSONS_RATIO, poissonsRatio);
    structure.ConstitutiveLawSetParameterDouble(material00, eConstitutiveParameter::TENSILE_STRENGTH, tensileStrength);
    structure.ConstitutiveLawSetParameterDouble(material00, eConstitutiveParameter::COMPRESSIVE_STRENGTH, compressiveStrength);
    structure.ConstitutiveLawSetParameterDouble(material00, eConstitutiveParameter::NONLOCAL_RADIUS, nonlocalRadius);
    structure.ConstitutiveLawSetParameterDouble(material00, eConstitutiveParameter::FRACTURE_ENERGY, fractureEnergy);

    structure.ElementTotalSetConstitutiveLaw(material00);



}
