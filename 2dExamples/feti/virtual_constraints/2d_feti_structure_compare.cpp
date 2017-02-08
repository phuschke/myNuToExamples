
#include <mpi.h>
#include <boost/mpi.hpp>

#include <iostream>
#include <string.h>
#include <vector>
#include "mechanics/feti/StructureFeti.h"
#include "mechanics/structures/unstructured/Structure.h"
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/SparseQR>
#include "ConjugateProjectedGradient.h"
#include <ctime>
#include <cstdlib>
#include <chrono>
#include "ImportMesh.h"
#include "mechanics/timeIntegration/NewmarkFeti.h"
#include "../../../EnumsAndTypedefs.h"

#include "mechanics/nodes/NodeBase.h"

#include "math/SparseMatrixCSR.h"
#include "math/SparseMatrixCSRGeneral.h"
#include "mechanics/dofSubMatrixStorage/BlockSparseMatrix.h"
#include "mechanics/structures/StructureOutputBlockMatrix.h"
#include "boost/filesystem.hpp"

constexpr int dim = 2;
using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::cout;
using std::endl;




constexpr   int         dimension                   = 2;
constexpr   double      thickness                   = 1.0;

// material
constexpr   double      youngsModulus               = 4.0e4;
constexpr   double      poissonsRatio               = 0.2;
constexpr   double      tensileStrength             = 3;
constexpr   double      compressiveStrength         = 30;
constexpr   double      fractureEnergy              = 0.1;

// integration
constexpr   bool        performLineSearch           = true;
constexpr   bool        automaticTimeStepping       = true;
constexpr   double      timeStep                    = 1e-0;
constexpr   double      minTimeStep                 = 1e-5;
constexpr   double      maxTimeStep                 =  1e-0;
constexpr   double      toleranceDisp              = 1e-6;
constexpr   double      simulationTime              = 1.0;
constexpr   double      loadFactor                  = -0.01;
constexpr   double      maxInterations              = 10;


const NuTo::FullVector<double, dimension> directionX    = NuTo::FullVector<double, dimension>::UnitX();
const NuTo::FullVector<double, dimension> directionY    = NuTo::FullVector<double, dimension>::UnitY();

void AssignSection(NuTo::StructureFETI& rStructure)
{
    std::cout << "***********************************" << std::endl;
    std::cout << "**      Section rank = " << rStructure.mRank << std::endl;
    std::cout << "***********************************" << std::endl;

    int section00 = rStructure.SectionCreate(NuTo::eSectionType::PLANE_STRESS);
    rStructure.SectionSetThickness(section00, thickness);

    rStructure.ElementTotalSetSection(section00);
}

void AssignMaterial(NuTo::StructureFETI& rStructure)
{
    std::cout << "***********************************" << std::endl;
    std::cout << "**      Material rank = " << rStructure.mRank << std::endl;
    std::cout << "***********************************" << std::endl;

    int material00 = rStructure.ConstitutiveLawCreate(eConstitutiveType::LINEAR_ELASTIC_ENGINEERING_STRESS);

    rStructure.ConstitutiveLawSetParameterDouble(material00, eConstitutiveParameter::YOUNGS_MODULUS, youngsModulus);
    rStructure.ConstitutiveLawSetParameterDouble(material00, eConstitutiveParameter::POISSONS_RATIO, poissonsRatio);


//    int material00 = rStructure.ConstitutiveLawCreate(eConstitutiveType::LOCAL_DAMAGE_MODEL);
//    rStructure.ConstitutiveLawSetParameterDouble(material00, eConstitutiveParameter::YOUNGS_MODULUS,       youngsModulus);
//    rStructure.ConstitutiveLawSetParameterDouble(material00, eConstitutiveParameter::POISSONS_RATIO,       poissonsRatio);
//    rStructure.ConstitutiveLawSetParameterDouble(material00, eConstitutiveParameter::TENSILE_STRENGTH,     tensileStrength);
//    rStructure.ConstitutiveLawSetParameterDouble(material00, eConstitutiveParameter::COMPRESSIVE_STRENGTH, compressiveStrength);
//    rStructure.ConstitutiveLawSetParameterDouble(material00, eConstitutiveParameter::FRACTURE_ENERGY,      fractureEnergy);

    rStructure.ElementTotalSetConstitutiveLaw(material00);

}


int main(int argc, char* argv[])
{
    boost::mpi::environment  env(argc, argv);
    boost::mpi::communicator world;

    const int rank = world.rank();

    std::string meshFile = "feti_beam_fine_2_subdomains.mesh" + std::to_string(rank);
    std::cout << meshFile << std::endl;

    NuTo::StructureFETI structure(dim);
    structure.SetNumTimeDerivatives(0);

    const int interpolationTypeId = structure.InterpolationTypeCreate(eShapeType::QUAD2D);
    structure.InterpolationTypeAdd(interpolationTypeId, eDof::COORDINATES,     eTypeOrder::EQUIDISTANT1);
    structure.InterpolationTypeAdd(interpolationTypeId, eDof::DISPLACEMENTS,   eTypeOrder::EQUIDISTANT1);


    structure.ImportMeshJson(meshFile,interpolationTypeId);

    structure.SetVerboseLevel(10);
    structure.SetShowTime(false);

    
    AssignMaterial(structure);
    AssignSection(structure);

    cout << "**********************************************" << endl;
    cout << "**  virtual constraints                     **" << endl;
    cout << "**********************************************" << endl;

    Eigen::VectorXd nodeCoords(2);


    // All subdomains are floating in total FETI
    structure.mNumRigidBodyModes = 3;

    if (structure.mRank == 0)
    {
        int groupNodesFakeConstraints00 = structure.GroupCreate(eGroupId::Nodes);
        int groupNodesFakeConstraints01 = structure.GroupCreate(eGroupId::Nodes);

        nodeCoords[0] = 28;
        nodeCoords[1] = 0;
        structure.GroupAddNodeRadiusRange(groupNodesFakeConstraints00, nodeCoords, 0, 1.e-6);
        structure.ConstraintLinearSetDisplacementNodeGroup(groupNodesFakeConstraints00, directionX, 0);


        nodeCoords[0] = 28;
        nodeCoords[1] = 10;
        structure.GroupAddNodeRadiusRange(groupNodesFakeConstraints01, nodeCoords, 0, 1.e-6);
        structure.ConstraintLinearSetDisplacementNodeGroup(groupNodesFakeConstraints01, directionX, 0);
        structure.ConstraintLinearSetDisplacementNodeGroup(groupNodesFakeConstraints01, directionY, 0);

    }

    if (structure.mRank == 1)
    {
        int groupNodesFakeConstraints00 = structure.GroupCreate(eGroupId::Nodes);
        int groupNodesFakeConstraints01 = structure.GroupCreate(eGroupId::Nodes);

        nodeCoords[0] = 32;
        nodeCoords[1] = 0;
        structure.GroupAddNodeRadiusRange(groupNodesFakeConstraints00, nodeCoords, 0, 1.e-6);
        structure.ConstraintLinearSetDisplacementNodeGroup(groupNodesFakeConstraints00, directionX, 0);


        nodeCoords[0] = 58;
        nodeCoords[1] = 10;
        structure.GroupAddNodeRadiusRange(groupNodesFakeConstraints01, nodeCoords, 0, 1.e-6);
        structure.ConstraintLinearSetDisplacementNodeGroup(groupNodesFakeConstraints01, directionX, 0);
        structure.ConstraintLinearSetDisplacementNodeGroup(groupNodesFakeConstraints01, directionY, 0);

    }

    structure.NodeBuildGlobalDofs(__PRETTY_FUNCTION__);

    cout << "**********************************************" << endl;
    cout << "**  real constraints                        **" << endl;
    cout << "**********************************************" << endl;

    if (structure.mRank == 0)
        structure.NodeInfo(10);

    int groupNodesLeftBoundary = structure.GroupCreate(eGroupId::Nodes);

    structure.GroupAddNodeCoordinateRange(groupNodesLeftBoundary,0,-1.e-6,+1.e-6);
    Eigen::VectorXi boundaryNodes = structure.GroupGetMemberIds(groupNodesLeftBoundary);

    for (int i = 0; i < boundaryNodes.rows(); ++i)
    {
        const int nodeId = boundaryNodes[i];
        std::vector<int> dofIds = structure.NodeGetDofIds(nodeId, eDof::DISPLACEMENTS);

        for (const auto& id : dofIds)
            structure.mBoundaryDofIds.push_back(id);
    }

    for (const auto& id : structure.mBoundaryDofIds)
        std::cout << "structure.mBoundaryDofIds \t" << id << std::endl;

    cout << "***************************************************************" << endl;
    cout << "**  determine the global ids for the constraints             **" << endl;
    cout << "***************************************************************" << endl;

    const int numProcesses = world.size();
    // recvCount:
    // Contais the number of elements that are received from each process.
    std::vector<int> recvCount;
    recvCount.resize(numProcesses, 0);

    boost::mpi::all_gather<int>(world,structure.mBoundaryDofIds.size(),recvCount);

    // displs:
    // Entry i specifies the displacement (relative to recvbuf) at which to place the incoming data from process i.
    std::vector<int> displs;
    displs.resize(numProcesses, 0);
    for (int i = 1; i < numProcesses; ++i)
        displs[i] = displs[i-1] + recvCount[i-1];


    if (structure.mRank == 0)
    {
        for (const auto& i : recvCount)
            std::cout << "recvCount \t" << i << std::endl;

        for (const auto& i : displs)
            std::cout << "displs \t" << i << std::endl;
    }

    const int numLocalBoundaryDofIds = structure.mBoundaryDofIds.size();
    MPI_Allreduce(&numLocalBoundaryDofIds,
                  &structure.mNumTotalBoundaryDofIds,
                  1,
                  MPI_INT,
                  MPI_SUM,
                  MPI_COMM_WORLD);


    std::cout << "structure.mBoundaryDofIds.size() \t" << structure.mBoundaryDofIds.size() << std::endl;
    std::cout << "structure.mNumTotalBoundaryDofIds \t" << structure.mNumTotalBoundaryDofIds << std::endl;


    structure.mGlobalStartIndexBoundaryDofIds = displs[structure.mRank];

    cout << "**********************************************" << endl;
    cout << "**  load                                    **" << endl;
    cout << "**********************************************" << endl;



    int loadNodeGroup = structure.GroupCreate(eGroupId::Nodes);

    nodeCoords[0] = 16;
    nodeCoords[1] = 10;
    structure.GroupAddNodeRadiusRange(loadNodeGroup, nodeCoords, 0, 1.e-6);

//    int loadId = structure.ConstraintLinearSetDisplacementNodeGroup(loadNodeGroup, directionY, 0);
    int loadId = structure.LoadCreateNodeGroupForce(0, loadNodeGroup, directionY, 100000);


    std::cout << "***********************************" << std::endl;
    std::cout << "**      Visualization rank = " << structure.mRank        << std::endl;
    std::cout << "***********************************" << std::endl;
    int groupAllElements = 9999;
    structure.GroupCreate(groupAllElements, eGroupId::Elements);
    structure.GroupAddElementsTotal(groupAllElements);
    structure.AddVisualizationComponent(groupAllElements, eVisualizeWhat::DISPLACEMENTS);
    structure.AddVisualizationComponent(groupAllElements, eVisualizeWhat::ENGINEERING_STRAIN);
    structure.AddVisualizationComponent(groupAllElements, eVisualizeWhat::ENGINEERING_STRESS);
//    structure.AddVisualizationComponent(groupAllElements, eVisualizeWhat::DAMAGE);


    cout << "**********************************************" << endl;
    cout << "**  integration sheme                       **" << endl;
    cout << "**********************************************" << endl;


    NuTo::NewmarkFeti myIntegrationScheme(&structure);
    boost::filesystem::path resultPath(std::string("/home/phuschke/results/feti/" + std::to_string(structure.mRank)));

    myIntegrationScheme.SetTimeStep                 ( timeStep                  );
    myIntegrationScheme.SetMaxNumIterations         ( maxInterations            );
    myIntegrationScheme.SetMinTimeStep              ( minTimeStep               );
    myIntegrationScheme.SetMaxTimeStep              ( maxTimeStep               );
    myIntegrationScheme.SetAutomaticTimeStepping    ( automaticTimeStepping     );
    myIntegrationScheme.SetResultDirectory          ( resultPath.string(), true );
    myIntegrationScheme.SetPerformLineSearch        ( performLineSearch         );
    myIntegrationScheme.SetToleranceResidual        ( eDof::DISPLACEMENTS, toleranceDisp );

    Eigen::Matrix2d dispRHS;
    dispRHS(0, 0) = 0;
    dispRHS(1, 0) = simulationTime;
    dispRHS(0, 1) = 0;
    dispRHS(1, 1) = loadFactor;

//    myIntegrationScheme.AddTimeDependentConstraint(loadId, dispRHS);
    myIntegrationScheme.SetTimeDependentLoadCase(loadId, dispRHS);

    cout << "***********************************" << std::endl;
    cout << "**      Solve                    **" << std::endl;
    cout << "***********************************" << std::endl;

    myIntegrationScheme.Solve(simulationTime);

}
