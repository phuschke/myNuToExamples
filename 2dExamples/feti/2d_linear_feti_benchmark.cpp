
#include <mpi.h>
#include <boost/mpi.hpp>
#include <iostream>
#include <string.h>
#include <vector>
#include "nuto/mechanics/structures/unstructured/StructureFETI.h"
#include "nuto/mechanics/structures/unstructured/Structure.h"
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/SparseQR>

#include <ctime>
#include <cstdlib>
#include <chrono>

#include "nuto/mechanics/timeIntegration/NewmarkFeti.h"
#include "../myNutoExamples/EnumsAndTypedefs.h"
#include "nuto/mechanics/constitutive/laws/PhaseField.h"
#include "nuto/mechanics/nodes/NodeBase.h"

#include "nuto/math/SparseMatrixCSR.h"
#include "nuto/math/SparseMatrixCSRGeneral.h"
#include "nuto/mechanics/dofSubMatrixStorage/BlockSparseMatrix.h"
#include "nuto/mechanics/structures/StructureOutputBlockMatrix.h"
#include "boost/filesystem.hpp"

constexpr int dim = 2;
using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::cout;
using std::endl;




constexpr   int         dimension                   = 2;
constexpr   double      thickness                   = 1.0;

// material
constexpr   double      youngsModulus               = 2.1e5;                // N/mm^2
constexpr   double      poissonsRatio               = 0.3;
constexpr   double      nonlocalRadius              = 1.0e-1;                   // mm
constexpr   double      fractureEnergy              = 100.0;                   // N/mm
constexpr   double      compressiveStrength         = 30000;                  // N/mm
constexpr   double      tensileStrength             = 3000;                  // N/mm


// integration
constexpr   double      timeStep                    = 1.0;
constexpr   double      toleranceDisp               = 1e-6;
constexpr   double      simulationTime              = 1.0;
constexpr   double      loadFactor                  = 1.0;


const NuTo::FullVector<double, dimension> directionX    = NuTo::FullVector<double, dimension>::UnitX();
const NuTo::FullVector<double, dimension> directionY    = NuTo::FullVector<double, dimension>::UnitY();

void Feti();

void SingleDomain();

int main(int argc, char* argv[])
{
    boost::mpi::environment     env;
    boost::mpi::communicator    world;

//    Feti();


    if(world.rank() == 0)
        SingleDomain();

}


void Feti()
{
    boost::mpi::communicator    world;

    const int rank = world.rank();

    NuTo::StructureFETI structure(dimension);
    structure.SetShowTime(false);

    std::string meshFile = "feti_beam_fine_tri.msh_" + std::to_string(rank);

    const int interpolationTypeId = structure.InterpolationTypeCreate(eShapeType::TRIANGLE2D);
    structure.InterpolationTypeAdd(interpolationTypeId, eDof::COORDINATES,     eTypeOrder::EQUIDISTANT1);
    structure.InterpolationTypeAdd(interpolationTypeId, eDof::DISPLACEMENTS,   eTypeOrder::EQUIDISTANT1);

    structure.ImportMeshJson(meshFile,interpolationTypeId);

    // section
    int sectionId = structure.SectionCreate(eSectionType::PLANE_STRESS);
    structure.SectionSetThickness(sectionId, thickness);
    structure.ElementTotalSetSection(sectionId);

    // material
    int materialId = structure.ConstitutiveLawCreate(eConstitutiveType::LINEAR_ELASTIC_ENGINEERING_STRESS);
    structure.ConstitutiveLawSetParameterDouble(materialId, eConstitutiveParameter::YOUNGS_MODULUS, youngsModulus);
    structure.ConstitutiveLawSetParameterDouble(materialId, eConstitutiveParameter::POISSONS_RATIO, poissonsRatio);
    structure.ElementTotalSetConstitutiveLaw(materialId);

    // constraints
    Eigen::VectorXd nodeCoords(2);

    int groupNodesLeftBoundary = structure.GroupCreate(eGroupId::Nodes);
    nodeCoords[0] = 0;
    nodeCoords[1] = 0;
    structure.GroupAddNodeRadiusRange(groupNodesLeftBoundary, nodeCoords, 0, 1.e-6);
    structure.ConstraintLinearSetDisplacementNodeGroup(groupNodesLeftBoundary, directionX, 0);
    structure.ConstraintLinearSetDisplacementNodeGroup(groupNodesLeftBoundary, directionY, 0);

    int groupNodesRightBoundary = structure.GroupCreate(eGroupId::Nodes);
    nodeCoords[0] = 60;
    nodeCoords[1] = 0;
    structure.GroupAddNodeRadiusRange(groupNodesRightBoundary, nodeCoords, 0, 1.e-6);
    structure.ConstraintLinearSetDisplacementNodeGroup(groupNodesRightBoundary, directionX, 0);
    structure.ConstraintLinearSetDisplacementNodeGroup(groupNodesRightBoundary, directionY, 0);

    // loads
    int loadNodeGroup = structure.GroupCreate(eGroupId::Nodes);

    nodeCoords[0] = 30;
    nodeCoords[1] = 10;
    structure.GroupAddNodeRadiusRange(loadNodeGroup, nodeCoords, 0, 1.e-6);

    int loadId = structure.ConstraintLinearSetDisplacementNodeGroup(loadNodeGroup, directionY, 0);

    // visualization
    int groupAllElements = 9999;
    structure.GroupCreate(groupAllElements, eGroupId::Elements);
    structure.GroupAddElementsTotal(groupAllElements);
    structure.AddVisualizationComponent(groupAllElements, eVisualizeWhat::DISPLACEMENTS);

    // time integration
    NuTo::NewmarkFeti timeIntegration(&structure);
    boost::filesystem::path resultPath("/home/phuschke/results/feti/benchmark/" + std::to_string(structure.mRank));

    timeIntegration.SetTimeStep                 ( timeStep                  );
    timeIntegration.SetResultDirectory          ( resultPath.string(), true );
    timeIntegration.SetToleranceResidual        ( eDof::DISPLACEMENTS, toleranceDisp );

    Eigen::Matrix2d dispRHS;
    dispRHS(0, 0) = 0;
    dispRHS(1, 0) = simulationTime;
    dispRHS(0, 1) = 0;
    dispRHS(1, 1) = loadFactor;

    timeIntegration.AddTimeDependentConstraint(loadId, dispRHS);
    timeIntegration.Solve(simulationTime);
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void SingleDomain()
{

    NuTo::Structure structure(dimension);
    structure.SetShowTime(false);
    structure.SetNumProcessors(1);

    std::string meshFile = "linear_benchmark.msh";

    const int interpolationTypeId = structure.InterpolationTypeCreate(eShapeType::TRIANGLE2D);
    structure.InterpolationTypeAdd(interpolationTypeId, eDof::COORDINATES,     eTypeOrder::EQUIDISTANT1);
    structure.InterpolationTypeAdd(interpolationTypeId, eDof::DISPLACEMENTS,   eTypeOrder::EQUIDISTANT1);
    structure.ElementTotalSetInterpolationType(interpolationTypeId);

    cout << "Start import" << endl;

    structure.ImportFromGmsh(meshFile);

    cout << "End import" << endl;


    cout << "Start set" << endl;
    structure.ElementTotalSetInterpolationType(interpolationTypeId);
    cout << "End set" << endl;

    cout << "Start convert" << endl;
    structure.ElementTotalConvertToInterpolationType();

    cout << "End concert" << endl;
    // section
    int sectionId = structure.SectionCreate(eSectionType::PLANE_STRESS);
    structure.SectionSetThickness(sectionId, thickness);
    structure.ElementTotalSetSection(sectionId);

    // material
    int materialId = structure.ConstitutiveLawCreate(eConstitutiveType::LINEAR_ELASTIC_ENGINEERING_STRESS);
    structure.ConstitutiveLawSetParameterDouble(materialId, eConstitutiveParameter::YOUNGS_MODULUS, youngsModulus);
    structure.ConstitutiveLawSetParameterDouble(materialId, eConstitutiveParameter::POISSONS_RATIO, poissonsRatio);
    structure.ElementTotalSetConstitutiveLaw(materialId);

    // constraints
    Eigen::VectorXd nodeCoords(2);

    int groupNodesLeftBoundary = structure.GroupCreate(eGroupId::Nodes);
    nodeCoords[0] = 0;
    nodeCoords[1] = 0;
    structure.GroupAddNodeRadiusRange(groupNodesLeftBoundary, nodeCoords, 0, 1.e-6);
    structure.ConstraintLinearSetDisplacementNodeGroup(groupNodesLeftBoundary, directionX, 0);
    structure.ConstraintLinearSetDisplacementNodeGroup(groupNodesLeftBoundary, directionY, 0);

    int groupNodesRightBoundary = structure.GroupCreate(eGroupId::Nodes);
    nodeCoords[0] = 60;
    nodeCoords[1] = 0;
    structure.GroupAddNodeRadiusRange(groupNodesRightBoundary, nodeCoords, 0, 1.e-6);
    structure.ConstraintLinearSetDisplacementNodeGroup(groupNodesRightBoundary, directionX, 0);
    structure.ConstraintLinearSetDisplacementNodeGroup(groupNodesRightBoundary, directionY, 0);

    // loads
    int loadNodeGroup = structure.GroupCreate(eGroupId::Nodes);

    nodeCoords[0] = 30;
    nodeCoords[1] = 10;
    structure.GroupAddNodeRadiusRange(loadNodeGroup, nodeCoords, 0, 1.e-6);

    int loadId = structure.ConstraintLinearSetDisplacementNodeGroup(loadNodeGroup, directionY, 0);

    // visualization
    int groupAllElements = 9999;
    structure.GroupCreate(groupAllElements, eGroupId::Elements);
    structure.GroupAddElementsTotal(groupAllElements);
    structure.AddVisualizationComponent(groupAllElements, eVisualizeWhat::DISPLACEMENTS);

    // time integration
    NuTo::NewmarkDirect timeIntegration(&structure);
    boost::filesystem::path resultPath(boost::filesystem::initial_path().string() + std::string("/resultFeti_reference"));

    timeIntegration.SetTimeStep                 ( timeStep                  );
    timeIntegration.SetResultDirectory          ( resultPath.string(), true );
    timeIntegration.SetToleranceResidual        ( eDof::DISPLACEMENTS, toleranceDisp );

    Eigen::Matrix2d dispRHS;
    dispRHS(0, 0) = 0;
    dispRHS(1, 0) = simulationTime;
    dispRHS(0, 1) = 0;
    dispRHS(1, 1) = loadFactor;

    timeIntegration.AddTimeDependentConstraint(loadId, dispRHS);

    cout << "Start solve: num dofs = " << structure.GetNumNodes() << endl;
    NuTo::Timer timer("Single domain solve");
    timeIntegration.Solve(simulationTime);
    timer.Reset();
    cout << "End solve" << endl;

}
