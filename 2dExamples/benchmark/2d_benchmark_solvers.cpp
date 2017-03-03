

#include <iostream>
#include <vector>
#include "mechanics/feti/StructureFeti.h"
#include <ctime>
#include <chrono>
#include "mechanics/timeIntegration/NewmarkDirect.h"
#include "../../EnumsAndTypedefs.h"

#include "mechanics/nodes/NodeBase.h"

#include "math/SparseMatrixCSR.h"
#include "math/SparseMatrixCSRGeneral.h"
#include "boost/filesystem.hpp"

#include "mechanics/mesh/MeshGenerator.h"
#include "mechanics/dofSubMatrixSolvers/SolverMUMPS.h"
#include "mechanics/dofSubMatrixSolvers/SolverEigen.h"
#include "base/Timer.h"

constexpr   int         dim                         = 2;
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
constexpr   double      minTimeStep                 = 1e-3;
constexpr   double      maxTimeStep                 =  1e-0;
constexpr   double      toleranceDisp              = 1e-6;
constexpr   double      simulationTime              = 1.0;
constexpr   double      loadFactor                  = -2;
constexpr   double      maxInterations              = 10;

void AssignSection(NuTo::Structure& structure);
void AssignMaterial(NuTo::Structure& structure);
int BuildStructure(NuTo::Structure& structure, int numElementsX, int numElementsY);
double Solve(NuTo::NewmarkDirect& newmarkDirect);
void Benchmark(std::unique_ptr<NuTo::SolverBase> solver, int i);

int main(int argc, char* argv[])
{

    std::ofstream file("benchmark.txt");
    file.close();

    for (int i = 10; i < 257; i+=10)
    {
        Benchmark(std::make_unique<NuTo::SolverEigen<Eigen::SparseLU<Eigen::SparseMatrix<double>,Eigen::COLAMDOrdering<int>>>>(), i);
        Benchmark(std::make_unique<NuTo::SolverEigen<Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>>>>(), i);
        Benchmark(std::make_unique<NuTo::SolverMUMPS>(), i);
    }


}


void Benchmark(std::unique_ptr<NuTo::SolverBase> solver, int i)
{
    int numElementsX = 10*i;
    int numElementsY = 100*i;

    NuTo::Structure structure(dim);
    int loadId = BuildStructure(structure, numElementsX, numElementsY);

    structure.GetLogger() << "**********************************************" << "\n";
    structure.GetLogger() << "**  integration scheme                       **" << "\n";
    structure.GetLogger() << "**********************************************" << "\n\n";

    NuTo::NewmarkDirect timeIntegration(&structure);
    timeIntegration.SetTimeStep(timeStep);
    timeIntegration.SetResultDirectory(boost::filesystem::initial_path().string());

    Eigen::Matrix2d dispRHS;
    dispRHS(0, 0) = 0;
    dispRHS(1, 0) = simulationTime;
    dispRHS(0, 1) = 0;
    dispRHS(1, 1) = loadFactor;

    timeIntegration.AddTimeDependentConstraint(loadId, dispRHS);



    structure.GetLogger() << "***********************************" << "\n";
    structure.GetLogger() << "**      Solve                    **" << "\n";
    structure.GetLogger() << "***********************************" << "\n\n";

    std::vector<double> timeToSolve;

    timeIntegration.SetSolver(std::move(solver));
    timeToSolve.push_back(Solve(timeIntegration));



    std::ofstream file("benchmark.txt", std::ofstream::app);
    for (const auto& time : timeToSolve)
    {
        file << structure.GetNumTotalDofs() << "\t" << time << "\t" << structure.GetNumElements() <<  "\n";
    }

    file.close();

}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double Solve(NuTo::NewmarkDirect& timeIntegration)
{
    try
    {
        NuTo::Timer timer("Time for NewmarkDirect");
        timeIntegration.Solve(simulationTime);
        return timer.GetTimeDifference();

    }
    catch(NuTo::MechanicsException e)
    {
        std::cout << e.ErrorMessage() << std::endl;
        return EXIT_FAILURE;
    }
    catch(NuTo::MathException e)
    {
        std::cout << e.ErrorMessage() << std::endl;
        return EXIT_FAILURE;
    }
    catch(...)
    {
        std::cout << "Some other exception" << std::endl;
        return EXIT_FAILURE;
    }

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int BuildStructure(NuTo::Structure& structure, int numElementsX, int numElementsY)
{

    std::vector<double> meshDimensions;
    meshDimensions.push_back(100.);
    meshDimensions.push_back(10.);

    std::vector<int> numCells;
    numCells.push_back(numElementsX);
    numCells.push_back(numElementsY);

    auto importContainer = NuTo::MeshGenerator::Grid(structure, meshDimensions, numCells, NuTo::Interpolation::eShapeType::QUAD2D);
    const int interpolationTypeId = importContainer.second;
    structure.InterpolationTypeAdd(interpolationTypeId, eDof::DISPLACEMENTS,   eTypeOrder::EQUIDISTANT1);
    structure.ElementTotalConvertToInterpolationType();

    AssignMaterial(structure);
    AssignSection(structure);

    structure.GetLogger() << "**********************************************" << "\n";
    structure.GetLogger() << "**  constraints                             **" << "\n";
    structure.GetLogger() << "**********************************************" << "\n\n";

    Eigen::VectorXd nodeCoords(2);
    nodeCoords[0] = 0;
    nodeCoords[1] = 0;

    int groupNodesLeftBoundary = structure.GroupCreate(eGroupId::Nodes);

    structure.GroupAddNodeRadiusRange(groupNodesLeftBoundary, nodeCoords, 0, 1.e-6);

    structure.ConstraintLinearSetDisplacementNodeGroup(groupNodesLeftBoundary, Vector2d::UnitX(), 0.0);
    structure.ConstraintLinearSetDisplacementNodeGroup(groupNodesLeftBoundary, Vector2d::UnitY(), 0.0);

    nodeCoords[0] = 100;
    nodeCoords[1] = 0;

    int groupNodesRightBoundary = structure.GroupCreate(eGroupId::Nodes);
    structure.GroupAddNodeRadiusRange(groupNodesRightBoundary, nodeCoords, 0, 1.e-6);

    structure.ConstraintLinearSetDisplacementNodeGroup(groupNodesRightBoundary, Vector2d::UnitX(), 0.0);
    structure.ConstraintLinearSetDisplacementNodeGroup(groupNodesRightBoundary, Vector2d::UnitY(), 0.0);

    structure.GetLogger() << "**********************************************" << "\n";
    structure.GetLogger() << "**  load                                    **" << "\n";
    structure.GetLogger() << "**********************************************" << "\n\n";

    int loadNodeGroup = structure.GroupCreate(eGroupId::Nodes);
    nodeCoords[0] = 20;
    nodeCoords[1] = 10;
    structure.GroupAddNodeRadiusRange(loadNodeGroup, nodeCoords, 0, 1.e-6);
    return structure.ConstraintLinearSetDisplacementNodeGroup(loadNodeGroup, Vector2d::UnitY(), 1);



}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void AssignSection(NuTo::Structure& structure)
{
    structure.GetLogger() << "***********************************" << "\n";
    structure.GetLogger() << "**      Section                  **" << "\n";
    structure.GetLogger() << "***********************************" << "\n\n";

    int section00 = structure.SectionCreate(NuTo::eSectionType::PLANE_STRESS);
    structure.SectionSetThickness(section00, thickness);

    structure.ElementTotalSetSection(section00);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void AssignMaterial(NuTo::Structure& structure)
{
    structure.GetLogger() << "***********************************" << "\n";
    structure.GetLogger() << "**      Material                 **" << "\n";
    structure.GetLogger() << "***********************************" << "\n\n";

    int material00 = structure.ConstitutiveLawCreate(eConstitutiveType::LINEAR_ELASTIC_ENGINEERING_STRESS);

    structure.ConstitutiveLawSetParameterDouble(material00, eConstitutiveParameter::YOUNGS_MODULUS, youngsModulus);
    structure.ConstitutiveLawSetParameterDouble(material00, eConstitutiveParameter::POISSONS_RATIO, poissonsRatio);

    structure.ElementTotalSetConstitutiveLaw(material00);

}
