

#include <iostream>
#include <vector>
#include "mechanics/feti/StructureFeti.h"
#include <ctime>
#include <chrono>
#include "mechanics/timeIntegration/NewmarkDirect.h"
#include "../../../EnumsAndTypedefs.h"

#include "mechanics/nodes/NodeBase.h"

#include "math/SparseMatrixCSR.h"
#include "math/SparseMatrixCSRGeneral.h"
#include "boost/filesystem.hpp"


#include "mechanics/dofSubMatrixSolvers/SolverMUMPS.h"
#include "mechanics/dofSubMatrixSolvers/SolverEigen.h"
#include "sys/types.h"
#include "sys/sysinfo.h"

#include "mechanics/sections/SectionPlane.h"


constexpr int dim = 2;
constexpr double thickness = 1.0;

// material
constexpr double youngsModulus = 2.1e5;
constexpr double poissonsRatio = 0.3;


// integration
constexpr bool performLineSearch = false;
constexpr bool automaticTimeStepping = false;
constexpr double timeStep = 1e-0;
constexpr double minTimeStep = 1e-1;
constexpr double maxTimeStep = 1e-0;
constexpr double toleranceDisp = 1e-6;
constexpr double simulationTime = 1.0;
constexpr double loadFactor = -10.0;
constexpr double maxIterations = 1;

const Eigen::Vector2d directionX = Eigen::Vector2d::UnitX();
const Eigen::Vector2d directionY = Eigen::Vector2d::UnitY();


int main(int argc, char* argv[])
{

    NuTo::Structure structure(dim);
    structure.SetNumTimeDerivatives(0);
    structure.SetVerboseLevel(10);
    structure.SetShowTime(false);
    structure.GetLogger().OpenFile("output_compare");

    auto importContainer = structure.ImportFromGmsh(argv[1]);

    const int interpolationTypeId = importContainer[0].second;
    structure.InterpolationTypeAdd(interpolationTypeId, eDof::DISPLACEMENTS, eTypeOrder::EQUIDISTANT1);
    structure.ElementTotalConvertToInterpolationType();

    const int materialId = structure.ConstitutiveLawCreate(eConstitutiveType::LINEAR_ELASTIC_ENGINEERING_STRESS);
    structure.ConstitutiveLawSetParameterDouble(materialId, eConstitutiveParameter::YOUNGS_MODULUS, youngsModulus);
    structure.ConstitutiveLawSetParameterDouble(materialId, eConstitutiveParameter::POISSONS_RATIO, poissonsRatio);

    structure.ElementTotalSetConstitutiveLaw(materialId);


    auto section = NuTo::SectionPlane::Create(thickness, true);
    structure.ElementTotalSetSection(section);


    structure.GetLogger() << "*********************************** \n"
                          << "**      constraints              ** \n"
                          << "*********************************** \n\n";

    Eigen::VectorXd coordinates(dim);


    int groupNodesLeftBoundary = structure.GroupCreate(eGroupId::Nodes);
    structure.GroupAddNodeCoordinateRange(groupNodesLeftBoundary, 0, -1.e-6, +1.e-6);
    structure.ConstraintLinearSetDisplacementNodeGroup(groupNodesLeftBoundary, Vector2d::UnitX(), 0.0);
    structure.ConstraintLinearSetDisplacementNodeGroup(groupNodesLeftBoundary, Vector2d::UnitY(), 0.0);

    structure.GetLogger() << "*********************************** \n"
                          << "**      load                     ** \n"
                          << "*********************************** \n\n";


    int loadNodeGroup = structure.GroupCreate(eGroupId::Nodes);

    coordinates[0] = 60;
    coordinates[1] = 0;
    structure.GroupAddNodeRadiusRange(loadNodeGroup, coordinates, 0, 1.e-6);
    int loadId = structure.ConstraintLinearSetDisplacementNodeGroup(loadNodeGroup, Vector2d::UnitY(), 1);


    structure.GetLogger() << "*********************************** \n"
                          << "**      visualization            ** \n"
                          << "*********************************** \n\n";

    int groupAllElements = 9999;
    structure.GroupCreate(groupAllElements, eGroupId::Elements);
    structure.GroupAddElementsTotal(groupAllElements);
    structure.AddVisualizationComponent(groupAllElements, eVisualizeWhat::DISPLACEMENTS);


    structure.GetLogger() << "*********************************** \n"
                          << "**      integration scheme       ** \n"
                          << "*********************************** \n\n";


    NuTo::NewmarkDirect myIntegrationScheme(&structure);
    boost::filesystem::path resultPath(boost::filesystem::path(getenv("HOME")).string() +
                                       std::string("/results/feti/compare/"));

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

    myIntegrationScheme.AddResultGroupNodeForce("forces", loadNodeGroup);

    coordinates[0] = 20;
    coordinates[1] = 10;
    int grpNodes_output_disp = structure.GroupCreate(eGroupId::Nodes);
    structure.GroupAddNodeRadiusRange(grpNodes_output_disp, coordinates, 0, 1.e-6);
    myIntegrationScheme.AddResultNodeDisplacements("displacements",
                                                   structure.GroupGetMemberIds(grpNodes_output_disp)[0]);

    myIntegrationScheme.AddTimeDependentConstraint(loadId, dispRHS);

    myIntegrationScheme.SetSolver(
            std::make_unique<
                    NuTo::SolverEigen<Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>>>>());
    //    myIntegrationScheme.SetSolver(std::make_unique<NuTo::SolverMUMPS>());


    structure.GetLogger() << "*********************************** \n"
                          << "**      solve                    ** \n"
                          << "*********************************** \n\n";

    myIntegrationScheme.Solve(simulationTime);

    structure.GetLogger() << "Total number of Dofs: \t" << structure.GetNumTotalDofs() << "\n\n";
}
