

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
constexpr   double      loadFactor                  = 10;
constexpr   double      maxIterations              = 10;



void AssignSection(NuTo::Structure& structure);
void AssignMaterial(NuTo::Structure &structure);



int main(int argc, char* argv[])
{


    struct sysinfo memInfo;
    sysinfo(&memInfo);
    long long totalVirtualMem = memInfo.totalram;
    totalVirtualMem += memInfo.totalswap;
    totalVirtualMem *= memInfo.mem_unit;

    long long virtualMemUsed = memInfo.totalram - memInfo.freeram;
    virtualMemUsed += memInfo.totalswap - memInfo.freeswap;
    virtualMemUsed *= memInfo.mem_unit;

    long long totalPhysMem = memInfo.totalram;
    totalPhysMem *= memInfo.mem_unit;

    long long physMemUsed = memInfo.totalram - memInfo.freeram;
    physMemUsed *= memInfo.mem_unit;

    std::cout << "totalVirtualMem:  \t" << totalVirtualMem << std::endl;
    std::cout << "virtualMemUsed:   \t" << virtualMemUsed << std::endl;
    std::cout << "totalPhysMem:     \t" << totalPhysMem << std::endl;
    std::cout << "physMemUsed:     \t" << physMemUsed << std::endl;

    NuTo::Structure structure(dim);
    structure.SetNumTimeDerivatives(0);
    structure.SetVerboseLevel(10);
    structure.SetShowTime(false);
    structure.GetLogger().OpenFile("output_compare");

    auto importContainer = structure.ImportFromGmsh(argv[1]);

    const int interpolationTypeId = importContainer[0].second;
    structure.InterpolationTypeAdd(interpolationTypeId, eDof::DISPLACEMENTS,   eTypeOrder::EQUIDISTANT1);
    structure.ElementTotalConvertToInterpolationType();

    AssignMaterial(structure);
    AssignSection(structure);

    structure.GetLogger() << "**********************************************" << "\n";
    structure.GetLogger() << "**  real constraints                        **" << "\n";
    structure.GetLogger() << "**********************************************" << "\n\n";

    Eigen::VectorXd nodeCoords(2);


    int groupNodesLeftBoundary = structure.GroupCreate(eGroupId::Nodes);
    structure.GroupAddNodeCoordinateRange(groupNodesLeftBoundary,0,-1.e-6,+1.e-6);
    structure.ConstraintLinearSetDisplacementNodeGroup(groupNodesLeftBoundary, Vector2d::UnitX(), 0.0);
    structure.ConstraintLinearSetDisplacementNodeGroup(groupNodesLeftBoundary, Vector2d::UnitY(), 0.0);

    structure.GetLogger() << "**********************************************" << "\n";
    structure.GetLogger() << "**  load                                    **" << "\n";
    structure.GetLogger() << "**********************************************" << "\n\n";




    int loadNodeGroup = structure.GroupCreate(eGroupId::Nodes);

    nodeCoords[0] = 60;
    nodeCoords[1] = 0;
    structure.GroupAddNodeRadiusRange(loadNodeGroup, nodeCoords, 0, 1.e-6);
    int loadId = structure.ConstraintLinearSetDisplacementNodeGroup(loadNodeGroup, Vector2d::UnitY(), 1);


    structure.GetLogger() << "***********************************" << "\n";
    structure.GetLogger() << "**      Visualization            **" << "\n";
    structure.GetLogger() << "***********************************" << "\n\n";

    int groupAllElements = 9999;
    structure.GroupCreate(groupAllElements, eGroupId::Elements);
    structure.GroupAddElementsTotal(groupAllElements);
    structure.AddVisualizationComponent(groupAllElements, eVisualizeWhat::DISPLACEMENTS);



    structure.GetLogger() << "**********************************************" << "\n";
    structure.GetLogger() << "**  integration sheme                       **" << "\n";
    structure.GetLogger() << "**********************************************" << "\n\n";




    NuTo::NewmarkDirect myIntegrationScheme(&structure);
    boost::filesystem::path resultPath(std::string("/home/phuschke/results/feti/compare"));

    myIntegrationScheme.SetTimeStep                 ( timeStep                  );
    myIntegrationScheme.SetMaxNumIterations         ( maxIterations            );
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

    myIntegrationScheme.AddResultGroupNodeForce("myforcebla", loadNodeGroup);

    nodeCoords[0] = 20;
    nodeCoords[1] = 10;
    int grpNodes_output_disp = structure.GroupCreate(eGroupId::Nodes);
    structure.GroupAddNodeRadiusRange(grpNodes_output_disp, nodeCoords, 0, 1.e-6);
    myIntegrationScheme.AddResultNodeDisplacements("mydisplacements", structure.GroupGetMemberIds(grpNodes_output_disp)[0]);

    myIntegrationScheme.AddTimeDependentConstraint(loadId, dispRHS);
//    myIntegrationScheme.SetTimeDependentLoadCase(loadId, dispRHS);


    myIntegrationScheme.SetSolver(std::make_unique<NuTo::SolverEigen<Eigen::SparseLU<Eigen::SparseMatrix<double>,Eigen::COLAMDOrdering<int>>>>());
//    myIntegrationScheme.SetSolver(std::make_unique<NuTo::SolverMUMPS>());


    structure.GetLogger() << "***********************************" << "\n";
    structure.GetLogger() << "**      Solve                    **" << "\n";
    structure.GetLogger() << "***********************************" << "\n\n";

    myIntegrationScheme.Solve(simulationTime);

    structure.GetLogger()   << "Total number of Dofs: \t"
                            << structure.GetNumTotalDofs() << "\n\n";

}


void AssignSection(NuTo::Structure& structure)
{
    auto section = NuTo::SectionPlane::Create(thickness,true);
    structure.ElementTotalSetSection(section);}

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
