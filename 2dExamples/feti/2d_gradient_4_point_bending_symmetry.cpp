#include "nuto/mechanics/structures/unstructured/Structure.h"
#include "nuto/mechanics/timeIntegration/NewmarkDirect.h"
#include "nuto/mechanics/nodes/NodeBase.h"
#include "../../EnumsAndTypedefs.h"
#include "nuto/mechanics/timeIntegration/NewmarkFeti.h"
#include <mpi.h>
#include <boost/filesystem.hpp>
#include <iostream>
#include <fstream>
#include <string>
#include <chrono>


using std::cout;
using std::endl;




constexpr int dimension = 2;

constexpr bool performLineSearch = true;
constexpr bool automaticTimeStepping = true;

constexpr double youngsModulus          = 4.0e4;   // concrete
constexpr double poissonsRatio          = 0.2;
constexpr double thickness              = 50;
constexpr double nonlocalRadius         = 4;
constexpr double tensileStrength        = 3;
constexpr double compressiveStrength    = 30;
constexpr double fractureEnergy         = 0.01;

constexpr double timeStep               = 1e-2;
constexpr double minTimeStep            = 1e-5;
constexpr double maxTimeStep            = 1e-1;
constexpr double toleranceDisp          = 1e-6;
constexpr int    maxInterations         = 10;
constexpr double toleranceNlEqStrain    = 1e-6;
constexpr double simulationTime         = 1.0;
constexpr double load                   = -0.4;


const NuTo::FullVector<double, dimension> directionX = NuTo::FullVector<double, dimension>::UnitX();
const NuTo::FullVector<double, dimension> directionY = NuTo::FullVector<double, dimension>::UnitY();


int main(int argc, char* argv[])
{

    MPI_Init(&argc, &argv);
    const int rank = MPI::COMM_WORLD.Get_rank();

    cout << "**********************************************" << endl;
    cout << "**  strucutre                               **" << endl;
    cout << "**********************************************" << endl;

    std::string meshFile = "fourPointBending.msh_" + std::to_string(rank);
    std::cout << meshFile << std::endl;

    NuTo::StructureFETI structure(dimension);
    structure.SetNumTimeDerivatives(0);

    const int interpolationTypeId = structure.InterpolationTypeCreate(eShapeType::QUAD2D);
    structure.InterpolationTypeAdd(interpolationTypeId, eDof::COORDINATES,     eTypeOrder::EQUIDISTANT1);
    structure.InterpolationTypeAdd(interpolationTypeId, eDof::DISPLACEMENTS,   eTypeOrder::EQUIDISTANT1);
    structure.InterpolationTypeAdd(interpolationTypeId, eDof::NONLOCALEQSTRAIN, eTypeOrder::EQUIDISTANT1);


    structure.ImportMeshJson(meshFile,interpolationTypeId);


    structure.SetVerboseLevel(10);
    structure.SetShowTime(false);


    cout << "**********************************************" << endl;
    cout << "**  section                                 **" << endl;
    cout << "**********************************************" << endl;

    int mySection = structure.SectionCreate(eSectionType::PLANE_STRESS);
    structure.SectionSetThickness(mySection, thickness);
    structure.ElementTotalSetSection(mySection);

    cout << "**********************************************" << endl;
    cout << "**  material                                **" << endl;
    cout << "**********************************************" << endl;


    int matrixMaterial = structure.ConstitutiveLawCreate(eConstitutiveType::GRADIENT_DAMAGE_ENGINEERING_STRESS);
    structure.ConstitutiveLawSetParameterDouble(matrixMaterial, eConstitutiveParameter::YOUNGS_MODULUS, youngsModulus);
    structure.ConstitutiveLawSetParameterDouble(matrixMaterial, eConstitutiveParameter::POISSONS_RATIO, poissonsRatio);
    structure.ConstitutiveLawSetParameterDouble(matrixMaterial, eConstitutiveParameter::TENSILE_STRENGTH, tensileStrength);
    structure.ConstitutiveLawSetParameterDouble(matrixMaterial, eConstitutiveParameter::COMPRESSIVE_STRENGTH, compressiveStrength);
    structure.ConstitutiveLawSetParameterDouble(matrixMaterial, eConstitutiveParameter::NONLOCAL_RADIUS, nonlocalRadius);
    structure.ConstitutiveLawSetParameterDouble(matrixMaterial, eConstitutiveParameter::FRACTURE_ENERGY, fractureEnergy);
    structure.ElementTotalSetConstitutiveLaw(matrixMaterial);


    cout << "**********************************************" << endl;
    cout << "**  bc                                      **" << endl;
    cout << "**********************************************" << endl;
    NuTo::FullVector<double, 2> center;
    constexpr double tol = 1.0e-6;
    // bottom left boundary
    center[0] = 25;
    center[1] = 0;
    int grpNodes_bottom_left = structure.GroupCreate(eGroupId::Nodes);
    structure.GroupAddNodeRadiusRange(grpNodes_bottom_left, center, 0, tol);
    structure.ConstraintLinearSetDisplacementNodeGroup(grpNodes_bottom_left, directionY, 0);

    int group_nodes_tmp_01 = structure.GroupCreate(NuTo::eGroupId::Nodes);
    structure.GroupAddNodeCoordinateRange(group_nodes_tmp_01, 0, 250 - tol, 250 + tol);

    int group_nodes_tmp_02 = structure.GroupCreate(NuTo::eGroupId::Nodes);
    structure.GroupAddNodeCoordinateRange(group_nodes_tmp_02, 1, 10 - tol, 100 + tol);




    structure.ConstraintLinearSetDisplacementNodeGroup(group_nodes_tmp_01, directionX, 0);

    cout << "**********************************************" << endl;
    cout << "**  load                                    **" << endl;
    cout << "**********************************************" << endl;

    // middle top load
    center[0] = 175;
    center[1] = 100;
    int grpNodes_load = structure.GroupCreate(eGroupId::Nodes);

    structure.GroupAddNodeRadiusRange(grpNodes_load, center, 0, 1e-6);


    int load = structure.ConstraintLinearSetDisplacementNodeGroup(grpNodes_load, directionY, 0);

    cout << "**********************************************" << endl;
    cout << "**  visualization                           **" << endl;
    cout << "**********************************************" << endl;

    int groupAllElements = 9999;
    structure.GroupCreate(groupAllElements, eGroupId::Elements);
    structure.GroupAddElementsTotal(groupAllElements);
    structure.AddVisualizationComponent(groupAllElements, eVisualizeWhat::DISPLACEMENTS);
    structure.AddVisualizationComponent(groupAllElements, eVisualizeWhat::DAMAGE);
    structure.AddVisualizationComponent(groupAllElements, eVisualizeWhat::NONLOCAL_EQ_STRAIN);
    structure.AddVisualizationComponent(groupAllElements, eVisualizeWhat::ENGINEERING_STRAIN);


    cout << "**********************************************" << endl;
    cout << "**  integration sheme                       **" << endl;
    cout << "**********************************************" << endl;


    NuTo::NewmarkFeti myIntegrationScheme(&structure);
    boost::filesystem::path resultPath(std::string("/home/phuschke/results/feti/fourPointBending/" + std::to_string(structure.mRank)));

    myIntegrationScheme.SetTimeStep                 ( timeStep                  );
    myIntegrationScheme.SetMaxNumIterations         ( maxInterations            );
    myIntegrationScheme.SetMinTimeStep              ( minTimeStep               );
    myIntegrationScheme.SetMaxTimeStep              ( maxTimeStep               );
    myIntegrationScheme.SetAutomaticTimeStepping    ( automaticTimeStepping     );
    myIntegrationScheme.SetResultDirectory          ( resultPath.string(), true );
    myIntegrationScheme.SetPerformLineSearch        ( performLineSearch         );
    myIntegrationScheme.SetToleranceResidual        ( eDof::DISPLACEMENTS, toleranceDisp );
    myIntegrationScheme.SetToleranceResidual        ( eDof::NONLOCALEQSTRAIN, toleranceNlEqStrain );


    cout << "**********************************************" << endl;
    cout << "**  solver                                  **" << endl;
    cout << "**********************************************" << endl;


//    center[0] = 175;
//    center[1] = 100;
//    int grpNodes_output = structure.GroupCreate(eGroupId::Nodes);
//    structure.GroupAddNodeRadiusRange(grpNodes_output, center, 0, 5e-1);

//    myIntegrationScheme.AddResultGroupNodeForce("myforce", grpNodes_output);

//    center[0] = 240;
//    center[1] = 0;
//    int grpNodes_output_disp = structure.GroupCreate(eGroupId::Nodes);
//    structure.GroupAddNodeRadiusRange(grpNodes_output_disp, center, 0, 7e-1);
//    myIntegrationScheme.AddResultNodeDisplacements("mydisplacements", structure.GroupGetMemberIds(grpNodes_output_disp).GetValue(0, 0));


    NuTo::FullMatrix<double, 2, 2> dispRHS;
    dispRHS(0, 0) = 0;
    dispRHS(1, 0) = simulationTime;
    dispRHS(0, 1) = 0;
    dispRHS(1, 1) = load;

    myIntegrationScheme.AddTimeDependentConstraint(load, dispRHS);

    structure.Info();

    myIntegrationScheme.Solve(simulationTime);

    cout << "**********************************************" << endl;
    cout << "**  end                                     **" << endl;
    cout << "**********************************************" << endl;

    MPI_Finalize();
}

