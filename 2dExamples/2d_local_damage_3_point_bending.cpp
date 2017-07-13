#include "nuto/mechanics/structures/unstructured/Structure.h"
#include "nuto/mechanics/timeIntegration/NewmarkDirect.h"


#include "../myNutoExamples/EnumsAndTypedefs.h"


#include "nuto/mechanics/nodes/NodeBase.h"

#include "nuto/math/SparseMatrixCSR.h"
#include "nuto/math/SparseMatrixCSRGeneral.h"
#include "nuto/mechanics/dofSubMatrixStorage/BlockSparseMatrix.h"
#include "nuto/mechanics/structures/StructureOutputBlockMatrix.h"


#include <boost/filesystem.hpp>
#include <fstream>
#include <string>
#include <chrono>


// geometry
constexpr int dimension = 2;

// material
constexpr double youngsModulus = 4.0e4;
constexpr double poissonsRatio = 0.2;
constexpr double tensileStrength = 3;
constexpr double compressiveStrength = 30;
constexpr double fractureEnergy = 0.01;

// section
constexpr double thickness = 1.0;

// integration
constexpr bool performLineSearch = false;
constexpr bool automaticTimeStepping = true;
constexpr double timeStep = 1e-2;
constexpr double minTimeStep = 1e-5;
constexpr double maxTimeStep = 1e-1;
constexpr double toleranceForce = 1e-6;
constexpr double simulationTime = 1.0;
constexpr double loadFactor = -0.001;

void WriteParameters(const std::string& subdirectory, boost::filesystem::path resultPath)
{
    std::ofstream file;
    file.open(std::string(resultPath.string() + "parameters.txt"));
    cout << "Writing simulation parameters to file: " << std::string(resultPath.string() + "parameters.txt")
         << std::endl;
    file << "subdirectory         " << subdirectory << std::endl;
    file << "dimension            " << dimension << std::endl;
    file << "performLineSearch    " << performLineSearch << std::endl;
    file << "automaticTimeStepping" << automaticTimeStepping << std::endl;
    file << "youngsModulus        " << youngsModulus << std::endl;
    file << "poissonsRatio        " << poissonsRatio << std::endl;
    file << "thickness            " << thickness << std::endl;
    file << "timeStep             " << timeStep << std::endl;
    file << "minTimeStep          " << minTimeStep << std::endl;
    file << "maxTimeStep          " << maxTimeStep << std::endl;
    file << "toleranceForce       " << toleranceForce << std::endl;
    file << "simulationTime       " << simulationTime << std::endl;
    file << "loadFactor           " << loadFactor << std::endl;

    file.close();
}

int main(int argc, char* argv[])
{

    if (argc != 2)
    {
        std::cout << "input arguments: subdirectory" << std::endl;
        return EXIT_FAILURE;
    }

    std::string subdirectory = argv[1];
    std::cout << subdirectory << endl;

    boost::filesystem::path resultPath(
            std::string("/home/phuschke/results/2d/2d_local_damage_model/" + subdirectory + "/"));
    cout << resultPath.string() << endl;

    const boost::filesystem::path meshFilePath("2d_miehe_symmetric_three_point_bending.msh");

    const NuTo::FullVector<double, dimension> directionX = NuTo::FullVector<double, dimension>::UnitX();
    const NuTo::FullVector<double, dimension> directionY = NuTo::FullVector<double, dimension>::UnitY();


    cout << "**********************************************" << endl;
    cout << "**  strucutre                               **" << endl;
    cout << "**********************************************" << endl;

    NuTo::Structure myStructure(dimension);
    myStructure.SetNumTimeDerivatives(0);
    myStructure.SetShowTime(false);

    cout << "Writing results to:" << endl;
    cout << resultPath.string() << endl;
    boost::filesystem::create_directory(resultPath);


    cout << "**********************************************" << endl;
    cout << "**  integration sheme                       **" << endl;
    cout << "**********************************************" << endl;

    NuTo::NewmarkDirect myIntegrationScheme(&myStructure);

    myIntegrationScheme.SetTimeStep(timeStep);
    myIntegrationScheme.SetMinTimeStep(minTimeStep);
    myIntegrationScheme.SetMaxTimeStep(maxTimeStep);
    myIntegrationScheme.SetToleranceForce(toleranceForce);
    myIntegrationScheme.SetAutomaticTimeStepping(automaticTimeStepping);
    myIntegrationScheme.SetPerformLineSearch(performLineSearch);
    myIntegrationScheme.SetResultDirectory(resultPath.string(), true);

    WriteParameters(subdirectory, resultPath);

    cout << "**********************************************" << endl;
    cout << "**  section                                 **" << endl;
    cout << "**********************************************" << endl;

    int mySection = myStructure.SectionCreate(NuTo::eSectionType::PLANE_STRAIN);
    myStructure.SectionSetThickness(mySection, thickness);

    cout << "**********************************************" << endl;
    cout << "**  material                                **" << endl;
    cout << "**********************************************" << endl;

    int myMaterial = myStructure.ConstitutiveLawCreate(eConstitutiveType::LOCAL_DAMAGE_MODEL);

    myStructure.ConstitutiveLawSetParameterDouble(myMaterial, eConstitutiveParameter::YOUNGS_MODULUS, youngsModulus);
    myStructure.ConstitutiveLawSetParameterDouble(myMaterial, eConstitutiveParameter::POISSONS_RATIO, poissonsRatio);
    myStructure.ConstitutiveLawSetParameterDouble(myMaterial, eConstitutiveParameter::TENSILE_STRENGTH,
                                                  tensileStrength);
    myStructure.ConstitutiveLawSetParameterDouble(myMaterial, eConstitutiveParameter::COMPRESSIVE_STRENGTH,
                                                  compressiveStrength);
    myStructure.ConstitutiveLawSetParameterDouble(myMaterial, eConstitutiveParameter::FRACTURE_ENERGY, fractureEnergy);


    cout << "**********************************************" << endl;
    cout << "**  geometry                                **" << endl;
    cout << "**********************************************" << endl;

    NuTo::FullMatrix<int, Eigen::Dynamic, Eigen::Dynamic> createdGroupIds = myStructure.ImportFromGmsh(
            meshFilePath.string(), eElementDataType::CONSTITUTIVELAWIP, eIpDataType::STATICDATA);
    int groupId = createdGroupIds.GetValue(0, 0);

    int myInterpolationType = myStructure.InterpolationTypeCreate(eShapeType::TRIANGLE2D);
    myStructure.InterpolationTypeAdd(myInterpolationType, eDof::COORDINATES, eTypeOrder::EQUIDISTANT1);
    myStructure.InterpolationTypeAdd(myInterpolationType, eDof::DISPLACEMENTS, eTypeOrder::EQUIDISTANT2);
    myStructure.ElementGroupSetInterpolationType(groupId, myInterpolationType);
    myStructure.InterpolationTypeSetIntegrationType(
            myInterpolationType, NuTo::eIntegrationType::IntegrationType2D3NGauss3Ip, eIpDataType::STATICDATA);

    myStructure.InterpolationTypeInfo(myInterpolationType);

    myStructure.ElementTotalConvertToInterpolationType(1.e-6, 10);

    myStructure.ElementTotalSetSection(mySection);
    myStructure.ElementTotalSetConstitutiveLaw(myMaterial);

    cout << "**********************************************" << endl;
    cout << "**  bc                                      **" << endl;
    cout << "**********************************************" << endl;
    NuTo::FullVector<double, 2> center;
    constexpr double tol = 1.0e-6;
    // bottom left boundary
    center[0] = 0.0;
    center[1] = 0.0;
    int grpNodes_bottom_left = myStructure.GroupCreate(NuTo::eGroupId::Nodes);
    myStructure.GroupAddNodeRadiusRange(grpNodes_bottom_left, center, 0, tol);
    myStructure.ConstraintLinearSetDisplacementNodeGroup(grpNodes_bottom_left, directionY, 0);
    myStructure.ConstraintLinearSetDisplacementNodeGroup(grpNodes_bottom_left, directionX, 0);

    // bottom right boundary
    center[0] = 8.0;
    center[1] = 0.0;

    int grpNodes_bottom_right = myStructure.GroupCreate(NuTo::eGroupId::Nodes);
    myStructure.GroupAddNodeRadiusRange(grpNodes_bottom_right, center, 0, tol);
    myStructure.ConstraintLinearSetDisplacementNodeGroup(grpNodes_bottom_right, directionY, 0);


    cout << "**********************************************" << endl;
    cout << "**  load                                    **" << endl;
    cout << "**********************************************" << endl;

    // middle top load
    center[0] = 4.0;
    center[1] = 2.0;
    int grpNodes_load = myStructure.GroupCreate(NuTo::eGroupId::Nodes);

    myStructure.GroupAddNodeRadiusRange(grpNodes_load, center, 0, 1e-6);


    int loadId = myStructure.ConstraintLinearSetDisplacementNodeGroup(grpNodes_load, directionY, 0);

    cout << "**********************************************" << endl;
    cout << "**  visualization                           **" << endl;
    cout << "**********************************************" << endl;

    myStructure.AddVisualizationComponent(groupId, NuTo::eVisualizeWhat::DISPLACEMENTS);
    myStructure.AddVisualizationComponent(groupId, NuTo::eVisualizeWhat::DAMAGE);

    cout << "**********************************************" << endl;
    cout << "**  solver                                  **" << endl;
    cout << "**********************************************" << endl;


    myStructure.NodeBuildGlobalDofs();
    myStructure.CalculateMaximumIndependentSets();


    //    center[0] = 175;
    //    center[1] = 100;
    //    int grpNodes_output = myStructure.GroupCreate(NuTo::Groups::Nodes);
    //    myStructure.GroupAddNodeRadiusRange(grpNodes_output, center, 0, 5e-1);
    //
    //    myIntegrationScheme.AddResultGroupNodeForce("myforce", grpNodes_output);
    //
    //    center[0] = 240;
    //    center[1] = 0;
    //    int grpNodes_output_disp = myStructure.GroupCreate(NuTo::Groups::Nodes);
    //    myStructure.GroupAddNodeRadiusRange(grpNodes_output_disp, center, 0, 7e-1);
    //    myIntegrationScheme.AddResultNodeDisplacements("mydisplacements",
    //    myStructure.GroupGetMemberIds(grpNodes_output_disp).GetValue(0, 0));

    NuTo::FullMatrix<double, 2, 2> dispRHS;
    dispRHS(0, 0) = 0;
    dispRHS(1, 0) = simulationTime;
    dispRHS(0, 1) = 0;
    dispRHS(1, 1) = loadFactor;

    myIntegrationScheme.AddTimeDependentConstraint(loadId, dispRHS);
    myStructure.ExportVtkDataFileElements(resultPath.string() + "bla.vtu", true);
    myIntegrationScheme.Solve(simulationTime);

    //    std::string command = "paste " +  OutputPath.string() + "myforce.dat " +  OutputPath.string() +
    //    "mydisplacements.dat > " +  OutputPath.string() + "forceDisp.dat";
    //            system(command.c_str());
    cout << "**********************************************" << endl;
    cout << "**  end                                     **" << endl;
    cout << "**********************************************" << endl;
}
