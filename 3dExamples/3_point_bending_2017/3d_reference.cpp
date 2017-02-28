
#include <boost/filesystem.hpp>

#include "mechanics/structures/unstructured/Structure.h"

#include "mechanics/timeIntegration/NewmarkDirect.h"

#include "mechanics/nodes/NodeEnum.h"
#include "mechanics/groups/GroupEnum.h"
#include "mechanics/sections/SectionEnum.h"
#include "mechanics/constitutive/ConstitutiveEnum.h"
#include "visualize/VisualizeEnum.h"
#include "mechanics/interpolationtypes/InterpolationTypeEnum.h"
#include "mechanics/elements/IpDataEnum.h"
#include "mechanics/integrationtypes/IntegrationTypeEnum.h"

#include "../../AuxFunctions.h"


using std::cout;
using std::endl;
using NuTo::Constitutive::eConstitutiveType;
using NuTo::Constitutive::eConstitutiveParameter;
using NuTo::Node::eDof;
using NuTo::Interpolation::eTypeOrder;
using NuTo::Interpolation::eShapeType;
using NuTo::eGroupId;
using NuTo::eVisualizeWhat;
using NuTo::eSectionType;



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char* argv[]) {

    // geometry
    const   int         dimension                   = 3;

    // material
    const   double      youngsModulusMatrix         = 6.5e3;
    const   double      poissonsRatioMatrix         = 0.2;
    const   double      tensileStrength             = std::stod(argv[3]);
    const   double      compressiveStrength         = std::stod(argv[4]);
    const   double      fractureEnergy              = std::stod(argv[5]);
    const   double      nonlocalRadius              = std::stod(argv[6]);

    const   double      youngsModulusFiber          = std::stod(argv[7]);
    const   double      crossSectionFiber           = std::stod(argv[8]);
    const   double      circumferenceFiber          = std::stod(argv[9]);

    const   double interfaceNormalStiffness         = std::stod(argv[10]);
    const   double alpha                            = std::stod(argv[11]);
    const   double maxBondStress                    = std::stod(argv[12]);
    const   double residualBondStress               = std::stod(argv[13]);
    const   double slipAtMaxBondStress              = std::stod(argv[14]);
    const   double slipAtResidualBondStress         = std::stod(argv[15]);


    // integration
    const   double      timeStep                    = 1.e-2;
    const   double      minTimeStep                 = 1.e-4;
    const   double      maxTimeStep                 = 1.e-1;
    const   double      automaticTimeStepping       = true;

    const   double      toleranceDisp               = 1e-6;
    const   double      simulationTime              = 1.0;
    const   double      loadFactor                  = -0.5;
    const   double      maxIterations              = 10;

    const   std::string subDirectory                = argv[16];
    boost::filesystem::path resultPath(boost::filesystem::path(getenv("HOME")).string() + std::string("/results/results_3_point_bending_reference/") + subDirectory);
    boost::filesystem::create_directory(resultPath);

    const auto directionX = Eigen::Matrix<double, dimension, 1>::UnitX();
    const auto directionY = Eigen::Matrix<double, dimension, 1>::UnitY();
    const auto directionZ = Eigen::Matrix<double, dimension, 1>::UnitZ();

    const int x_component = 0;
    const int y_component = 1;
    const int z_component = 2;

    std::map<std::string, std::string> parameters;
    parameters.emplace("dimension",                    std::to_string(dimension));
    parameters.emplace("youngsModulusMatrix",          std::to_string(youngsModulusMatrix));
    parameters.emplace("poissonsRatioMatrix",          std::to_string(poissonsRatioMatrix));
    parameters.emplace("tensileStrength",              std::to_string(tensileStrength));
    parameters.emplace("compressiveStrength",          std::to_string(compressiveStrength));
    parameters.emplace("fractureEnergy",               std::to_string(fractureEnergy));
    parameters.emplace("nonlocalRadius",               std::to_string(nonlocalRadius));
    parameters.emplace("youngsModulusFiber",           std::to_string(youngsModulusFiber));
    parameters.emplace("crossSectionFiber",            std::to_string(crossSectionFiber));
    parameters.emplace("circumferenceFiber",           std::to_string(circumferenceFiber));
    parameters.emplace("youngsModulusFiber",           std::to_string(youngsModulusFiber));
    parameters.emplace("interfaceNormalStiffness",     std::to_string(interfaceNormalStiffness));
    parameters.emplace("alpha",                        std::to_string(alpha));
    parameters.emplace("maxBondStress",                std::to_string(maxBondStress));
    parameters.emplace("residualBondStress",           std::to_string(residualBondStress));
    parameters.emplace("slipAtMaxBondStress",          std::to_string(slipAtMaxBondStress));
    parameters.emplace("slipAtResidualBondStress",     std::to_string(slipAtResidualBondStress));
    parameters.emplace("resultPath",                   resultPath.string());
    NuTo::WriteSimulationParameters(parameters, resultPath.string());

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // structure
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    NuTo::Structure structure(dimension);
    structure.SetShowTime(true);
    structure.SetVerboseLevel(10);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Matrix
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    auto eleGroupAndInterpolationTypeMatrix = structure.ImportFromGmsh(argv[1]);

    const int groupEleDamage = eleGroupAndInterpolationTypeMatrix[0].first;
    const int interpolationTypeIdDamage = eleGroupAndInterpolationTypeMatrix[0].second;
    structure.InterpolationTypeAdd(interpolationTypeIdDamage, eDof::DISPLACEMENTS,    eTypeOrder::EQUIDISTANT2);
    structure.InterpolationTypeAdd(interpolationTypeIdDamage, eDof::NONLOCALEQSTRAIN, eTypeOrder::EQUIDISTANT1);

    const int groupEleLinear = eleGroupAndInterpolationTypeMatrix[1].first;
    const int interpolationTypeIdLinear = eleGroupAndInterpolationTypeMatrix[1].second;
    structure.InterpolationTypeAdd(interpolationTypeIdLinear, eDof::DISPLACEMENTS,    eTypeOrder::EQUIDISTANT2);
//    structure.InterpolationTypeAdd(interpolationTypeIdLinear, eDof::NONLOCALEQSTRAIN, eTypeOrder::EQUIDISTANT1);

    // section
    int sectionId = structure.SectionCreate(eSectionType::VOLUME);
    structure.ElementTotalSetSection(sectionId);

    // material
    const int materialIdMatrixDamage = structure.ConstitutiveLawCreate(eConstitutiveType::GRADIENT_DAMAGE_ENGINEERING_STRESS);
    structure.ConstitutiveLawSetParameterDouble(materialIdMatrixDamage, eConstitutiveParameter::YOUNGS_MODULUS, youngsModulusMatrix);
    structure.ConstitutiveLawSetParameterDouble(materialIdMatrixDamage, eConstitutiveParameter::POISSONS_RATIO, poissonsRatioMatrix);
    structure.ConstitutiveLawSetParameterDouble(materialIdMatrixDamage, eConstitutiveParameter::TENSILE_STRENGTH, tensileStrength);
    structure.ConstitutiveLawSetParameterDouble(materialIdMatrixDamage, eConstitutiveParameter::COMPRESSIVE_STRENGTH, compressiveStrength);
    structure.ConstitutiveLawSetParameterDouble(materialIdMatrixDamage, eConstitutiveParameter::NONLOCAL_RADIUS, nonlocalRadius);
    structure.ConstitutiveLawSetParameterDouble(materialIdMatrixDamage, eConstitutiveParameter::FRACTURE_ENERGY, fractureEnergy);
    structure.ElementGroupSetConstitutiveLaw(groupEleDamage, materialIdMatrixDamage);

    const int materialIdMatrixLinear = structure.ConstitutiveLawCreate(eConstitutiveType::LINEAR_ELASTIC_ENGINEERING_STRESS);
    structure.ConstitutiveLawSetParameterDouble(materialIdMatrixLinear, eConstitutiveParameter::YOUNGS_MODULUS, youngsModulusMatrix);
    structure.ConstitutiveLawSetParameterDouble(materialIdMatrixLinear, eConstitutiveParameter::POISSONS_RATIO, poissonsRatioMatrix);
    structure.ElementGroupSetConstitutiveLaw(groupEleLinear, materialIdMatrixLinear);

    structure.ElementTotalConvertToInterpolationType();

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Fibers
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//    // section
//    const int sectionIdFiber = structure.SectionCreate(eSectionType::TRUSS);
//    structure.SectionSetArea(sectionIdFiber, crossSectionFiber);
//
//    // material
//    const int materialIdFiber = structure.ConstitutiveLawCreate(NuTo::Constitutive::eConstitutiveType::LINEAR_ELASTIC_ENGINEERING_STRESS);
//    structure.ConstitutiveLawSetParameterDouble(materialIdFiber, NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, youngsModulusFiber);
//
//    auto eleGroupAndInterpolationTypeFibers = structure.ImportFromGmsh(argv[2]);
//
//    const int groupEleFibers = eleGroupAndInterpolationTypeFibers[0].first;
//    const int interpolationTypeIdFibers = eleGroupAndInterpolationTypeFibers[0].second;
//    structure.InterpolationTypeAdd(interpolationTypeIdFibers, eDof::DISPLACEMENTS,    eTypeOrder::EQUIDISTANT1);
//
//    structure.ElementGroupSetSection(groupEleFibers,           sectionIdFiber);
//    structure.ElementGroupSetConstitutiveLaw(groupEleFibers,   materialIdFiber);
//    structure.ElementConvertToInterpolationType(groupEleFibers);
//
//    structure.AddVisualizationComponent(groupEleFibers, eVisualizeWhat::DISPLACEMENTS);
//    structure.AddVisualizationComponent(groupEleFibers, eVisualizeWhat::CONSTITUTIVE);
//
//    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//    // constraints
//    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//    int groupConstraintNodes = structure.GroupCreate(eGroupId::Nodes);
//    structure.GroupAddNodesFromElements(groupConstraintNodes, groupEleFibers);
//
//    int numNearestNeighbours = 1;
//
//    auto nodeIds = structure.GroupGetMemberIds(groupConstraintNodes);
//    for (const int iNode : nodeIds)
//        structure.ConstraintLinearEquationNodeToElementCreate(iNode, groupEleDamage, NuTo::Node::eDof::DISPLACEMENTS);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Interface
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//    const int materialIdBond = structure.ConstitutiveLawCreate(NuTo::Constitutive::eConstitutiveType::FIBRE_MATRIX_BOND_STRESS_SLIP);
//    structure.ConstitutiveLawSetParameterDouble(materialIdBond, NuTo::Constitutive::eConstitutiveParameter::NORMAL_STIFFNESS, interfaceNormalStiffness);
//    structure.ConstitutiveLawSetParameterDouble(materialIdBond, NuTo::Constitutive::eConstitutiveParameter::ALPHA, alpha);
//    structure.ConstitutiveLawSetParameterDouble(materialIdBond, NuTo::Constitutive::eConstitutiveParameter::MAX_BOND_STRESS, maxBondStress);
//    structure.ConstitutiveLawSetParameterDouble(materialIdBond, NuTo::Constitutive::eConstitutiveParameter::RESIDUAL_BOND_STRESS, residualBondStress);
//    structure.ConstitutiveLawSetParameterDouble(materialIdBond, NuTo::Constitutive::eConstitutiveParameter::SLIP_AT_MAX_BOND_STRESS, slipAtMaxBondStress);
//    structure.ConstitutiveLawSetParameterDouble(materialIdBond, NuTo::Constitutive::eConstitutiveParameter::SLIP_AT_RESIDUAL_BOND_STRESS, slipAtResidualBondStress);
//
//    const int sectionIdBond = structure.SectionCreate(eSectionType::FIBRE_MATRIX_BOND);
//    structure.SectionSetCircumference(sectionIdBond, circumferenceFiber);
//
//    const int interpolationTypeIdInterface = structure.InterpolationTypeCreate(NuTo::Interpolation::eShapeType::INTERFACE);
//    structure.InterpolationTypeAdd(interpolationTypeIdInterface, eDof::COORDINATES,   eTypeOrder::EQUIDISTANT1);
//    structure.InterpolationTypeAdd(interpolationTypeIdInterface, eDof::DISPLACEMENTS, eTypeOrder::EQUIDISTANT1);
//
//    auto pairGroupFiberGroupBond = structure.InterfaceElementsCreate(groupEleFibers, interpolationTypeIdInterface, interpolationTypeIdFibers);
//
//    const int groupEleFiber = pairGroupFiberGroupBond.first;
//    const int groupEleBond = pairGroupFiberGroupBond.second;
//
//    structure.ElementGroupSetSection(groupEleFiber,   sectionIdFiber);
//    structure.ElementGroupSetSection(groupEleBond,    sectionIdBond);
//
//    structure.ElementGroupSetConstitutiveLaw(groupEleFiber,   materialIdFiber);
//    structure.ElementGroupSetConstitutiveLaw(groupEleBond,    materialIdBond);
//
//    structure.AddVisualizationComponent(groupEleFiber, eVisualizeWhat::DISPLACEMENTS);
//    structure.AddVisualizationComponent(groupEleFiber, eVisualizeWhat::CONSTITUTIVE);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // boundary conditions
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    Eigen::VectorXd nodeCoords(dimension);

    Eigen::VectorXd nodeCoordsCenter(dimension);
    nodeCoordsCenter << 30, 0, 0;
    int groupNodesLeftBoundary = structure.GroupCreate(eGroupId::Nodes);
    structure.GroupAddNodeCylinderRadiusRange(groupNodesLeftBoundary, nodeCoordsCenter, directionZ, 0, 1.e-3);
    structure.ConstraintLinearSetDisplacementNodeGroup(groupNodesLeftBoundary, directionX, 0);
    structure.ConstraintLinearSetDisplacementNodeGroup(groupNodesLeftBoundary, directionY, 0);


    nodeCoordsCenter << 130, 0, 0;
    int groupNodesRightBoundary = structure.GroupCreate(eGroupId::Nodes);
    structure.GroupAddNodeCylinderRadiusRange(groupNodesRightBoundary, nodeCoordsCenter, directionZ, 0, 1.e-3);
    structure.ConstraintLinearSetDisplacementNodeGroup(groupNodesRightBoundary, directionY, 0);

    int groupNodesSymmetryInZ = structure.GroupCreate(eGroupId::Nodes);
    structure.GroupAddNodeCoordinateRange(groupNodesSymmetryInZ, z_component, 20 -1.e-3, 20 +1.e-3);

    structure.ConstraintLinearSetDisplacementNodeGroup(groupNodesSymmetryInZ, directionZ, 0);


    structure.ConstraintInfo(1);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // loads
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    nodeCoordsCenter << 80,40,0;
    int loadNodeGroup = structure.GroupCreate(eGroupId::Nodes);
    structure.GroupAddNodeCylinderRadiusRange(loadNodeGroup, nodeCoordsCenter,directionZ,0,1.e-3);

    int loadId = structure.ConstraintLinearSetDisplacementNodeGroup(loadNodeGroup, directionY, 0);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // time integration
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    NuTo::NewmarkDirect timeIntegration(&structure);


    timeIntegration.SetTimeStep                 ( timeStep                  );
    timeIntegration.SetMinTimeStep              ( minTimeStep                  );
    timeIntegration.SetMaxTimeStep              ( maxTimeStep                  );
    timeIntegration.SetAutomaticTimeStepping    ( automaticTimeStepping                );
    timeIntegration.SetMaxNumIterations         ( maxIterations            );
    timeIntegration.SetResultDirectory          ( resultPath.string(), true );
    timeIntegration.SetToleranceResidual        ( eDof::DISPLACEMENTS, toleranceDisp );

    Eigen::Matrix2d dispRHS;
    dispRHS(0, 0) = 0;
    dispRHS(1, 0) = simulationTime;
    dispRHS(0, 1) = 0;
    dispRHS(1, 1) = loadFactor;



    timeIntegration.AddResultGroupNodeForce("force", loadNodeGroup);

    nodeCoordsCenter << 80,40,10;
    int grpNodes_output_disp = structure.GroupCreate(eGroupId ::Nodes);
    structure.GroupAddNodeRadiusRange(grpNodes_output_disp, nodeCoordsCenter, 0, 7e-1);
    timeIntegration.AddResultNodeDisplacements("displacements", structure.GroupGetMemberIds(grpNodes_output_disp)[0]);

    timeIntegration.AddTimeDependentConstraint(loadId, dispRHS);


    try{
        timeIntegration.Solve(simulationTime);

    }
    catch(NuTo::MechanicsException e)
    {
        std::cout << e.ErrorMessage() << std::endl;
        std::cout << "structure.GetNumTotalDofs(): \t" << structure.GetNumTotalDofs() << std::endl;
        return EXIT_FAILURE;
    }

    // visualization
//    structure.AddVisualizationComponent(groupEleDamage, eVisualizeWhat::DISPLACEMENTS);

    structure.AddVisualizationComponent(groupEleDamage, eVisualizeWhat::DAMAGE);
    structure.SetVisualizationType(groupEleDamage, NuTo::eVisualizationType::POINTS);

//    structure.AddVisualizationComponent(groupEleLinear, eVisualizeWhat::DISPLACEMENTS);
//    structure.AddVisualizationComponent(groupEleLinear, eVisualizeWhat::DAMAGE);


    std::cout << "structure.GetNumTotalDofs(): \t" << structure.GetNumTotalDofs() << std::endl;

}
