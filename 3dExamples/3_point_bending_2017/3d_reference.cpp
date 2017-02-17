
#include <boost/filesystem.hpp>


#include "mechanics/structures/unstructured/Structure.h"
#include "mechanics/structures/StructureBase.h"

#include "mechanics/timeIntegration/NewmarkDirect.h"

#include "mechanics/nodes/NodeEnum.h"
#include "mechanics/groups/GroupEnum.h"
#include "mechanics/sections/SectionEnum.h"
#include "mechanics/constitutive/ConstitutiveEnum.h"
#include "visualize/VisualizeEnum.h"
#include "mechanics/interpolationtypes/InterpolationTypeEnum.h"
#include "mechanics/elements/IpDataEnum.h"
#include "mechanics/integrationtypes/IntegrationTypeEnum.h"

#include "mechanics/sections/SectionEnum.h"



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

// geometry
constexpr   int         dimension                   = 3;
constexpr   double      thickness                   = 1.0;

// material
constexpr   double      youngsModulus               = 4.0e4;
constexpr   double      poissonsRatio               = 0.2;
constexpr   double      tensileStrength             = 3;
constexpr   double      compressiveStrength         = 200;
constexpr   double      fractureEnergy              = 0.005;
constexpr   double      nonlocalRadius              = 0.5;

// integration
constexpr   double      timeStep                    = 1.e-2;
constexpr   double      minTimeStep                 = 1.e-4;
constexpr   double      maxTimeStep                 = 1.e-1;
constexpr   double      automaticTimeStepping       = true;

constexpr   double      toleranceDisp               = 1e-6;
constexpr   double      simulationTime              = 1.0;
constexpr   double      loadFactor                  = -0.5;
constexpr   double      maxInterations              = 10;

const auto directionX = Eigen::Matrix<double, dimension, 1>::UnitX();
const auto directionY = Eigen::Matrix<double, dimension, 1>::UnitY();
const auto directionZ = Eigen::Matrix<double, dimension, 1>::UnitZ();

constexpr int x_component = 0;
constexpr int y_component = 1;
constexpr int z_component = 2;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char* argv[]) {
    NuTo::Structure structure(dimension);
    structure.SetShowTime(false);

    std::string meshFile = argv[1];


    auto bla = structure.ImportFromGmsh(meshFile);


    const int interpolationTypeId = bla[0].second;
    structure.InterpolationTypeAdd(interpolationTypeId, eDof::DISPLACEMENTS,    eTypeOrder::EQUIDISTANT2);
    structure.InterpolationTypeAdd(interpolationTypeId, eDof::NONLOCALEQSTRAIN, eTypeOrder::EQUIDISTANT1);

    structure.ElementTotalConvertToInterpolationType(1.e-6, 10);

    structure.SetVerboseLevel(10);


    // section
    int sectionId = structure.SectionCreate(eSectionType::VOLUME);
    structure.ElementTotalSetSection(sectionId);

    // material

    int materialId = structure.ConstitutiveLawCreate(eConstitutiveType::GRADIENT_DAMAGE_ENGINEERING_STRESS);
    structure.ConstitutiveLawSetParameterDouble(materialId, eConstitutiveParameter::YOUNGS_MODULUS, youngsModulus);
    structure.ConstitutiveLawSetParameterDouble(materialId, eConstitutiveParameter::POISSONS_RATIO, poissonsRatio);
    structure.ConstitutiveLawSetParameterDouble(materialId, eConstitutiveParameter::TENSILE_STRENGTH, tensileStrength);
    structure.ConstitutiveLawSetParameterDouble(materialId, eConstitutiveParameter::COMPRESSIVE_STRENGTH,
                                                compressiveStrength);
    structure.ConstitutiveLawSetParameterDouble(materialId, eConstitutiveParameter::NONLOCAL_RADIUS, nonlocalRadius);
    structure.ConstitutiveLawSetParameterDouble(materialId, eConstitutiveParameter::FRACTURE_ENERGY, fractureEnergy);




    structure.ElementTotalSetConstitutiveLaw(materialId);

    // visualization
    int groupAllElements = 9999;
    structure.GroupCreate(groupAllElements, eGroupId::Elements);
    structure.GroupAddElementsTotal(groupAllElements);
    structure.AddVisualizationComponent(groupAllElements, eVisualizeWhat::DISPLACEMENTS);
    structure.AddVisualizationComponent(groupAllElements, eVisualizeWhat::DAMAGE);
    structure.AddVisualizationComponent(groupAllElements, eVisualizeWhat::NONLOCAL_EQ_STRAIN);


    // constraints
    Eigen::VectorXd nodeCoords(dimension);


    Eigen::VectorXd nodeCoordsCenter(dimension);
    nodeCoordsCenter << 30, 0, 0;
    int groupNodesLeftBoundary = structure.GroupCreate(eGroupId::Nodes);
    structure.GroupAddNodeCylinderRadiusRange(groupNodesLeftBoundary, nodeCoordsCenter, directionZ, 0, 1.e-3);
    structure.ConstraintLinearSetDisplacementNodeGroup(groupNodesLeftBoundary, directionX, 0);
    structure.ConstraintLinearSetDisplacementNodeGroup(groupNodesLeftBoundary, directionY, 0);
    structure.ConstraintLinearSetDisplacementNodeGroup(groupNodesLeftBoundary, directionZ, 0);

    nodeCoordsCenter << 130, 0, 0;
    int groupNodesRightBoundary = structure.GroupCreate(eGroupId::Nodes);
    structure.GroupAddNodeCylinderRadiusRange(groupNodesRightBoundary, nodeCoordsCenter, directionZ, 0, 1.e-3);
    structure.ConstraintLinearSetDisplacementNodeGroup(groupNodesRightBoundary, directionY, 0);
    structure.ConstraintLinearSetDisplacementNodeGroup(groupNodesRightBoundary, directionZ, 0);
    structure.ConstraintLinearSetDisplacementNodeGroup(groupNodesRightBoundary, directionX, 0);

    structure.ConstraintInfo(1);

    // loads
    nodeCoordsCenter << 80,40,0;
    int loadNodeGroup = structure.GroupCreate(eGroupId::Nodes);
    structure.GroupAddNodeCylinderRadiusRange(loadNodeGroup, nodeCoordsCenter,directionZ,0,1.e-3);

    int loadId = structure.ConstraintLinearSetDisplacementNodeGroup(loadNodeGroup, directionY, 0);

    // time integration
    NuTo::NewmarkDirect timeIntegration(&structure);
    boost::filesystem::path resultPath(boost::filesystem::initial_path().string() + std::string("/results_3_point_bending_reference"));

    timeIntegration.SetTimeStep                 ( timeStep                  );
    timeIntegration.SetMinTimeStep              ( minTimeStep                  );
    timeIntegration.SetMaxTimeStep              ( maxTimeStep                  );
    timeIntegration.SetAutomaticTimeStepping    ( automaticTimeStepping                );
    timeIntegration.SetMaxNumIterations         ( maxInterations            );
    timeIntegration.SetResultDirectory          ( resultPath.string(), true );
    timeIntegration.SetToleranceResidual        ( eDof::DISPLACEMENTS, toleranceDisp );

    Eigen::Matrix2d dispRHS;
    dispRHS(0, 0) = 0;
    dispRHS(1, 0) = simulationTime;
    dispRHS(0, 1) = 0;
    dispRHS(1, 1) = loadFactor;



    timeIntegration.AddResultGroupNodeForce("force", loadNodeGroup);

    nodeCoordsCenter << 80,40,0;
    int grpNodes_output_disp = structure.GroupCreate(eGroupId ::Nodes);
    structure.GroupAddNodeRadiusRange(grpNodes_output_disp, nodeCoordsCenter, 0, 7e-1);
    timeIntegration.AddResultNodeDisplacements("displacements", structure.GroupGetMemberIds(grpNodes_output_disp)[0]);

    timeIntegration.AddTimeDependentConstraint(loadId, dispRHS);

    timeIntegration.Solve(simulationTime);
    try{

    }
    catch(...)
    {
        std::cout << "structure.GetNumTotalDofs(): \t" << structure.GetNumTotalDofs() << std::endl;
    }


}
