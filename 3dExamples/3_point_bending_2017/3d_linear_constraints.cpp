
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
constexpr   double      compressiveStrength         = 30;
constexpr   double      fractureEnergy              = 0.1;

// integration
constexpr   double      timeStep                    = 1.0;
constexpr   double      toleranceDisp               = 1e-6;
constexpr   double      simulationTime              = 1.0;
constexpr   double      loadFactor                  = -10;
constexpr   double      maxInterations              = 10;

const auto directionX = Eigen::Matrix<double, dimension, 1>::UnitX();
const auto directionY = Eigen::Matrix<double, dimension, 1>::UnitY();
const auto directionZ = Eigen::Matrix<double, dimension, 1>::UnitZ();

constexpr int x_component = 0;
constexpr int y_component = 1;
constexpr int z_component = 2;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char* argv[])
{
    NuTo::Structure structure(dimension);
    structure.SetShowTime(false);

    std::string meshFile = argv[1];


    structure.ImportFromGmsh(meshFile);



    const int interpolationTypeId = structure.InterpolationTypeCreate(eShapeType::TETRAHEDRON3D);
    structure.InterpolationTypeAdd(interpolationTypeId, eDof::COORDINATES,     eTypeOrder::EQUIDISTANT1);
    structure.InterpolationTypeAdd(interpolationTypeId, eDof::DISPLACEMENTS,   eTypeOrder::EQUIDISTANT1);
    structure.ElementTotalSetInterpolationType(interpolationTypeId);

    structure.ElementTotalSetInterpolationType(interpolationTypeId);
    structure.ElementTotalConvertToInterpolationType(1.e-6, 10);

    structure.SetVerboseLevel(10);
    structure.Info();
    structure.ElementTotalConvertToInterpolationType();

    // section
    int sectionId = structure.SectionCreate(eSectionType::PLANE_STRESS);
    structure.SectionSetThickness(sectionId, thickness);
    structure.ElementTotalSetSection(sectionId);

    // material
    int materialId = structure.ConstitutiveLawCreate(eConstitutiveType::LINEAR_ELASTIC_ENGINEERING_STRESS);
    structure.ConstitutiveLawSetParameterDouble(materialId, eConstitutiveParameter::YOUNGS_MODULUS, youngsModulus);
    structure.ConstitutiveLawSetParameterDouble(materialId, eConstitutiveParameter::POISSONS_RATIO, poissonsRatio);
    structure.ElementTotalSetConstitutiveLaw(materialId);

    // visualization
    int groupAllElements = 9999;
    structure.GroupCreate(groupAllElements, eGroupId::Elements);
    structure.GroupAddElementsTotal(groupAllElements);
    structure.AddVisualizationComponent(groupAllElements, eVisualizeWhat::DISPLACEMENTS);


    // constraints
    Eigen::VectorXd nodeCoords(dimension);



    Eigen::VectorXd nodeCoordsCenter(dimension);
    nodeCoordsCenter << 30,0,0;
    int groupNodesLeftBoundary = structure.GroupCreate(eGroupId::Nodes);
    structure.GroupAddNodeCylinderRadiusRange(groupNodesLeftBoundary, nodeCoordsCenter,directionZ,0,1.e-3);
    structure.ConstraintLinearSetDisplacementNodeGroup(groupNodesLeftBoundary,directionX,0);
    structure.ConstraintLinearSetDisplacementNodeGroup(groupNodesLeftBoundary,directionY,0);
    structure.ConstraintLinearSetDisplacementNodeGroup(groupNodesLeftBoundary,directionZ,0);

//    nodeCoordsCenter << 130,0,0;
//    int groupNodesRightBoundary = structure.GroupCreate(eGroupId::Nodes);
//    structure.GroupAddNodeCylinderRadiusRange(groupNodesRightBoundary, nodeCoordsCenter,directionZ,0,1.e-3);
//    structure.ConstraintLinearSetDisplacementNodeGroup(groupNodesRightBoundary,directionY,0);
//    structure.ConstraintLinearSetDisplacementNodeGroup(groupNodesRightBoundary,directionZ,0);


    int groupNodesLeftInterface = structure.GroupCreate(eGroupId::Nodes);
    structure.GroupAddNodeCoordinateRange(groupNodesLeftInterface, 0, 80-1.e-3,80+1.e-3);



    int groupNodesRightInterface = structure.GroupCreate(eGroupId::Nodes);
    structure.GroupAddNodeCoordinateRange(groupNodesRightInterface, 0, 120-1.e-3,120+1.e-3);

    int constraintId = 1337;
    for (const auto& nodeId : structure.GroupGetMemberIds(groupNodesRightInterface))
    {

        structure.NodeGetCoordinates(nodeId, nodeCoords);
        nodeCoords[0] += 40;


        int groupTmp = structure.GroupCreate(eGroupId::Nodes);
        structure.GroupAddNodeRadiusRange(groupTmp, nodeCoords, 0, 1.e-6);
        int nodeIdTmp = structure.GroupGetMemberIds(groupTmp)[0];

        structure.ConstraintLinearEquationCreate(constraintId, nodeId,eDof::DISPLACEMENTS,x_component,1,0);
        structure.ConstraintLinearEquationAddTerm(constraintId, nodeIdTmp, eDof::DISPLACEMENTS,x_component,-1.0);
        ++constraintId;
        structure.ConstraintLinearEquationCreate(constraintId, nodeId,eDof::DISPLACEMENTS,y_component,1,0);
        structure.ConstraintLinearEquationAddTerm(constraintId, nodeIdTmp, eDof::DISPLACEMENTS,y_component,-1.0);
        ++constraintId;
        structure.ConstraintLinearEquationCreate(constraintId, nodeId,eDof::DISPLACEMENTS,z_component,1,0);
        structure.ConstraintLinearEquationAddTerm(constraintId, nodeIdTmp, eDof::DISPLACEMENTS,z_component,-1.0);
        ++constraintId;
    }


//    nodeCoords << 80,0,0;
//    int groupNodeLeft = structure.GroupCreate(eGroupId::Nodes);
//    structure.GroupAddNodeRadiusRange(groupNodeLeft, nodeCoords, 0, 1.e-6);
//
//
//
//    nodeCoords << 120,0,0;
//    int groupNodeRight = structure.GroupCreate(eGroupId::Nodes);
//    structure.GroupAddNodeRadiusRange(groupNodeRight, nodeCoords, 0, 1.e-6);
//
//
//    int nodeIdLeft  = structure.GroupGetMemberIds(groupNodeLeft)[0];
//    int nodeIdRight = structure.GroupGetMemberIds(groupNodeRight)[0];
//    int constraintId = 1337000;
//    structure.ConstraintLinearEquationCreate( constraintId,  nodeIdLeft,  eDof::DISPLACEMENTS, x_component,  1,0);
//    structure.ConstraintLinearEquationAddTerm(constraintId,  nodeIdRight, eDof::DISPLACEMENTS, x_component,  -1.0);
//    constraintId++;
//    structure.ConstraintLinearEquationCreate( constraintId,  nodeIdLeft,  eDof::DISPLACEMENTS, y_component,  1,0);
//    structure.ConstraintLinearEquationAddTerm(constraintId,  nodeIdRight, eDof::DISPLACEMENTS, y_component,  -1.0);
//    constraintId++;
//    structure.ConstraintLinearEquationCreate( constraintId,  nodeIdLeft,  eDof::DISPLACEMENTS, z_component,  1,0);
//    structure.ConstraintLinearEquationAddTerm(constraintId,  nodeIdRight, eDof::DISPLACEMENTS, z_component,  -1.0);
//



    // loads
    nodeCoordsCenter << 170,40,0;
    int loadNodeGroup = structure.GroupCreate(eGroupId::Nodes);
    structure.GroupAddNodeCylinderRadiusRange(loadNodeGroup, nodeCoordsCenter,directionZ,0,1.e-3);

    std::cout << "loadNodeGroup: \t" << loadNodeGroup << std::endl;
    structure.GroupInfo(19);
    int loadId = structure.ConstraintLinearSetDisplacementNodeGroup(loadNodeGroup, directionY, 0);

    // time integration
    NuTo::NewmarkDirect timeIntegration(&structure);
    boost::filesystem::path resultPath(boost::filesystem::initial_path().string() + std::string("/results_3_point_bending"));

    timeIntegration.SetTimeStep                 ( timeStep                  );
    timeIntegration.SetMaxNumIterations         ( maxInterations            );
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
