
#include <iostream>
#include <string.h>
#include <vector>
#include "nuto/mechanics/structures/unstructured/StructureFETI.h"
#include "nuto/mechanics/structures/unstructured/Structure.h"
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/SparseQR>
#include "ConjugateProjectedGradient.h"
#include <ctime>
#include <cstdlib>
#include <chrono>
#include "ImportMesh.h"
#include "FetiSolver.h"
#include "nuto/mechanics/timeIntegration/NewmarkDirect.h"
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
constexpr   double      youngsModulus               = 4.0;                // N/mm^2
constexpr   double      poissonsRatio               = 0.2;
constexpr   double      nonlocalRadius              = 1;                   // mm
constexpr   double      fractureEnergy              = 1e-6;                   // N/mm
constexpr   double      compressiveStrength         = 30.e-4;                  // N/mm
constexpr   double      tensileStrength             = 3.e-4;                  // N/mm


// integration
constexpr   bool        performLineSearch           = false;
constexpr   bool        automaticTimeStepping       = true;
constexpr   double      timeStep                    = 1e-2;
constexpr   double      minTimeStep                 = 1e-5;
constexpr   double      maxTimeStep                 = 1e-1;
constexpr   double      toleranceDisp               = 1e-6;
constexpr   double      toleranceNlEqStrain         = 1e-6;
constexpr   double      simulationTime              = 1.0;
constexpr   double      loadFactor                  = -0.3;
constexpr   double      maxInterations              = 10;

const NuTo::FullVector<double, dimension> directionX    = NuTo::FullVector<double, dimension>::UnitX();
const NuTo::FullVector<double, dimension> directionY    = NuTo::FullVector<double, dimension>::UnitY();

void AssignSection(NuTo::Structure& rStructure)
{

    int section00 = rStructure.SectionCreate(NuTo::eSectionType::PLANE_STRAIN);
    rStructure.SectionSetThickness(section00, thickness);

    rStructure.ElementTotalSetSection(section00);
}

void AssignMaterial(NuTo::Structure& rStructure)
{

    int material00 = rStructure.ConstitutiveLawCreate(eConstitutiveType::GRADIENT_DAMAGE_ENGINEERING_STRESS);
    rStructure.ConstitutiveLawSetParameterDouble(material00, eConstitutiveParameter::YOUNGS_MODULUS,       youngsModulus);
    rStructure.ConstitutiveLawSetParameterDouble(material00, eConstitutiveParameter::POISSONS_RATIO,       poissonsRatio);
    rStructure.ConstitutiveLawSetParameterDouble(material00, eConstitutiveParameter::TENSILE_STRENGTH,     tensileStrength);
    rStructure.ConstitutiveLawSetParameterDouble(material00, eConstitutiveParameter::COMPRESSIVE_STRENGTH, compressiveStrength);
    rStructure.ConstitutiveLawSetParameterDouble(material00, eConstitutiveParameter::NONLOCAL_RADIUS,      nonlocalRadius);
    rStructure.ConstitutiveLawSetParameterDouble(material00, eConstitutiveParameter::FRACTURE_ENERGY,      fractureEnergy);

    rStructure.ElementTotalSetConstitutiveLaw(material00);

}


int main(int argc, char* argv[])
{

    try{

        std::string meshFile = "feti_beam.msh";
        std::cout << meshFile << std::endl;

        NuTo::Structure structure(dim);
        structure.SetNumTimeDerivatives(0);
        structure.SetVerboseLevel(10);
        structure.SetShowTime(false);


        structure.ImportFromGmsh(meshFile,eElementDataType::CONSTITUTIVELAWIP,eIpDataType::STATICDATA);


        int interpolationType = structure.InterpolationTypeCreate(NuTo::Interpolation::eShapeType::QUAD2D);
        structure.InterpolationTypeAdd(interpolationType, NuTo::Node::eDof::COORDINATES, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
        structure.InterpolationTypeAdd(interpolationType, NuTo::Node::eDof::DISPLACEMENTS, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
        structure.InterpolationTypeAdd(interpolationType, NuTo::Node::eDof::NONLOCALEQSTRAIN, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
        structure.ElementGroupSetInterpolationType(4001, interpolationType);

        structure.ElementTotalConvertToInterpolationType();
        structure.NodeBuildGlobalDofs();
        AssignMaterial(structure);
        AssignSection(structure);

        cout << "**********************************************" << endl;
        cout << "**  constraints                             **" << endl;
        cout << "**********************************************" << endl;
        Eigen::VectorXd nodeCoords(2);

        int groupNodesLeftBoundary = structure.GroupCreate(eGroupId::Nodes);

//            structure.GroupAddNodeCoordinateRange(groupNodesLeftBoundary,0,-1.e-6,+1.e-6);

        nodeCoords[0] = 0;
        nodeCoords[1] = 0;
        structure.GroupAddNodeRadiusRange(groupNodesLeftBoundary, nodeCoords, 0, 1.e-6);

        structure.ConstraintLinearSetDisplacementNodeGroup(groupNodesLeftBoundary, directionX, 0);
        structure.ConstraintLinearSetDisplacementNodeGroup(groupNodesLeftBoundary, directionY, 0);

        int groupNodesRightBoundary = structure.GroupCreate(eGroupId::Nodes);
//            structure.GroupAddNodeCoordinateRange(groupNodesRightBoundary,0,60-1.e-6,60+1.e-6);

        nodeCoords[0] = 60;
        nodeCoords[1] = 0;
        structure.GroupAddNodeRadiusRange(groupNodesRightBoundary, nodeCoords, 0, 1.e-6);

        structure.ConstraintLinearSetDisplacementNodeGroup(groupNodesRightBoundary, directionX, 0);
        structure.ConstraintLinearSetDisplacementNodeGroup(groupNodesRightBoundary, directionY, 0);

        cout << "**********************************************" << endl;
        cout << "**  load                                    **" << endl;
        cout << "**********************************************" << endl;



        int loadNodeGroup = structure.GroupCreate(eGroupId::Nodes);

        nodeCoords[0] = 30;
        nodeCoords[1] = 10;
        structure.GroupAddNodeRadiusRange(loadNodeGroup, nodeCoords, 0, 1.e-6);

        //    nodeCoords[0] = 29;
        //    nodeCoords[1] = 10;
        //    structure.GroupAddNodeRadiusRange(loadNodeGroup, nodeCoords, 0, 1.e-6);


        int loadId = structure.ConstraintLinearSetDisplacementNodeGroup(loadNodeGroup, directionY, 0);
        //    int loadId = structure.LoadCreateNodeGroupForce(0, loadNodeGroup, directionY, 1);


        std::cout << "***********************************" << std::endl;
        std::cout << "**      Visualization rank = " << std::endl;
        std::cout << "***********************************" << std::endl;
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


        NuTo::NewmarkDirect myIntegrationScheme(&structure);
        boost::filesystem::path resultPath(std::string("/home/phuschke/results/feti/compare"));

        myIntegrationScheme.SetTimeStep                 ( timeStep                  );
        myIntegrationScheme.SetMaxNumIterations         ( maxInterations            );
        myIntegrationScheme.SetMinTimeStep              ( minTimeStep               );
        myIntegrationScheme.SetMaxTimeStep              ( maxTimeStep               );
        myIntegrationScheme.SetAutomaticTimeStepping    ( automaticTimeStepping     );
        myIntegrationScheme.SetResultDirectory          ( resultPath.string(), true );
        myIntegrationScheme.SetPerformLineSearch        ( performLineSearch         );
        myIntegrationScheme.SetToleranceResidual        ( eDof::DISPLACEMENTS, toleranceDisp );
        myIntegrationScheme.SetToleranceResidual        ( eDof::NONLOCALEQSTRAIN, toleranceNlEqStrain );

        myIntegrationScheme.AddResultGroupNodeForce("myforce", loadNodeGroup);
        nodeCoords[0] = 30;
        nodeCoords[1] = 10;
        int grpNodes_output_disp = structure.GroupCreate(eGroupId::Nodes);
        structure.GroupAddNodeRadiusRange(grpNodes_output_disp, nodeCoords, 0, 1.e-6);
        myIntegrationScheme.AddResultNodeDisplacements("mydisplacements", structure.GroupGetMemberIds(grpNodes_output_disp).GetValue(0, 0));

        Eigen::Matrix2d dispRHS;        
        dispRHS(0, 0) = 0;
        dispRHS(1, 0) = simulationTime;
        dispRHS(0, 1) = 0;
        dispRHS(1, 1) = loadFactor;

        myIntegrationScheme.AddTimeDependentConstraint(loadId, dispRHS);
        //    myIntegrationScheme.SetTimeDependentLoadCase(loadId, dispRHS);

        cout << "***********************************" << std::endl;
        cout << "**      Solve                    **" << std::endl;
        cout << "***********************************" << std::endl;

        myIntegrationScheme.Solve(simulationTime);

        myIntegrationScheme.Solve(simulationTime);
        std::string command = "paste " +  resultPath.string() + "/myforce.dat " +  resultPath.string() + "/mydisplacements.dat > " +  resultPath.string() + "/forceDisp.dat";
        system(command.c_str());


    }
    catch(...)
    {
        cout << "Something bad happened" << endl;
    }

    std::cout << "END" << std::endl;



}
