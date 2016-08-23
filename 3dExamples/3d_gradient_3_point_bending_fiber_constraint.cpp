#include "nuto/mechanics/structures/unstructured/Structure.h"
#include "nuto/mechanics/timeIntegration/NewmarkDirect.h"

#include <boost/filesystem.hpp>
#include <iostream>
#include <fstream>
#include <string>
constexpr unsigned int dimension = 3;

class Parameters
{
public:

    static const int mDimension = dimension;

    static const bool mPerformLineSearch = true;
    static const bool mAutomaticTimeStepping = true;

    static constexpr double mMatrixYoungsModulus = 4.0e4;   // concrete
    static constexpr double mMatrixPoissonsRatio = 0.2;
    static constexpr double mMatrixNonlocalRadius = 0.1;
    static constexpr double mMatrixTensileStrength = 3;
    static constexpr double mMatrixCompressiveStrength = 30;
    static constexpr double mMatrixFractureEnergy = 0.01;

    static constexpr double mFibreYoungsModulus = 2.1e5;   // steel
    static constexpr double mFibrePoissonsRatio = 0.2;
    static constexpr double mFibreCrossSection = 0.1;
    static constexpr double mFibreCircumference = 0.6;

    static constexpr double mInterfaceNormalStiffness = 1e6;
    static constexpr double mAlpha = 1;
    static constexpr double mMaxBondStress = 3.e0;
    static constexpr double mResidualBondStress = 0.0;
    static constexpr double mSlipAtMaxBondStress = 0.01;
    static constexpr double mSlipAtResidualBondStress = 1;

    static constexpr double mTimeStep = 1e-1;
    static constexpr double mMinTimeStep = 1e-3;
    static constexpr double mMaxTimeStep = 1e-1;
    static constexpr double mToleranceForce = 1e-4;
    static constexpr double mSimulationTime = 1.0;
    static constexpr double mLoad = -0.04;

    static const std::string mOutputPath;
    static const std::string mMeshFilePathMatrix;
    static const std::string mMeshFilePathFibre;

    static const NuTo::FullVector<double, dimension> mDirectionX;
    static const NuTo::FullVector<double, dimension> mDirectionY;
    static const NuTo::FullVector<double, dimension> mDirectionZ;
};

const std::string Parameters::mOutputPath("/home/phuschke/3d_gradient_3_point_bending_fiber_constraint_no_fiber/");
const std::string Parameters::mMeshFilePathMatrix("/home/phuschke/meshFiles/3d/3d_uniaxial_matrix.msh");
const std::string Parameters::mMeshFilePathFibre("/home/phuschke/meshFiles/3d/trusses.msh");



const NuTo::FullVector<double, dimension> Parameters::mDirectionX = NuTo::FullVector<double, dimension>::UnitX();
const NuTo::FullVector<double, dimension> Parameters::mDirectionY = NuTo::FullVector<double, dimension>::UnitY();
const NuTo::FullVector<double, dimension> Parameters::mDirectionZ = NuTo::FullVector<double, dimension>::UnitZ();

int main(int argc, char* argv[])
{

    try
    {
        std::cout << "***********************************" << std::endl;
        std::cout << "**      Structure                **" << std::endl;
        std::cout << "***********************************" << std::endl;

        NuTo::Structure myStructure(Parameters::mDimension);
        myStructure.SetVerboseLevel(10);
        myStructure.SetShowTime(false);
        myStructure.SetNumProcessors(2);
        std::cout << "***********************************" << std::endl;
        std::cout << "**      Integration Scheme       **" << std::endl;
        std::cout << "***********************************" << std::endl;

        NuTo::NewmarkDirect myIntegrationScheme(&myStructure);
        myIntegrationScheme.SetTimeStep(Parameters::mTimeStep);
        myIntegrationScheme.SetMinTimeStep(Parameters::mMinTimeStep);
        myIntegrationScheme.SetMaxTimeStep(Parameters::mMaxTimeStep);
        myIntegrationScheme.SetToleranceForce(Parameters::mToleranceForce);
        myIntegrationScheme.SetAutomaticTimeStepping(Parameters::mAutomaticTimeStepping);
        myIntegrationScheme.SetPerformLineSearch(Parameters::mPerformLineSearch);
        myIntegrationScheme.SetResultDirectory(Parameters::mOutputPath, true);


        std::cout << "***********************************" << std::endl;
        std::cout << "**      Section                  **" << std::endl;
        std::cout << "***********************************" << std::endl;

        int matrixSection = myStructure.SectionCreate(NuTo::Section::VOLUME);

        int fibreSection = myStructure.SectionCreate(NuTo::Section::TRUSS);
        myStructure.SectionSetArea(fibreSection, Parameters::mFibreCrossSection);

        int bondSection = myStructure.SectionCreate(NuTo::Section::FIBRE_MATRIX_BOND);
        myStructure.SectionSetCircumference(bondSection, Parameters::mFibreCircumference);

        std::cout << "***********************************" << std::endl;
        std::cout << "**      Material                 **" << std::endl;
        std::cout << "***********************************" << std::endl;

//        int matrixMaterial = myStructure.ConstitutiveLawCreate(NuTo::Constitutive::eConstitutiveType::LINEAR_ELASTIC_ENGINEERING_STRESS);
//        myStructure.ConstitutiveLawSetParameterDouble(matrixMaterial, NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, Parameters::mMatrixYoungsModulus);
//        myStructure.ConstitutiveLawSetParameterDouble(matrixMaterial, NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO, Parameters::mMatrixPoissonsRatio);

        int matrixMaterial = myStructure.ConstitutiveLawCreate(NuTo::Constitutive::eConstitutiveType::GRADIENT_DAMAGE_ENGINEERING_STRESS);
        myStructure.ConstitutiveLawSetParameterDouble(matrixMaterial, NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, Parameters::mMatrixYoungsModulus);
        myStructure.ConstitutiveLawSetParameterDouble(matrixMaterial, NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO, Parameters::mMatrixPoissonsRatio);
        myStructure.ConstitutiveLawSetParameterDouble(matrixMaterial, NuTo::Constitutive::eConstitutiveParameter::NONLOCAL_RADIUS, Parameters::mMatrixNonlocalRadius);
        myStructure.ConstitutiveLawSetParameterDouble(matrixMaterial, NuTo::Constitutive::eConstitutiveParameter::TENSILE_STRENGTH, Parameters::mMatrixTensileStrength);
        myStructure.ConstitutiveLawSetParameterDouble(matrixMaterial, NuTo::Constitutive::eConstitutiveParameter::COMPRESSIVE_STRENGTH, Parameters::mMatrixCompressiveStrength);
        myStructure.ConstitutiveLawSetParameterDouble(matrixMaterial, NuTo::Constitutive::eConstitutiveParameter::FRACTURE_ENERGY, Parameters::mMatrixFractureEnergy);

        int fibreMaterial = myStructure.ConstitutiveLawCreate(NuTo::Constitutive::eConstitutiveType::LINEAR_ELASTIC_ENGINEERING_STRESS);
        myStructure.ConstitutiveLawSetParameterDouble(fibreMaterial, NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, Parameters::mFibreYoungsModulus);
        myStructure.ConstitutiveLawSetParameterDouble(fibreMaterial, NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO, Parameters::mFibrePoissonsRatio);

        int interfaceMaterial = myStructure.ConstitutiveLawCreate(NuTo::Constitutive::eConstitutiveType::FIBRE_MATRIX_BOND_STRESS_SLIP);
        myStructure.ConstitutiveLawSetParameterDouble(interfaceMaterial, NuTo::Constitutive::eConstitutiveParameter::NORMAL_STIFFNESS, Parameters::mInterfaceNormalStiffness);
        myStructure.ConstitutiveLawSetParameterDouble(interfaceMaterial, NuTo::Constitutive::eConstitutiveParameter::ALPHA, Parameters::mAlpha);
        myStructure.ConstitutiveLawSetParameterDouble(interfaceMaterial, NuTo::Constitutive::eConstitutiveParameter::MAX_BOND_STRESS, Parameters::mMaxBondStress);
        myStructure.ConstitutiveLawSetParameterDouble(interfaceMaterial, NuTo::Constitutive::eConstitutiveParameter::RESIDUAL_BOND_STRESS, Parameters::mResidualBondStress);
        myStructure.ConstitutiveLawSetParameterDouble(interfaceMaterial, NuTo::Constitutive::eConstitutiveParameter::SLIP_AT_MAX_BOND_STRESS, Parameters::mSlipAtMaxBondStress);
        myStructure.ConstitutiveLawSetParameterDouble(interfaceMaterial, NuTo::Constitutive::eConstitutiveParameter::SLIP_AT_RESIDUAL_BOND_STRESS, Parameters::mSlipAtResidualBondStress);

        std::cout << "***********************************" << std::endl;
        std::cout << "**      Interpolation Type       **" << std::endl;
        std::cout << "***********************************" << std::endl;

        int matrixInterpolationType = myStructure.InterpolationTypeCreate(NuTo::Interpolation::eShapeType::TETRAHEDRON3D);
        myStructure.InterpolationTypeAdd(matrixInterpolationType, NuTo::Node::COORDINATES, NuTo::Interpolation::EQUIDISTANT1);
        myStructure.InterpolationTypeAdd(matrixInterpolationType, NuTo::Node::DISPLACEMENTS, NuTo::Interpolation::EQUIDISTANT2);
        myStructure.InterpolationTypeAdd(matrixInterpolationType, NuTo::Node::NONLOCALEQSTRAIN, NuTo::Interpolation::EQUIDISTANT1);

        int fiberInterpolationType = myStructure.InterpolationTypeCreate(NuTo::Interpolation::eShapeType::TRUSSXD);
        myStructure.InterpolationTypeAdd(fiberInterpolationType, NuTo::Node::COORDINATES, NuTo::Interpolation::EQUIDISTANT1);
        myStructure.InterpolationTypeAdd(fiberInterpolationType, NuTo::Node::DISPLACEMENTS, NuTo::Interpolation::EQUIDISTANT2);

        int interfaceInterpolationType = myStructure.InterpolationTypeCreate(NuTo::Interpolation::eShapeType::INTERFACE);
        myStructure.InterpolationTypeAdd(interfaceInterpolationType, NuTo::Node::COORDINATES, NuTo::Interpolation::EQUIDISTANT1);
        myStructure.InterpolationTypeAdd(interfaceInterpolationType, NuTo::Node::DISPLACEMENTS, NuTo::Interpolation::EQUIDISTANT2);

        std::cout << "***********************************" << std::endl;
        std::cout << "**      Matrix                   **" << std::endl;
        std::cout << "***********************************" << std::endl;

        NuTo::FullMatrix<int, Eigen::Dynamic, Eigen::Dynamic> createdGroupIdMatrix = myStructure.ImportFromGmsh(Parameters::mMeshFilePathMatrix, NuTo::ElementData::CONSTITUTIVELAWIP, NuTo::IpData::eIpDataType::STATICDATA);
        int groupEleMatrix = createdGroupIdMatrix.GetValue(0, 0);

        myStructure.ElementGroupSetSection(groupEleMatrix, matrixSection);
        myStructure.ElementGroupSetConstitutiveLaw(groupEleMatrix, matrixMaterial);
        myStructure.ElementTotalConvertToInterpolationType(1.e-6, 10);
        myStructure.ElementGroupSetInterpolationType(groupEleMatrix, matrixInterpolationType);
        myStructure.InterpolationTypeSetIntegrationType(matrixInterpolationType, NuTo::IntegrationType::IntegrationType3D4NGauss4Ip, NuTo::IpData::eIpDataType::STATICDATA);
        myStructure.ElementTotalConvertToInterpolationType(1.e-6, 10);

        std::cout << "***********************************" << std::endl;
        std::cout << "**      Boundary Conditions      **" << std::endl;
        std::cout << "***********************************" << std::endl;

        NuTo::FullVector<double, 3> nodeCoords;
        nodeCoords(0,0) = 0;
        nodeCoords(1,0) = 0;
        nodeCoords(2,0) = 0;

        int groupNodeBCLeft = myStructure.GroupCreate(NuTo::Groups::eGroupId::Nodes);
        myStructure.GroupAddNodeCylinderRadiusRange(groupNodeBCLeft, nodeCoords, Parameters::mDirectionZ, 0, 1.e-6);

        myStructure.ConstraintLinearSetDisplacementNodeGroup(groupNodeBCLeft, Parameters::mDirectionX, 0);
        myStructure.ConstraintLinearSetDisplacementNodeGroup(groupNodeBCLeft, Parameters::mDirectionY, 0);

        int groupSingleNodeLeft = myStructure.GroupCreate(NuTo::Groups::eGroupId::Nodes);
        myStructure.GroupAddNodeRadiusRange(groupSingleNodeLeft, nodeCoords, 0, 1.e-6);

        myStructure.ConstraintLinearSetDisplacementNodeGroup(groupSingleNodeLeft, Parameters::mDirectionZ, 0);


        nodeCoords(0,0) = 50;
        nodeCoords(1,0) = 0;
        nodeCoords(2,0) = 0;

        int groupNodeBCRight = myStructure.GroupCreate(NuTo::Groups::eGroupId::Nodes);
        myStructure.GroupAddNodeCylinderRadiusRange(groupNodeBCRight, nodeCoords, Parameters::mDirectionZ, 0, 1.e-6);

        myStructure.ConstraintLinearSetDisplacementNodeGroup(groupNodeBCRight, Parameters::mDirectionY, 0);

        int groupSingleNodeRight = myStructure.GroupCreate(NuTo::Groups::eGroupId::Nodes);
        myStructure.GroupAddNodeRadiusRange(groupSingleNodeRight, nodeCoords, 0, 1.e-6);

        myStructure.ConstraintLinearSetDisplacementNodeGroup(groupSingleNodeRight, Parameters::mDirectionZ, 0);



        std::cout << "***********************************" << std::endl;
        std::cout << "**      Fibre                    **" << std::endl;
        std::cout << "***********************************" << std::endl;

        NuTo::FullMatrix<int, Eigen::Dynamic, Eigen::Dynamic> createdGroupIdFibre = myStructure.ImportFromGmsh(Parameters::mMeshFilePathFibre, NuTo::ElementData::CONSTITUTIVELAWIP, NuTo::IpData::eIpDataType::STATICDATA);
        int groupIdFiber = createdGroupIdFibre.GetValue(0, 0);

        myStructure.ElementGroupSetInterpolationType(groupIdFiber, fiberInterpolationType);
        myStructure.ElementConvertToInterpolationType(groupIdFiber);
        myStructure.ElementGroupSetConstitutiveLaw(groupIdFiber, fibreMaterial);
        myStructure.ElementGroupSetSection(groupIdFiber,fibreSection);

        std::cout << "***********************************" << std::endl;
        std::cout << "**      Constraints              **" << std::endl;
        std::cout << "***********************************" << std::endl;


        int groupMatrixNodes = myStructure.GroupCreate(NuTo::Groups::Nodes);
        myStructure.GroupAddNodeFromElementGroupCoordinateRange(groupMatrixNodes, groupEleMatrix, 0, 0, 10);

        int groupMatrixElements = myStructure.GroupCreate(NuTo::Groups::Elements);
        myStructure.GroupAddElementsFromNodes(groupMatrixElements, groupMatrixNodes, false);

        int groupConstraintNodes = myStructure.GroupCreate(NuTo::Groups::Nodes);
        myStructure.GroupAddNodeFromElementGroupCoordinateRange(groupConstraintNodes, groupIdFiber, 0, 0, 10);

        myStructure.GroupAddNodesFromElements(groupConstraintNodes, groupIdFiber);

        int numNearestNeighbours = 1;

        auto nodeIds = myStructure.GroupGetMemberIds(groupConstraintNodes);

        for (int iNode = 0; iNode < nodeIds.rows(); ++iNode)
        {
            myStructure.ConstraintLinearEquationNodeToElementCreate(nodeIds.at(iNode, 0), groupEleMatrix, NuTo::Node::eDof::DISPLACEMENTS, numNearestNeighbours);
        }
//
//        std::cout << "***********************************" << std::endl;
//        std::cout << "**      Interface                **" << std::endl;
//        std::cout << "***********************************" << std::endl;
//
//        auto pairGroupFiberGroupBond = myStructure.InterfaceElementsCreate(groupIdFiber, interfaceInterpolationType, fiberInterpolationType);
//
//        int groupEleFiber = pairGroupFiberGroupBond.first;
//        std::cout << "groupEleFiber" << groupEleFiber << std::endl;
//        int groupEleBond = pairGroupFiberGroupBond.second;
//        std::cout << "groupEleBond" << groupEleBond << std::endl;

        std::cout << "***********************************" << std::endl;
        std::cout << "**      Loads                    **" << std::endl;
        std::cout << "***********************************" << std::endl;


        nodeCoords(0,0) = 25;
        nodeCoords(1,0) = 10;
        nodeCoords(2,0) = 0;

        int groupLoad = myStructure.GroupCreate(NuTo::Groups::eGroupId::Nodes);
        myStructure.GroupAddNodeCylinderRadiusRange(groupLoad, nodeCoords, Parameters::mDirectionZ, 0, 1.e-6);

        int timeDependentConstraint = myStructure.ConstraintLinearSetDisplacementNodeGroup(groupLoad, Parameters::mDirectionY, 1);


        std::cout << "***********************************" << std::endl;
        std::cout << "**      Visualization            **" << std::endl;
        std::cout << "***********************************" << std::endl;

        myStructure.AddVisualizationComponent(groupEleMatrix, NuTo::VisualizeBase::DISPLACEMENTS);
        myStructure.AddVisualizationComponent(groupEleMatrix, NuTo::VisualizeBase::DAMAGE);
        myStructure.AddVisualizationComponent(groupEleMatrix, NuTo::VisualizeBase::ENGINEERING_STRAIN);
        myStructure.AddVisualizationComponent(groupEleMatrix, NuTo::VisualizeBase::ENGINEERING_STRESS);

//        myStructure.AddVisualizationComponent(groupEleFiber, NuTo::VisualizeBase::DISPLACEMENTS);

//        myStructure.AddVisualizationComponent(groupNewInterfaceElements, NuTo::VisualizeBase::BOND_STRESS);

        std::cout << "***********************************" << std::endl;
        std::cout << "**      Solver                   **" << std::endl;
        std::cout << "***********************************" << std::endl;

        myStructure.NodeBuildGlobalDofs();

//        myStructure.ElementCheckCoefficientMatrix_0(1e-8);
//        myStructure.CheckCoefficientMatrix_0(1e-8, true);
//        myStructure.Info();

        myStructure.CalculateMaximumIndependentSets();


        myIntegrationScheme.AddResultGroupNodeForce("myforce", groupLoad);
        auto nodeID = myStructure.GroupGetMemberIds(groupLoad);
        myIntegrationScheme.AddResultNodeDisplacements("mydisplacements", nodeID.at(0,0));


        NuTo::FullMatrix<double, 2, 2> timeDependentLoad;
        timeDependentLoad(0, 0) = 0;
        timeDependentLoad(1, 0) = 1.0 * Parameters::mSimulationTime;

        timeDependentLoad(0, 1) = 0;
        timeDependentLoad(1, 1) = 1.0 * Parameters::mLoad;

        myIntegrationScheme.AddTimeDependentConstraint(timeDependentConstraint, timeDependentLoad);

//        myStructure.ElementCheckHessian0(768,1e-6,-1,true);
        myIntegrationScheme.Solve(Parameters::mSimulationTime);

        std::string command = "paste " + Parameters::mOutputPath + "myforce.dat " + Parameters::mOutputPath + "mydisplacements.dat > " + Parameters::mOutputPath + "forceDisp.dat";
              system(command.c_str());

        std::cout << "***********************************" << std::endl;
        std::cout << "**      END                      **" << std::endl;
        std::cout << "***********************************" << std::endl;

    } catch (NuTo::MechanicsException& e)
    {
        std::cout << e.ErrorMessage();

    } catch (NuTo::MathException& e)
    {
        std::cout << e.ErrorMessage();

    }
    catch (...)
    {
        std::cout << "Damn it!" << std::endl;

    }

}

