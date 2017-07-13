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
    static const bool mAutomaticTimeStepping = false;

    static constexpr double mMatrixYoungsModulus = 4.0e4; // concrete
    static constexpr double mMatrixPoissonsRatio = 0.2;
    static constexpr double mMatrixThickness = 10;
    static constexpr double mMatrixNonlocalRadius = 2;
    static constexpr double mMatrixTensileStrength = 3;
    static constexpr double mMatrixCompressiveStrength = 30;
    static constexpr double mMatrixFractureEnergy = 1;

    static constexpr double mFibreYoungsModulus = 2.1e7; // steel
    static constexpr double mFibrePoissonsRatio = 0.2;
    static constexpr double mFibreCrossSection = 0.1;

    static constexpr double mInterfaceNormalStiffness = 1e8;
    static constexpr double mAlpha = 1;
    static constexpr double mMaxBondStress = 4e4;
    static constexpr double mResidualBondStress = 1e3;
    static constexpr double mSlipAtMaxBondStress = 1;
    static constexpr double mSlipAtResidualBondStress = 5;

    static constexpr double mTimeStep = 1e-2;
    static constexpr double mMinTimeStep = 1e-5;
    static constexpr double mMaxTimeStep = 1e-2;
    static constexpr double mToleranceForce = 1e-6;
    static constexpr double mSimulationTime = 1.0;
    static constexpr double mLoad = 3.0;

    static const boost::filesystem::path mOutputPath;
    static const boost::filesystem::path mMeshFilePathMatrix;
    static const boost::filesystem::path mMeshFilePathFibre;

    static const NuTo::FullVector<double, dimension> mDirectionX;
    static const NuTo::FullVector<double, dimension> mDirectionY;
    static const NuTo::FullVector<double, dimension> mDirectionZ;
};

const boost::filesystem::path Parameters::mOutputPath("/home/phuschke/3d_fibre_pullout_interface/");
const boost::filesystem::path Parameters::mMeshFilePathMatrix(
        "/home/phuschke/develop/nuto/myNutoExamples/MeshFiles/3d/3d_fibre_pullout_matrix.msh");

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
        myIntegrationScheme.SetResultDirectory(Parameters::mOutputPath.string(), true);

        std::cout << "***********************************" << std::endl;
        std::cout << "**      Section                  **" << std::endl;
        std::cout << "***********************************" << std::endl;

        int matrixSection = myStructure.SectionCreate(NuTo::Section::VOLUME);

        int fibreSection = myStructure.SectionCreate(NuTo::Section::TRUSS);
        myStructure.SectionSetArea(fibreSection, Parameters::mFibreCrossSection);

        std::cout << "***********************************" << std::endl;
        std::cout << "**      Material                 **" << std::endl;
        std::cout << "***********************************" << std::endl;

        int matrixMaterial = myStructure.ConstitutiveLawCreate(
                NuTo::Constitutive::eConstitutiveType::LINEAR_ELASTIC_ENGINEERING_STRESS);
        myStructure.ConstitutiveLawSetParameterDouble(matrixMaterial,
                                                      NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS,
                                                      Parameters::mMatrixYoungsModulus);
        myStructure.ConstitutiveLawSetParameterDouble(matrixMaterial,
                                                      NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO,
                                                      Parameters::mMatrixPoissonsRatio);

        int fibreMaterial = myStructure.ConstitutiveLawCreate(
                NuTo::Constitutive::eConstitutiveType::LINEAR_ELASTIC_ENGINEERING_STRESS);
        myStructure.ConstitutiveLawSetParameterDouble(fibreMaterial,
                                                      NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS,
                                                      Parameters::mFibreYoungsModulus);
        myStructure.ConstitutiveLawSetParameterDouble(fibreMaterial,
                                                      NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO,
                                                      Parameters::mFibrePoissonsRatio);

        int interfaceMaterial =
                myStructure.ConstitutiveLawCreate(NuTo::Constitutive::eConstitutiveType::FIBRE_MATRIX_BOND_STRESS_SLIP);
        myStructure.ConstitutiveLawSetParameterDouble(interfaceMaterial,
                                                      NuTo::Constitutive::eConstitutiveParameter::NORMAL_STIFFNESS,
                                                      Parameters::mInterfaceNormalStiffness);
        myStructure.ConstitutiveLawSetParameterDouble(
                interfaceMaterial, NuTo::Constitutive::eConstitutiveParameter::ALPHA, Parameters::mAlpha);
        myStructure.ConstitutiveLawSetParameterDouble(interfaceMaterial,
                                                      NuTo::Constitutive::eConstitutiveParameter::MAX_BOND_STRESS,
                                                      Parameters::mMaxBondStress);
        myStructure.ConstitutiveLawSetParameterDouble(interfaceMaterial,
                                                      NuTo::Constitutive::eConstitutiveParameter::RESIDUAL_BOND_STRESS,
                                                      Parameters::mResidualBondStress);
        myStructure.ConstitutiveLawSetParameterDouble(
                interfaceMaterial, NuTo::Constitutive::eConstitutiveParameter::SLIP_AT_MAX_BOND_STRESS,
                Parameters::mSlipAtMaxBondStress);
        myStructure.ConstitutiveLawSetParameterDouble(
                interfaceMaterial, NuTo::Constitutive::eConstitutiveParameter::SLIP_AT_RESIDUAL_BOND_STRESS,
                Parameters::mSlipAtResidualBondStress);

        std::cout << "***********************************" << std::endl;
        std::cout << "**      Interpolation Type       **" << std::endl;
        std::cout << "***********************************" << std::endl;

        int matrixInterpolationType = myStructure.InterpolationTypeCreate(NuTo::Interpolation::eShapeType::BRICK3D);
        myStructure.InterpolationTypeAdd(matrixInterpolationType, NuTo::Node::COORDINATES,
                                         NuTo::Interpolation::EQUIDISTANT1);
        myStructure.InterpolationTypeAdd(matrixInterpolationType, NuTo::Node::DISPLACEMENTS,
                                         NuTo::Interpolation::EQUIDISTANT1);

        int fibreInterpolationType = myStructure.InterpolationTypeCreate(NuTo::Interpolation::eShapeType::TRUSSXD);
        myStructure.InterpolationTypeAdd(fibreInterpolationType, NuTo::Node::COORDINATES,
                                         NuTo::Interpolation::EQUIDISTANT1);
        myStructure.InterpolationTypeAdd(fibreInterpolationType, NuTo::Node::DISPLACEMENTS,
                                         NuTo::Interpolation::EQUIDISTANT1);

        int interfaceInterpolationType =
                myStructure.InterpolationTypeCreate(NuTo::Interpolation::eShapeType::INTERFACE);
        myStructure.InterpolationTypeAdd(interfaceInterpolationType, NuTo::Node::COORDINATES,
                                         NuTo::Interpolation::EQUIDISTANT1);
        myStructure.InterpolationTypeAdd(interfaceInterpolationType, NuTo::Node::DISPLACEMENTS,
                                         NuTo::Interpolation::EQUIDISTANT1);

        std::cout << "***********************************" << std::endl;
        std::cout << "**      Matrix                   **" << std::endl;
        std::cout << "***********************************" << std::endl;

        NuTo::FullMatrix<int, Eigen::Dynamic, Eigen::Dynamic> createdGroupIdMatrix =
                myStructure.ImportFromGmsh(Parameters::mMeshFilePathMatrix.string(),
                                           NuTo::ElementData::CONSTITUTIVELAWIP, NuTo::IpData::eIpDataType::STATICDATA);
        int groupIdMatrix = createdGroupIdMatrix.GetValue(0, 0);

        myStructure.ElementGroupSetInterpolationType(groupIdMatrix, matrixInterpolationType);
        myStructure.ElementGroupSetSection(groupIdMatrix, matrixSection);
        myStructure.ElementGroupSetConstitutiveLaw(groupIdMatrix, matrixMaterial);

        myStructure.ElementTotalConvertToInterpolationType(1e-6, 10);


        std::cout << "***********************************" << std::endl;
        std::cout << "**      Boundary Conditions      **" << std::endl;
        std::cout << "***********************************" << std::endl;

        int groupNodeBCLeft = myStructure.GroupCreate(NuTo::Groups::eGroupId::Nodes);
        myStructure.GroupAddNodeCoordinateRange(groupNodeBCLeft, 0, -1e-6, +1e-6);
        myStructure.ConstraintLinearSetDisplacementNodeGroup(groupNodeBCLeft, Parameters::mDirectionX, 0);
        myStructure.ConstraintLinearSetDisplacementNodeGroup(groupNodeBCLeft, Parameters::mDirectionY, 0);
        myStructure.ConstraintLinearSetDisplacementNodeGroup(groupNodeBCLeft, Parameters::mDirectionZ, 0);


        std::cout << "***********************************" << std::endl;
        std::cout << "**      Fibre                    **" << std::endl;
        std::cout << "***********************************" << std::endl;

        NuTo::FullVector<double, Eigen::Dynamic> nodeCoordinates(dimension);
        std::set<NuTo::Node::eAttributes> dofs;
        dofs.insert(NuTo::Node::COORDINATES);
        dofs.insert(NuTo::Node::DISPLACEMENTS);

        nodeCoordinates[0] = 5;
        nodeCoordinates[1] = 5;
        nodeCoordinates[2] = 5;
        int node01 = myStructure.NodeCreate(nodeCoordinates, dofs);

        nodeCoordinates[0] = 10;
        nodeCoordinates[1] = 5;
        nodeCoordinates[2] = 5;
        int node02 = myStructure.NodeCreate(nodeCoordinates, dofs);

        NuTo::FullVector<int, Eigen::Dynamic> nodeIndicesFibre(2);
        nodeIndicesFibre[0] = node01;
        nodeIndicesFibre[1] = node02;

        int fibreElementID = myStructure.ElementCreate(fibreInterpolationType, nodeIndicesFibre,
                                                       NuTo::ElementData::eElementDataType::CONSTITUTIVELAWIP,
                                                       NuTo::IpData::eIpDataType::NOIPDATA);
        myStructure.ElementSetSection(fibreElementID, fibreSection);
        myStructure.ElementSetConstitutiveLaw(fibreElementID, fibreMaterial);

        std::cout << "***********************************" << std::endl;
        std::cout << "**      Interface                **" << std::endl;
        std::cout << "***********************************" << std::endl;


        NuTo::FullVector<int, Eigen::Dynamic> nodeIndicesInterface(4);
        nodeIndicesInterface[0] = 26;
        nodeIndicesInterface[1] = 22;
        nodeIndicesInterface[2] = node02;
        nodeIndicesInterface[3] = node01;
        int interfaceElementID = myStructure.ElementCreate(interfaceInterpolationType, nodeIndicesInterface,
                                                           NuTo::ElementData::eElementDataType::CONSTITUTIVELAWIP,
                                                           NuTo::IpData::eIpDataType::STATICDATA);
        myStructure.ElementSetConstitutiveLaw(interfaceElementID, interfaceMaterial);

        // myStructure.InterfaceElementsCreate(groupIdFibre, interfaceInterpolationType, interfaceMaterial,
        // fibreInterpolationType, fibreMaterial, fibreSection);


        std::cout << "***********************************" << std::endl;
        std::cout << "**      Loads                	   **" << std::endl;
        std::cout << "***********************************" << std::endl;

        //        int groupNodeBCRight = myStructure.GroupCreate(NuTo::Groups::eGroupId::Nodes);
        //        myStructure.GroupAddNodeCoordinateRange(groupNodeBCRight, 0, 1-1e-6, 1+1e-6);
        //
        //        int timeDependentConstraint = myStructure.ConstraintLinearSetDisplacementNodeGroup(groupNodeBCRight,
        //        Parameters::mDirectionX, 1);

        int timeDependentConstraint =
                myStructure.ConstraintLinearSetDisplacementNode(node02, Parameters::mDirectionX, 1);

        std::cout << "***********************************" << std::endl;
        std::cout << "**      Visualization            **" << std::endl;
        std::cout << "***********************************" << std::endl;

        myStructure.AddVisualizationComponentDisplacements();
        myStructure.AddVisualizationComponentConstitutive();

        std::cout << "***********************************" << std::endl;
        std::cout << "**      Solver                   **" << std::endl;
        std::cout << "***********************************" << std::endl;

        myStructure.NodeBuildGlobalDofs();
        myStructure.Info();

        myStructure.ElementCheckCoefficientMatrix_0(1e-8);
        myStructure.CheckCoefficientMatrix_0(1e-8, true);

        myStructure.NodeBuildGlobalDofs();

        myStructure.CalculateMaximumIndependentSets();


        NuTo::FullMatrix<double, 2, 2> timeDependentLoad;
        timeDependentLoad(0, 0) = 0;
        timeDependentLoad(1, 0) = 1.0 * Parameters::mSimulationTime;

        timeDependentLoad(0, 1) = 0;
        timeDependentLoad(1, 1) = 1.0 * Parameters::mLoad;

        myIntegrationScheme.SetTimeDependentConstraint(timeDependentConstraint, timeDependentLoad);

        myIntegrationScheme.Solve(Parameters::mSimulationTime);
    }
    catch (NuTo::MechanicsException& e)
    {
        std::cout << e.ErrorMessage();
    }
    catch (NuTo::MathException& e)
    {
        std::cout << e.ErrorMessage();
    }

    std::cout << "***********************************" << std::endl;
    std::cout << "**      END                      **" << std::endl;
    std::cout << "***********************************" << std::endl;
}
