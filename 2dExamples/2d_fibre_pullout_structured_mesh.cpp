#include "nuto/math/FullMatrix.h"
#include "nuto/mechanics/structures/unstructured/Structure.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStrain3D.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStress3D.h"
#include "nuto/mechanics/constitutive/mechanics/GradientDamageEngineeringStress.h"
#include "nuto/mechanics/MechanicsException.h"
#include "nuto/mechanics/constitutive/ConstitutiveStaticDataBase.h"
#include "nuto/mechanics/constitutive/mechanics/ConstitutiveStaticDataGradientDamage1D.h"

#include "nuto/math/SparseMatrixCSRVector2General.h"

#include "nuto/mechanics/timeIntegration/NewmarkDirect.h"
#include "nuto/base/Logger.h"

#include <boost-1_55/boost/filesystem.hpp>
#include <boost-1_55/boost/lexical_cast.hpp>
#include <boost-1_55/boost/progress.hpp>

#include <iostream>
#include <fstream>
#include <string>
#include <array>
#include <cmath>

#include "ANN/ANN.h"
void PointInTet4N()
{
    // this function only works fpr linear tetrahedra
    const int numNodes = 4;
    const int dimension = 3;

    Eigen::Vector3d point;
    point << 0,0.0,0;

    Eigen::MatrixXd nodes(dimension,numNodes);
    nodes << 0,1,0,0,
             0,0,2,0,
             0,0,0,3;

    Eigen::MatrixXd matrices(numNodes,numNodes);
    matrices.col(3) = Eigen::Vector4d::Ones();
    matrices.block<4,3>(0,0) = nodes.transpose();
    const double det = matrices.determinant();

    bool pointInside = true;

    // check if all determinantes have the same sign
    for (int i = 0; i < numNodes and pointInside; ++i)
    {
        auto mat(matrices);
        mat.block<1, dimension>(i, 0) = point.transpose();

        if (mat.determinant() != 0.0)
            pointInside = (std::signbit(mat.determinant()) == std::signbit(det));

    }



}




constexpr unsigned int dimension = 2;
class Parameters
{
public:

    static const int mDimension = dimension;

    static const bool mPerformLineSearch = true;
    static const bool mAutomaticTimeStepping = true;

    static constexpr double mMatrixYoungsModulus = 3.0e4;   // concrete
    static constexpr double mMatrixPoissonsRatio = 0.2;
    static constexpr double mMatrixThickness = 0.2;
    static constexpr double mMatrixNonlocalRadius = 2;
    static constexpr double mMatrixTensileStrength = 3;
    static constexpr double mMatrixCompressiveStrength = 50;
    static constexpr double mMatrixFractureEnergy = 0.01;

    static constexpr double mFibreYoungsModulus = 2.1e7;   // steel
    static constexpr double mFibrePoissonsRatio = 0.2;
    static constexpr double mFibreCrossSection = 0.1;
    static constexpr double mFibreCircumference = 0.6;

    static constexpr double mInterfaceNormalStiffness = 1e6;
    static constexpr double mAlpha = 1;
    static constexpr double mMaxBondStress = 1e6;
    static constexpr double mResidualBondStress = 1e6;
    static constexpr double mSlipAtMaxBondStress = 0.1;
    static constexpr double mSlipAtResidualBondStress = 1;

    static constexpr double mTimeStep = 1e-1;
    static constexpr double mMinTimeStep = 1e-7;
    static constexpr double mMaxTimeStep = 1e-1;
    static constexpr double mToleranceForce = 1e-6;
    static constexpr double mSimulationTime = 1.0;
    static constexpr double mLoad =1;

    static const boost::filesystem::path mOutputPath;
    static const boost::filesystem::path mMeshFilePath;

    static const NuTo::FullVector<double, dimension> mDirectionX;
    static const NuTo::FullVector<double, dimension> mDirectionY;
};

const boost::filesystem::path Parameters::mOutputPath("/home/phuschke/2d_fibre_pullout_strucutred_mesh/");
const boost::filesystem::path Parameters::mMeshFilePath("/home/phuschke/develop/nuto/myNutoExamples/MeshFiles/interface/2d_fiber_structured_mesh.msh");

const NuTo::FullVector<double, dimension> Parameters::mDirectionX = NuTo::FullVector<double, dimension>::UnitX();
const NuTo::FullVector<double, dimension> Parameters::mDirectionY = NuTo::FullVector<double, dimension>::UnitY();

int main(int argc, char* argv[])
{




    std::cout << "hello" << std::endl;
    PointInTet4N();
    std::cout << "hello" << std::endl;
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

        int matrixSection = myStructure.SectionCreate(NuTo::Section::PLANE_STRESS);
        myStructure.SectionSetThickness(matrixSection, Parameters::mMatrixThickness);

        int fibreSection = myStructure.SectionCreate(NuTo::Section::TRUSS);
        myStructure.SectionSetArea(fibreSection, Parameters::mFibreCrossSection);

        int fibreMatrixBond = myStructure.SectionCreate(NuTo::Section::FIBRE_MATRIX_BOND);
        myStructure.SectionSetCircumference(fibreMatrixBond, Parameters::mFibreCircumference);

        std::cout << "***********************************" << std::endl;
        std::cout << "**      Material                 **" << std::endl;
        std::cout << "***********************************" << std::endl;

        int matrixMaterial = myStructure.ConstitutiveLawCreate(NuTo::Constitutive::eConstitutiveType::LINEAR_ELASTIC_ENGINEERING_STRESS);
        myStructure.ConstitutiveLawSetParameterDouble(matrixMaterial, NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, Parameters::mMatrixYoungsModulus);
        myStructure.ConstitutiveLawSetParameterDouble(matrixMaterial, NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO, Parameters::mMatrixPoissonsRatio);

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

        NuTo::FullMatrix<int, Eigen::Dynamic, Eigen::Dynamic> gmshImportedGroupIds = myStructure.ImportFromGmsh(Parameters::mMeshFilePath.string(), NuTo::ElementData::CONSTITUTIVELAWIP, NuTo::IpData::eIpDataType::STATICDATA);

        int groupIdMatrix = gmshImportedGroupIds.GetValue(0, 0);

        int matrixInterpolationType = myStructure.InterpolationTypeCreate(NuTo::Interpolation::eShapeType::TRIANGLE2D);
        myStructure.InterpolationTypeAdd(matrixInterpolationType, NuTo::Node::COORDINATES, NuTo::Interpolation::EQUIDISTANT1);
        myStructure.InterpolationTypeAdd(matrixInterpolationType, NuTo::Node::DISPLACEMENTS, NuTo::Interpolation::EQUIDISTANT2);


        myStructure.InterpolationTypeInfo(matrixInterpolationType);
        myStructure.IntegrationTypeInfo(10);

        int fibreInterpolationType = myStructure.InterpolationTypeCreate(NuTo::Interpolation::eShapeType::TRUSSXD);
        myStructure.InterpolationTypeAdd(fibreInterpolationType, NuTo::Node::COORDINATES, NuTo::Interpolation::EQUIDISTANT1);
        myStructure.InterpolationTypeAdd(fibreInterpolationType, NuTo::Node::DISPLACEMENTS, NuTo::Interpolation::EQUIDISTANT1);

        int interfaceInterpolationType = myStructure.InterpolationTypeCreate(NuTo::Interpolation::eShapeType::INTERFACE);
        myStructure.InterpolationTypeAdd(interfaceInterpolationType, NuTo::Node::COORDINATES, NuTo::Interpolation::EQUIDISTANT1);
        myStructure.InterpolationTypeAdd(interfaceInterpolationType, NuTo::Node::DISPLACEMENTS, NuTo::Interpolation::EQUIDISTANT1);
        myStructure.InterpolationTypeSetIntegrationType(interfaceInterpolationType, NuTo::IntegrationType::IntegrationType1D2NGauss2Ip, NuTo::IpData::eIpDataType::STATICDATA);

        std::cout << "***********************************" << std::endl;
        std::cout << "**      Matrix                   **" << std::endl;
        std::cout << "***********************************" << std::endl;




        myStructure.ElementGroupSetSection(groupIdMatrix, matrixSection);
        myStructure.ElementGroupSetConstitutiveLaw(groupIdMatrix, matrixMaterial);
        myStructure.ElementTotalConvertToInterpolationType(1.e-6, 10);
        myStructure.ElementGroupSetInterpolationType(groupIdMatrix, matrixInterpolationType);
        myStructure.InterpolationTypeSetIntegrationType(matrixInterpolationType, NuTo::IntegrationType::IntegrationType2D3NGauss3Ip, NuTo::IpData::eIpDataType::NOIPDATA);
        myStructure.ElementTotalConvertToInterpolationType(1.e-6, 10);


        std::cout << "***********************************" << std::endl;
        std::cout << "**      Create Extra Nodes       **" << std::endl;
        std::cout << "***********************************" << std::endl;

        NuTo::FullVector<double, dimension> nodeCoords;


        std::set<NuTo::Node::eAttributes> dofs;
        dofs.insert(NuTo::Node::COORDINATES);
        dofs.insert(NuTo::Node::DISPLACEMENTS);

        nodeCoords[0] = 10.0;
        nodeCoords[1] = -10.0;
        int node00 = myStructure.NodeCreate(nodeCoords, dofs);

        nodeCoords[0] = 10.0;
        nodeCoords[1] = 1.0;
        int node01 = myStructure.NodeCreate(nodeCoords, dofs);

        nodeCoords[0] = 70.0;
        nodeCoords[1] = -10.0;
        int node02 = myStructure.NodeCreate(nodeCoords, dofs);

        nodeCoords[0] = 70.0;
        nodeCoords[1] = 1.0;
        int node03 = myStructure.NodeCreate(nodeCoords, dofs);

        std::cout << "***********************************" << std::endl;
        std::cout << "**      Fibre                    **" << std::endl;
        std::cout << "***********************************" << std::endl;

        NuTo::FullVector<int, 2> fiberNodeIds;
        fiberNodeIds[0] = node00;
        fiberNodeIds[1] = node02;
        int fibreElement = myStructure.ElementCreate(fibreInterpolationType, fiberNodeIds, NuTo::ElementData::eElementDataType::CONSTITUTIVELAWIP, NuTo::IpData::eIpDataType::NOIPDATA);

        myStructure.ElementSetSection(fibreElement, fibreSection);
        myStructure.ElementSetConstitutiveLaw(fibreElement, fibreMaterial);


        std::cout << "***********************************" << std::endl;
        std::cout << "**      Interface                **" << std::endl;
        std::cout << "***********************************" << std::endl;

        NuTo::FullVector<int, 4> bondNodeIds;
        bondNodeIds[0] = node01;
        bondNodeIds[1] = node03;
        bondNodeIds[2] = node02;
        bondNodeIds[3] = node00;
        int bondElement = myStructure.ElementCreate(interfaceInterpolationType, bondNodeIds, NuTo::ElementData::eElementDataType::CONSTITUTIVELAWIP, NuTo::IpData::eIpDataType::STATICDATA);

        myStructure.ElementSetSection(bondElement, fibreMatrixBond);
        myStructure.ElementSetConstitutiveLaw(bondElement, interfaceMaterial);



        std::cout << "***********************************" << std::endl;
        std::cout << "**      Boundary Conditions      **" << std::endl;
        std::cout << "***********************************" << std::endl;

        int groupNodeBCLeft = myStructure.GroupCreate(NuTo::Groups::eGroupId::Nodes);
        myStructure.GroupAddNodeCoordinateRange(groupNodeBCLeft, 0, -1e-6, +1e-6);

        myStructure.ConstraintLinearSetDisplacementNodeGroup(groupNodeBCLeft, Parameters::mDirectionX, 0);


        nodeCoords[0] = 0.0;
        nodeCoords[1] = 0.0;
        int nodeLeft = myStructure.NodeGetIdAtCoordinate(nodeCoords, 1e-6);

        myStructure.ConstraintLinearSetDisplacementNode(nodeLeft, Parameters::mDirectionY, 0);

        myStructure.ConstraintLinearEquationNodeToElementCreate(node01,groupIdMatrix,NuTo::Node::eAttributes::DISPLACEMENTS, 0, 5);
        myStructure.ConstraintLinearEquationNodeToElementCreate(node01,groupIdMatrix,NuTo::Node::eAttributes::DISPLACEMENTS, 1, 5);

        myStructure.ConstraintLinearEquationNodeToElementCreate(node03,groupIdMatrix,NuTo::Node::eAttributes::DISPLACEMENTS, 0, 5);
        myStructure.ConstraintLinearEquationNodeToElementCreate(node03,groupIdMatrix,NuTo::Node::eAttributes::DISPLACEMENTS, 1, 5);


//        myStructure.ConstraintLinearEquationCreate(444,node01, NuTo::Node::eAttributes::DISPLACEMENTS,0,1.0,0.0);
//        myStructure.ConstraintLinearEquationAddTerm(444, 20, NuTo::Node::eAttributes::DISPLACEMENTS, 0, -1.0);
//
//        myStructure.ConstraintLinearEquationCreate(555,node01, NuTo::Node::eAttributes::DISPLACEMENTS,1,1.0,0.0);
//        myStructure.ConstraintLinearEquationAddTerm(555, 20, NuTo::Node::eAttributes::DISPLACEMENTS, 1, -1.0);
//
//
//        myStructure.ConstraintLinearEquationCreate(666,node03, NuTo::Node::eAttributes::DISPLACEMENTS,0,1.0,0.0);
//        myStructure.ConstraintLinearEquationAddTerm(666, 21, NuTo::Node::eAttributes::DISPLACEMENTS, 0, -1.0);
//
//        myStructure.ConstraintLinearEquationCreate(777,node03, NuTo::Node::eAttributes::DISPLACEMENTS,1,1.0,0.0);
//        myStructure.ConstraintLinearEquationAddTerm(777, 21, NuTo::Node::eAttributes::DISPLACEMENTS, 1, -1.0);


        std::cout << "***********************************" << std::endl;
        std::cout << "**      Loads                    **" << std::endl;
        std::cout << "***********************************" << std::endl;

        int groupNodeBCRight = myStructure.GroupCreate(NuTo::Groups::eGroupId::Nodes);
        myStructure.GroupAddNodeCoordinateRange(groupNodeBCRight, 0, 160 - 1e-6, 160 + 1e-6);

        int timeDependentConstraint = myStructure.ConstraintLinearSetDisplacementNodeGroup(groupNodeBCRight, Parameters::mDirectionX, 1);

        myStructure.Info();



        std::cout << "***********************************" << std::endl;
        std::cout << "**      Visualization            **" << std::endl;
        std::cout << "***********************************" << std::endl;




        myStructure.AddVisualizationComponentDisplacements();
        myStructure.AddVisualizationComponentConstitutive();
        std::cout << "***********************************" << std::endl;
        std::cout << "**      Solver                   **" << std::endl;
        std::cout << "***********************************" << std::endl;


        myStructure.NodeBuildGlobalDofs();


        myStructure.CalculateMaximumIndependentSets();






        NuTo::FullMatrix<double, 2, 2> timeDependentLoad;
        timeDependentLoad(0, 0) = 0;
        timeDependentLoad(1, 0) = Parameters::mSimulationTime;
        timeDependentLoad(0, 1) = 0;
        timeDependentLoad(1, 1) = Parameters::mLoad;

//        myIntegrationScheme.SetTimeDependentLoadCase(loadCase, timeDependentLoad);
        myIntegrationScheme.SetTimeDependentConstraint(timeDependentConstraint, timeDependentLoad);

        myIntegrationScheme.Solve(Parameters::mSimulationTime);


        std::cout << "\n\n\n Results written to " + Parameters::mOutputPath.string() << std::endl;







    } catch (NuTo::MechanicsException& e)
    {
        std::cout << e.ErrorMessage();

    }
    catch (NuTo::MathException& e)
    {
        std::cout << e.ErrorMessage();

    }
    catch (...)
    {
        std::cout << "Something else went wrong" << std::endl;

    }

}

