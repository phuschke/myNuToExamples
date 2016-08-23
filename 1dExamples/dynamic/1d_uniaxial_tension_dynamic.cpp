#include "nuto/mechanics/structures/unstructured/Structure.h"
#include "nuto/mechanics/timeIntegration/NewmarkDirect.h"

#include <boost/filesystem.hpp>
#include <iostream>
#include <fstream>
#include <string>

class Parameters
{
public:

    static constexpr int mDimension = 1;

    static constexpr bool mPerformLineSearch = true;
    static constexpr bool mAutomaticTimeStepping = true;

    static constexpr double mMatrixYoungsModulus = 4.0e4;   // concrete
    static constexpr double mMatrixPoissonsRatio = 0.2;
    static constexpr double mMatrixCrossSection = 1;
    static constexpr double mMatrixNonlocalRadius = 2;
    static constexpr double mMatrixTensileStrength = 3;
    static constexpr double mMatrixCompressiveStrength = 30;
    static constexpr double mMatrixFractureEnergy = 1;
    static constexpr double mMatrixDensity= 1e+2;   // concrete

    static constexpr double mTimeStep = 1e-1;
    static constexpr double mMinTimeStep = 1e-2;
    static constexpr double mMaxTimeStep = 1e-1;
    static constexpr double mToleranceForce = 1e-8;
    static constexpr double mSimulationTime = 1.0;
    static constexpr double mLoad = 5;

    static const boost::filesystem::path mOutputPath;
    static const boost::filesystem::path mMeshFilePathMatrix;
    static const boost::filesystem::path mMeshFilePathFibre;

    static const NuTo::FullVector<double, Parameters::mDimension> mDirectionX;

};

const boost::filesystem::path Parameters::mOutputPath("/home/phuschke/1d_uniaxial_tension_dynamic/");
const NuTo::FullVector<double, Parameters::mDimension> Parameters::mDirectionX = NuTo::FullVector<double, Parameters::mDimension>::UnitX();

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
        myStructure.SetNumTimeDerivatives(2);
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
        //myIntegrationScheme.SetCheckCoefficientMatrix(true);

        std::cout << "***********************************" << std::endl;
        std::cout << "**      Section                  **" << std::endl;
        std::cout << "***********************************" << std::endl;

        int matrixSection = myStructure.SectionCreate(NuTo::Section::TRUSS);
        myStructure.SectionSetArea(matrixSection, Parameters::mMatrixCrossSection);


        std::cout << "***********************************" << std::endl;
        std::cout << "**      Material                 **" << std::endl;
        std::cout << "***********************************" << std::endl;

        int matrixMaterial = myStructure.ConstitutiveLawCreate(NuTo::Constitutive::eConstitutiveType::LINEAR_ELASTIC_ENGINEERING_STRESS);
        myStructure.ConstitutiveLawSetParameterDouble(matrixMaterial, NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, Parameters::mMatrixYoungsModulus);
        myStructure.ConstitutiveLawSetParameterDouble(matrixMaterial, NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO, Parameters::mMatrixPoissonsRatio);
        myStructure.ConstitutiveLawSetParameterDouble(matrixMaterial, NuTo::Constitutive::eConstitutiveParameter::DENSITY, Parameters::mMatrixDensity);

//        int matrixMaterial = myStructure.ConstitutiveLawCreate(NuTo::Constitutive::eConstitutiveType::GRADIENT_DAMAGE_ENGINEERING_STRESS);
//        myStructure.ConstitutiveLawSetParameterDouble(matrixMaterial, NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, Parameters::mMatrixYoungsModulus);
//        myStructure.ConstitutiveLawSetParameterDouble(matrixMaterial, NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO, Parameters::mMatrixPoissonsRatio);
//        myStructure.ConstitutiveLawSetParameterDouble(matrixMaterial, NuTo::Constitutive::eConstitutiveParameter::TENSILE_STRENGTH, Parameters::mMatrixTensileStrength);
//        myStructure.ConstitutiveLawSetParameterDouble(matrixMaterial, NuTo::Constitutive::eConstitutiveParameter::COMPRESSIVE_STRENGTH, Parameters::mMatrixCompressiveStrength);
//        myStructure.ConstitutiveLawSetParameterDouble(matrixMaterial, NuTo::Constitutive::eConstitutiveParameter::NONLOCAL_RADIUS, Parameters::mMatrixNonlocalRadius);
//        myStructure.ConstitutiveLawSetParameterDouble(matrixMaterial, NuTo::Constitutive::eConstitutiveParameter::FRACTURE_ENERGY, Parameters::mMatrixFractureEnergy);
//        myStructure.ConstitutiveLawSetParameterDouble(matrixMaterial, NuTo::Constitutive::eConstitutiveParameter::DENSITY, Parameters::mMatrixDensity);

        std::cout << "***********************************" << std::endl;
        std::cout << "**      Interpolation Type       **" << std::endl;
        std::cout << "***********************************" << std::endl;

        int myInterpolationType = myStructure.InterpolationTypeCreate("Truss1D");
        myStructure.InterpolationTypeAdd(myInterpolationType, NuTo::Node::COORDINATES, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
        myStructure.InterpolationTypeAdd(myInterpolationType, NuTo::Node::DISPLACEMENTS, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);

        std::cout << "***********************************" << std::endl;
        std::cout << "**      Matrix                   **" << std::endl;
        std::cout << "***********************************" << std::endl;

        //create nodes
        int numElements = 5;
        int numNodes = numElements + 1;
        double length = 10;
        double deltaX = length / (numElements);

        for (int iNode = 0; iNode < numNodes; iNode++)
        {
            NuTo::FullVector<double, Eigen::Dynamic> coordinates(1);
            coordinates(0) = iNode * deltaX;
            myStructure.NodeCreate(iNode, coordinates);

        }



        //create elements
        for (int iEle = 0; iEle < numElements; iEle++)
        {
            NuTo::FullVector<int, Eigen::Dynamic> nodes(2);
            nodes(0) = iEle;
            nodes(1) = iEle + 1;
            myStructure.ElementCreate(myInterpolationType, nodes);
        }

        myStructure.SetVerboseLevel(10);
        myStructure.ElementTotalConvertToInterpolationType();
        myStructure.ElementTotalSetSection(matrixSection);
        myStructure.ElementTotalSetConstitutiveLaw(matrixMaterial);

        std::cout << "***********************************" << std::endl;
        std::cout << "**      Boundary Conditions      **" << std::endl;
        std::cout << "***********************************" << std::endl;

        int groupNodeBCLeft = myStructure.GroupCreate(NuTo::Groups::eGroupId::Nodes);
        myStructure.GroupAddNodeCoordinateRange(groupNodeBCLeft, 0, 0.0 - 1.0e-6, 0.0 + 1.0e-6);
        myStructure.ConstraintLinearSetDisplacementNodeGroup(groupNodeBCLeft, Parameters::mDirectionX, 0);

        std::cout << "***********************************" << std::endl;
        std::cout << "**      Loads                    **" << std::endl;
        std::cout << "***********************************" << std::endl;

        int groupNodeBCRight = myStructure.GroupCreate(NuTo::Groups::eGroupId::Nodes);
        myStructure.GroupAddNodeCoordinateRange(groupNodeBCRight, 0, length - 1.0e-6, length + 1.0e-6);

        int timeDependentConstraint = myStructure.ConstraintLinearSetDisplacementNodeGroup(groupNodeBCRight, Parameters::mDirectionX, 1);

        std::cout << "***********************************" << std::endl;
        std::cout << "**      Visualization            **" << std::endl;
        std::cout << "***********************************" << std::endl;

        int visualizationGroup = myStructure.GroupCreate(NuTo::Groups::eGroupId::Elements);
        myStructure.GroupAddElementsTotal(visualizationGroup);

        myStructure.AddVisualizationComponent(visualizationGroup, NuTo::VisualizeBase::DISPLACEMENTS);
        myStructure.AddVisualizationComponent(visualizationGroup, NuTo::VisualizeBase::VELOCITY);
        myStructure.AddVisualizationComponent(visualizationGroup, NuTo::VisualizeBase::ACCELERATION);
        myStructure.AddVisualizationComponent(visualizationGroup, NuTo::VisualizeBase::ENGINEERING_STRAIN);
        myStructure.AddVisualizationComponent(visualizationGroup, NuTo::VisualizeBase::ENGINEERING_STRESS);
        myStructure.AddVisualizationComponent(visualizationGroup, NuTo::VisualizeBase::DAMAGE);

        std::cout << "***********************************" << std::endl;
        std::cout << "**      Solver                   **" << std::endl;
        std::cout << "***********************************" << std::endl;


myIntegrationScheme.SetUseLumpedMass(false);
        myStructure.NodeBuildGlobalDofs();

        myStructure.CalculateMaximumIndependentSets();

//        myIntegrationScheme.AddResultElementIpData("strain", 1, NuTo::IpData::ENGINEERING_STRAIN);
//        myIntegrationScheme.AddResultElementIpData("stress", 1, NuTo::IpData::ENGINEERING_STRESS);
//        myIntegrationScheme.AddResultElementIpData("damage", 1, NuTo::IpData::DAMAGE);


        myStructure.Info();
        //myStructure.CheckCoefficientMatrix_0(1e-6,true);

        NuTo::SparseMatrixCSRGeneral<double> massMatrix;
        NuTo::FullVector<double> someVector;

        myStructure.BuildGlobalCoefficientMatrix2(massMatrix, someVector);
        std::cout << "massMatrix.Info() \t"<< std::endl;
        massMatrix.Info();

        NuTo::FullMatrix<double, 2, 2> timeDependentLoad;
        timeDependentLoad(0, 0) = 0;
        timeDependentLoad(1, 0) = 1*Parameters::mSimulationTime;
        timeDependentLoad(0, 1) = 0;
        timeDependentLoad(1, 1) = Parameters::mLoad;

        myIntegrationScheme.SetTimeDependentConstraint(timeDependentConstraint, timeDependentLoad);

        myIntegrationScheme.Solve(Parameters::mSimulationTime);



        std::cout << "\n\n\n Results written to " + Parameters::mOutputPath.string() << std::endl;




    } catch (NuTo::MechanicsException& e)
    {
        std::cout << e.ErrorMessage();

    } catch (NuTo::MathException& e)
    {
        std::cout << e.ErrorMessage();

    }

}

