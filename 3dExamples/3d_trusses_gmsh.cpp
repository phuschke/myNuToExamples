#include "nuto/math/FullMatrix.h"
#include "nuto/mechanics/structures/unstructured/Structure.h"

#include "nuto/mechanics/MechanicsException.h"

#include "nuto/mechanics/timeIntegration/NewmarkDirect.h"
#include "nuto/base/Logger.h"

#include <boost-1_55/boost/filesystem.hpp>
#include <boost-1_55/boost/lexical_cast.hpp>
#include <boost-1_55/boost/progress.hpp>

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
    static constexpr double mMatrixThickness = 10;

    static constexpr double mFibreYoungsModulus = 2.1e9;   // steel
    static constexpr double mFibrePoissonsRatio = 0.2;
    static constexpr double mFibreCrossSection = 0.01;

    static constexpr double mNonlocalRadius = 4;
    static constexpr double mTensileStrength = 3;
    static constexpr double mCompressiveStrength = 30;
    static constexpr double mFractureEnergy = 0.01;
    static constexpr double mTimeStep = 1e-3;
    static constexpr double mMinTimeStep = 1e-5;
    static constexpr double mMaxTimeStep = 1e-1;
    static constexpr double mToleranceForce = 1e-6;
    static constexpr double mSimulationTime = 1.0;
    static constexpr double mLoad = 2.0;

    static const boost::filesystem::path mOutputPath;
    static const boost::filesystem::path mMeshFilePath;

    static const NuTo::FullVector<double, dimension> mDirectionX;
    static const NuTo::FullVector<double, dimension> mDirectionY;
    static const NuTo::FullVector<double, dimension> mDirectionZ;
};

const boost::filesystem::path Parameters::mOutputPath("/home/phuschke/3d_trusses_gmsh/");
const boost::filesystem::path Parameters::mMeshFilePath("/home/phuschke/develop/nuto/myNutoExamples/3d_truss.msh");

const NuTo::FullVector<double, dimension> Parameters::mDirectionX = NuTo::FullVector<double, dimension>::UnitX();
const NuTo::FullVector<double, dimension> Parameters::mDirectionY = NuTo::FullVector<double, dimension>::UnitY();
const NuTo::FullVector<double, dimension> Parameters::mDirectionZ = NuTo::FullVector<double, dimension>::UnitZ();

void run()
{

    //**********************************************
    //          Structure
    //**********************************************

    NuTo::Structure myStructure(Parameters::mDimension);
    myStructure.SetVerboseLevel(10);
    myStructure.SetShowTime(false);

    //**********************************************
    //         Integration Scheme
    //**********************************************

    NuTo::NewmarkDirect myIntegrationScheme(&myStructure);
    myIntegrationScheme.SetTimeStep(Parameters::mTimeStep);
    myIntegrationScheme.SetMinTimeStep(Parameters::mMinTimeStep);
    myIntegrationScheme.SetMaxTimeStep(Parameters::mMaxTimeStep);
    myIntegrationScheme.SetToleranceForce(Parameters::mToleranceForce);
    myIntegrationScheme.SetAutomaticTimeStepping(Parameters::mAutomaticTimeStepping);
    myIntegrationScheme.SetPerformLineSearch(Parameters::mPerformLineSearch);

    //**********************************************
    //          Section
    //**********************************************


    int fibreSection = myStructure.SectionCreate(NuTo::Section::TRUSS);
    myStructure.SectionSetArea(fibreSection, Parameters::mFibreCrossSection);

    //**********************************************
    //          Material
    //**********************************************
    int fibreMaterial = myStructure.ConstitutiveLawCreate(NuTo::Constitutive::eConstitutiveType::LINEAR_ELASTIC_ENGINEERING_STRESS);
    myStructure.ConstitutiveLawSetParameterDouble(fibreMaterial, NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, Parameters::mFibreYoungsModulus);
    myStructure.ConstitutiveLawSetParameterDouble(fibreMaterial, NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO, Parameters::mFibrePoissonsRatio);

    //**********************************************
    //          Interpolationtype
    //**********************************************

    int fibreInterpolationType = myStructure.InterpolationTypeCreate(NuTo::Interpolation::eShapeType::TRUSSXD);
    myStructure.InterpolationTypeAdd(fibreInterpolationType, NuTo::Node::COORDINATES, NuTo::Interpolation::EQUIDISTANT1);
    myStructure.InterpolationTypeAdd(fibreInterpolationType, NuTo::Node::DISPLACEMENTS, NuTo::Interpolation::EQUIDISTANT1);


    //**********************************************
    //          Matrix Geometry
    //**********************************************

    NuTo::FullMatrix<int, Eigen::Dynamic, Eigen::Dynamic> createdGroupIds = myStructure.ImportFromGmsh(Parameters::mMeshFilePath.string(), NuTo::ElementData::CONSTITUTIVELAWIP, NuTo::IpData::eIpDataType::STATICDATA);
    int groupIdMatrix = createdGroupIds.GetValue(0, 0);

    myStructure.ElementGroupSetInterpolationType(groupIdMatrix, fibreInterpolationType);
    myStructure.ElementGroupSetSection(groupIdMatrix, fibreSection);
    myStructure.ElementGroupSetConstitutiveLaw(groupIdMatrix, fibreMaterial);

    myStructure.Info();

    myStructure.ElementTotalConvertToInterpolationType(1.e-6, 10);
    myStructure.NodeBuildGlobalDofs();
    myStructure.Info();
    myStructure.CheckCoefficientMatrix_0(1e-6, true);

    //**********************************************
    //          Boundary Conditions
    //**********************************************

    NuTo::FullVector<double, dimension> directionAligned;
    NuTo::FullVector<double, dimension> directionOrthogonal0;
    NuTo::FullVector<double, dimension> directionOrthogonal1;


    directionAligned << 1, 1, 1;
    directionOrthogonal0 << -1, -1, 2;
    directionOrthogonal1 << 2, -1, -1;

    myStructure.ConstraintLinearSetDisplacementNode(0, Parameters::mDirectionX,0);
    myStructure.ConstraintLinearSetDisplacementNode(0, Parameters::mDirectionY,0);
    myStructure.ConstraintLinearSetDisplacementNode(0, Parameters::mDirectionZ,0);

    myStructure.ConstraintLinearSetDisplacementNode(1, directionOrthogonal0,0);
    myStructure.ConstraintLinearSetDisplacementNode(1, directionOrthogonal1,0);

    //**********************************************
    //          Loads
    //**********************************************


    myStructure.LoadCreateNodeForce(0, 1, directionAligned, 1e7);



    //**********************************************
    //          Visualisation
    //**********************************************

    myStructure.AddVisualizationComponentDisplacements();
    myStructure.AddVisualizationComponentEngineeringStrain();
    myStructure.AddVisualizationComponentEngineeringStress();
    myStructure.AddVisualizationComponentSection();
    myStructure.AddVisualizationComponentConstitutive();

    //**********************************************
    //          Solver
    //**********************************************

    myStructure.NodeBuildGlobalDofs();
    myStructure.CheckCoefficientMatrix_0(1e-6, true);
    myStructure.Info();

    boost::filesystem::remove_all(Parameters::mOutputPath);
    boost::filesystem::create_directory(Parameters::mOutputPath);

    myStructure.NodeBuildGlobalDofs();
    myStructure.CalculateMaximumIndependentSets();

    NuTo::FullMatrix<double, 2, 2> timeDependentLoad;
    timeDependentLoad(0, 0) = 0;
    timeDependentLoad(1, 0) = Parameters::mSimulationTime;
    timeDependentLoad(0, 1) = 0;
    timeDependentLoad(1, 1) = Parameters::mLoad;

    //myIntegrationScheme.SetTimeDependentLoadCase(loadCase, timeDependentLoad);
    myIntegrationScheme.SetTimeDependentLoadCase(0, timeDependentLoad);

    myIntegrationScheme.SetResultDirectory(Parameters::mOutputPath.string(), false);

    myIntegrationScheme.Solve(Parameters::mSimulationTime);



    myStructure.Info();


}



int main()
{
    try
    {
     run();
    } catch (NuTo::MechanicsException& e)
    {
        std::cout << e.ErrorMessage();
        return -1;
    } catch (...)
    {
        std::cout << "Something else went wrong." << std::endl;
        return -1;
    }

    std::cout << " ===> End <=== " << std::endl;
}
