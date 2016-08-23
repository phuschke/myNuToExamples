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

class Parameters
{
public:

    static const int mDimension = 2;

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

    static const NuTo::FullVector<double, 2> mDirectionX;
    static const NuTo::FullVector<double, 2> mDirectionY;
};

const boost::filesystem::path Parameters::mOutputPath("/home/phuschke/2d_trusses_gmsh/");
const boost::filesystem::path Parameters::mMeshFilePath("/home/phuschke/develop/nuto/myNutoExamples/MeshFiles/2d_truss.msh");

const NuTo::FullVector<double, 2> Parameters::mDirectionX = NuTo::FullVector<double, 2>::UnitX();
const NuTo::FullVector<double, 2> Parameters::mDirectionY = NuTo::FullVector<double, 2>::UnitY();

void run()
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

    std::cout << "***********************************" << std::endl;
    std::cout << "**      Section                  **" << std::endl;
    std::cout << "***********************************" << std::endl;


    int fibreSection = myStructure.SectionCreate(NuTo::Section::TRUSS);
    myStructure.SectionSetArea(fibreSection, Parameters::mFibreCrossSection);

    std::cout << "***********************************" << std::endl;
    std::cout << "**      Material                 **" << std::endl;
    std::cout << "***********************************" << std::endl;

    int fibreMaterial = myStructure.ConstitutiveLawCreate(NuTo::Constitutive::eConstitutiveType::LINEAR_ELASTIC_ENGINEERING_STRESS);
    myStructure.ConstitutiveLawSetParameterDouble(fibreMaterial, NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, Parameters::mFibreYoungsModulus);
    myStructure.ConstitutiveLawSetParameterDouble(fibreMaterial, NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO, Parameters::mFibrePoissonsRatio);



    std::cout << "***********************************" << std::endl;
    std::cout << "**      Interpolation Type       **" << std::endl;
    std::cout << "***********************************" << std::endl;

    int fibreInterpolationType = myStructure.InterpolationTypeCreate(NuTo::Interpolation::eShapeType::TRUSSXD);
    myStructure.InterpolationTypeAdd(fibreInterpolationType, NuTo::Node::COORDINATES, NuTo::Interpolation::EQUIDISTANT2);
    myStructure.InterpolationTypeAdd(fibreInterpolationType, NuTo::Node::DISPLACEMENTS, NuTo::Interpolation::EQUIDISTANT2);


    std::cout << "***********************************" << std::endl;
    std::cout << "**      Matrix                   **" << std::endl;
    std::cout << "***********************************" << std::endl;

    NuTo::FullMatrix<int, Eigen::Dynamic, Eigen::Dynamic> createdGroupIds = myStructure.ImportFromGmsh(Parameters::mMeshFilePath.string(), NuTo::ElementData::CONSTITUTIVELAWIP, NuTo::IpData::eIpDataType::STATICDATA);
    int groupIdTruss = createdGroupIds.GetValue(0, 0);

    myStructure.ElementGroupSetInterpolationType(groupIdTruss, fibreInterpolationType);
    myStructure.ElementTotalConvertToInterpolationType(1.e-6, 10);
    myStructure.ElementGroupSetSection(groupIdTruss, fibreSection);
    myStructure.ElementGroupSetConstitutiveLaw(groupIdTruss, fibreMaterial);

    myStructure.Info();

    myStructure.ElementTotalConvertToInterpolationType(1.e-6, 10);
    myStructure.NodeBuildGlobalDofs();
    myStructure.Info();
    myStructure.CheckCoefficientMatrix_0(1e-6, true);

    std::cout << "***********************************" << std::endl;
    std::cout << "**      Boundary Conditions      **" << std::endl;
    std::cout << "***********************************" << std::endl;

    NuTo::FullVector<double, 2> directionAligned;
    directionAligned << 1, 1;
    NuTo::FullVector<double, 2> directionOrthogonal;
    directionOrthogonal << -1, 1;

    NuTo::FullVector<double, 2> nodeCoords;
    nodeCoords[0] = 0;
    nodeCoords[1] = 0;
    int groupNodeAll = myStructure.GroupCreate(NuTo::Groups::eGroupId::Nodes);
    myStructure.GroupAddNodeRadiusRange(groupNodeAll, nodeCoords, 0, 20);


    myStructure.ConstraintLinearSetDisplacementNode(0, directionAligned,0);


    myStructure.ConstraintLinearSetDisplacementNodeGroup(groupNodeAll, directionOrthogonal,0);

    std::cout << "***********************************" << std::endl;
    std::cout << "**      Loads                    **" << std::endl;
    std::cout << "***********************************" << std::endl;

    nodeCoords[0] = 10;
    nodeCoords[1] = 10;
    int groupNodeLoad = myStructure.GroupCreate(NuTo::Groups::eGroupId::Nodes);
    myStructure.GroupAddNodeRadiusRange(groupNodeLoad, nodeCoords, 0, 1e-6);

    myStructure.LoadCreateNodeGroupForce(0, groupNodeLoad, directionAligned, 1e7);



    std::cout << "***********************************" << std::endl;
    std::cout << "**      Visualization            **" << std::endl;
    std::cout << "***********************************" << std::endl;

    myStructure.AddVisualizationComponentDisplacements();
    myStructure.AddVisualizationComponentEngineeringStrain();
    myStructure.AddVisualizationComponentEngineeringStress();
    myStructure.AddVisualizationComponentSection();
    myStructure.AddVisualizationComponentConstitutive();

    std::cout << "***********************************" << std::endl;
    std::cout << "**      Solver                   **" << std::endl;
    std::cout << "***********************************" << std::endl;

    myStructure.Info();
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

    myIntegrationScheme.AddResultElementIpStress("stress", 0);
    myIntegrationScheme.AddResultElementIpValue("strain", 0);

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
