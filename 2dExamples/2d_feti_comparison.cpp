#include "nuto/mechanics/structures/unstructured/Structure.h"
#include "nuto/mechanics/timeIntegration/NewmarkDirect.h"
#include "nuto/math/SparseDirectSolverMUMPS.h"
#include <boost/filesystem.hpp>
#include <iostream>
#include <fstream>
#include <string>
#include <ctime>
#include<chrono>
constexpr unsigned int dimension = 2;

class Parameters
{
public:

    static const int mDimension = dimension;

    static const bool mPerformLineSearch = false;
    static const bool mAutomaticTimeStepping = false;

    static constexpr double mMatrixYoungsModulus = 200;
    static constexpr double mMatrixPoissonsRatio = 0.2;
    static constexpr double mMatrixThickness = 1;

    static constexpr double mMatrixNonlocalRadius = 2;
    static constexpr double mMatrixTensileStrength = 3;
    static constexpr double mMatrixCompressiveStrength = 30;
    static constexpr double mMatrixFractureEnergy = 1;

    static constexpr double mTimeStep = 1;
    static constexpr double mMinTimeStep = 1;
    static constexpr double mMaxTimeStep = 1;
    static constexpr double mToleranceForce = 1e-6;
    static constexpr double mSimulationTime = 1.0;
    static constexpr double mLoad = 1.0;

    static const boost::filesystem::path mOutputPath;
    static const boost::filesystem::path mMeshFilePathMatrix;

    static const NuTo::FullVector<double, dimension> mDirectionX;
    static const NuTo::FullVector<double, dimension> mDirectionY;

};

const boost::filesystem::path Parameters::mOutputPath("/home/phuschke/2d_feti/");
const boost::filesystem::path Parameters::mMeshFilePathMatrix("/home/phuschke/develop/nuto/myNutoExamples/MeshFiles/feti00_comparison.msh");


const NuTo::FullVector<double, dimension> Parameters::mDirectionX = NuTo::FullVector<double, dimension>::UnitX();
const NuTo::FullVector<double, dimension> Parameters::mDirectionY = NuTo::FullVector<double, dimension>::UnitY();


int main(int argc, char* argv[])
{

    Eigen::VectorXd::RealScalar bla;
    bla = 3.12;

    std::cout << "bla*bla" << bla*bla << std::endl;
    try
    {
        std::cout << "***********************************" << std::endl;
        std::cout << "**      Structure                **" << std::endl;
        std::cout << "***********************************" << std::endl;

        NuTo::Structure structure00(dimension);
        structure00.SetVerboseLevel(10);
        structure00.SetShowTime(false);
        structure00.SetNumProcessors(1);


        std::cout << "***********************************" << std::endl;
        std::cout << "**      Integration Scheme       **" << std::endl;
        std::cout << "***********************************" << std::endl;

        NuTo::NewmarkDirect myIntegrationScheme(&structure00);
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

        int section00 = structure00.SectionCreate(NuTo::Section::PLANE_STRESS);
        structure00.SectionSetThickness(section00, Parameters::mMatrixThickness);

        std::cout << "***********************************" << std::endl;
        std::cout << "**      Material                 **" << std::endl;
        std::cout << "***********************************" << std::endl;

        int material00 = structure00.ConstitutiveLawCreate(NuTo::Constitutive::eConstitutiveType::LINEAR_ELASTIC_ENGINEERING_STRESS);
        structure00.ConstitutiveLawSetParameterDouble(material00, NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, Parameters::mMatrixYoungsModulus);
        structure00.ConstitutiveLawSetParameterDouble(material00, NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO, Parameters::mMatrixPoissonsRatio);


        std::cout << "***********************************" << std::endl;
        std::cout << "**      Interpolation Type       **" << std::endl;
        std::cout << "***********************************" << std::endl;

        int interpolationType00 = structure00.InterpolationTypeCreate(NuTo::Interpolation::eShapeType::QUAD2D);
        structure00.InterpolationTypeAdd(interpolationType00, NuTo::Node::COORDINATES, NuTo::Interpolation::EQUIDISTANT1);
        structure00.InterpolationTypeAdd(interpolationType00, NuTo::Node::DISPLACEMENTS, NuTo::Interpolation::EQUIDISTANT2);


        std::cout << "***********************************" << std::endl;
        std::cout << "**      Matrix                   **" << std::endl;
        std::cout << "***********************************" << std::endl;


        boost::filesystem::path meshPath00("/home/phuschke/develop/nuto/myNutoExamples/MeshFiles/feti00_comparison.msh");



        NuTo::FullMatrix<int, Eigen::Dynamic, Eigen::Dynamic> createdGroupIdMatrix = structure00.ImportFromGmsh(meshPath00.string(), NuTo::ElementData::CONSTITUTIVELAWIP, NuTo::IpData::eIpDataType::STATICDATA);
        int subdomain00 = createdGroupIdMatrix.GetValue(0, 0);


        structure00.ElementTotalSetSection(section00);
        structure00.ElementTotalSetConstitutiveLaw(material00);
        structure00.ElementGroupSetInterpolationType(subdomain00, interpolationType00);
        structure00.ElementTotalConvertToInterpolationType(1.e-6, 10);

        std::cout << "***********************************" << std::endl;
        std::cout << "**      Boundary Conditions      **" << std::endl;
        std::cout << "***********************************" << std::endl;

        int groupNodeBCLeft = structure00.GroupCreate(NuTo::Groups::eGroupId::Nodes);
        structure00.GroupAddNodeCoordinateRange(groupNodeBCLeft, 0, -1e-6, +1e-6);

        structure00.ConstraintLinearSetDisplacementNodeGroup(groupNodeBCLeft, Parameters::mDirectionX, 0);
        structure00.ConstraintLinearSetDisplacementNodeGroup(groupNodeBCLeft, Parameters::mDirectionY, 0);


        std::cout << "***********************************" << std::endl;
        std::cout << "**      Visualization            **" << std::endl;
        std::cout << "***********************************" << std::endl;

        int groupAllElements = structure00.GroupCreate(NuTo::Groups::Elements);
        structure00.GroupAddElementsTotal(groupAllElements);
        structure00.AddVisualizationComponent(groupAllElements, NuTo::VisualizeBase::DISPLACEMENTS);
        structure00.AddVisualizationComponent(groupAllElements, NuTo::VisualizeBase::ENGINEERING_STRAIN);
        structure00.AddVisualizationComponent(groupAllElements, NuTo::VisualizeBase::ENGINEERING_STRESS);


        std::cout << "***********************************" << std::endl;
        std::cout << "**      Loads                    **" << std::endl;
        std::cout << "***********************************" << std::endl;

        NuTo::FullVector<double,2> nodeCoords;
        nodeCoords[0] = 60;
        nodeCoords[1] = 0;

        int loadNodeGroup = structure00.GroupCreate(NuTo::Groups::Nodes);
        structure00.GroupAddNodeCoordinateRange(loadNodeGroup, 0, 60 - 1e-6, 60 + 1e6);
//        structure00.GroupAddNodeRadiusRange(loadNodeGroup, nodeCoords, 0, 1.0e-6);

        structure00.LoadCreateNodeGroupForce(0,loadNodeGroup, Parameters::mDirectionY, 0.1);


        std::cout << "***********************************" << std::endl;
        std::cout << "**      Solve            **" << std::endl;
        std::cout << "***********************************" << std::endl;


        structure00.NodeBuildGlobalDofs();

        structure00.CalculateMaximumIndependentSets();


        NuTo::FullMatrix<double, 2, 2> timeDependentLoad;                                          // Gradient damage model
        timeDependentLoad(0, 0) = 0;                                                               // Gradient damage model
        timeDependentLoad(1, 0) = 1*Parameters::mSimulationTime;                                   // Gradient damage model

        timeDependentLoad(0, 1) = 0;                                                               // Gradient damage model
        timeDependentLoad(1, 1) = Parameters::mLoad;                                               // Gradient damage model

        myIntegrationScheme.SetTimeDependentLoadCase(0, timeDependentLoad);



        NuTo::SparseMatrixCSRVector2General<double> stiffnessMatrixCSRVector2;
        NuTo::FullVector<double,Eigen::Dynamic> dispForceVector, displacementVector, residualVector;

        structure00.BuildGlobalCoefficientMatrix0(stiffnessMatrixCSRVector2, dispForceVector);
        NuTo::SparseMatrixCSRGeneral<double> stiffnessMatrix (stiffnessMatrixCSRVector2);


        // solve
        NuTo::SparseDirectSolverMUMPS mySolver;
        stiffnessMatrix.SetOneBasedIndexing();



        std::clock_t begin = std::clock();
        std::chrono::time_point<std::chrono::system_clock> start, end;
        start = std::chrono::system_clock::now();
                mySolver.Solve(stiffnessMatrix, dispForceVector, displacementVector);
//        myIntegrationScheme.Solve(Parameters::mSimulationTime);
        end = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds = end-start;
        std::cout << "elapsed time: " << elapsed_seconds.count() << "s\n";

//        mySolver.Solve(stiffnessMatrix, dispForceVector, displacementVector);
        std::cout << "time to solve: " << ( std::clock() - begin ) / (double)CLOCKS_PER_SEC << std::endl;


        std::cout << "\n\n\n Results written to " + Parameters::mOutputPath.string() << std::endl;








    } catch (NuTo::MechanicsException& e)
    {
        std::cout << e.ErrorMessage();

    } catch (NuTo::MathException& e)
    {
        std::cout << e.ErrorMessage();

    }
    return 0;
}

