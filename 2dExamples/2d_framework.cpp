#include "nuto/math/FullMatrix.h"
#include "nuto/mechanics/structures/unstructured/Structure.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStrain3D.h"
#include "nuto/mechanics/constitutive/mechanics/EngineeringStress3D.h"
#include "nuto/mechanics/constitutive/mechanics/GradientDamageEngineeringStress.h"
#include "nuto/mechanics/MechanicsException.h"
#include "nuto/mechanics/constitutive/ConstitutiveStaticDataBase.h"
#include "nuto/mechanics/constitutive/mechanics/ConstitutiveStaticDataGradientDamage1D.h"
#include "nuto/math/SparseMatrixCSRGeneral.h"
#include "nuto/math/SparseDirectSolverMUMPS.h"
#include "nuto/math/SparseMatrixCSRVector2General.h"

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

    static constexpr double mFibreYoungsModulus = 10;
    static constexpr double mFibrePoissonsRatio = 0.0;
    static constexpr double mFibreCrossSection = 0.1;

    static constexpr double mTimeStep = 1e-0;
    static constexpr double mMinTimeStep = 1e-5;
    static constexpr double mMaxTimeStep = 1e-0;
    static constexpr double mToleranceForce = 1e-6;
    static constexpr double mSimulationTime = 1.0;
    static constexpr double mLoadFinalDisplacement = 1.0;

    static const boost::filesystem::path mOutputPath;
    static const boost::filesystem::path mMeshFilePath;

    static const NuTo::FullVector<double, 2> mDirectionX;
    static const NuTo::FullVector<double, 2> mDirectionY;
};

const boost::filesystem::path Parameters::mOutputPath("/home/phuschke/2D_Framework/");

const NuTo::FullVector<double, 2> Parameters::mDirectionX({1, 0});
const NuTo::FullVector<double, 2> Parameters::mDirectionY({0, 1});

int main(int argc, char* argv[])
{

    boost::filesystem::remove_all(Parameters::mOutputPath);
    boost::filesystem::create_directory(Parameters::mOutputPath);

    //**********************************************stiffnessMatrix.Info()
    //          Structure
    //**********************************************

    NuTo::Structure myStructure(Parameters::mDimension);

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
    myIntegrationScheme.SetResultDirectory(Parameters::mOutputPath.string(), false);

    //**********************************************
    //          Section
    //**********************************************

    int fibreSection = myStructure.SectionCreate(NuTo::Section::TRUSS);
    myStructure.SectionSetArea(fibreSection, Parameters::mFibreCrossSection);

    //**********************************************
    //          Material
    //**********************************************

    int fibreMaterial =
            myStructure.ConstitutiveLawCreate(NuTo::Constitutive::eConstitutiveType::LINEAR_ELASTIC_ENGINEERING_STRESS);
    myStructure.ConstitutiveLawSetYoungsModulus(fibreMaterial, Parameters::mFibreYoungsModulus);
    myStructure.ConstitutiveLawSetPoissonsRatio(fibreMaterial, Parameters::mFibrePoissonsRatio);

    //**********************************************
    //          Interpolation
    //**********************************************

    int fibreInterpolationType = myStructure.InterpolationTypeCreate(NuTo::Interpolation::eShapeType::TRUSSXD);
    myStructure.InterpolationTypeAdd(fibreInterpolationType, NuTo::Node::COORDINATES,
                                     NuTo::Interpolation::EQUIDISTANT1);
    myStructure.InterpolationTypeAdd(fibreInterpolationType, NuTo::Node::DISPLACEMENTS,
                                     NuTo::Interpolation::EQUIDISTANT1);
    myStructure.InterpolationTypeSetIntegrationType(
            fibreInterpolationType, NuTo::IntegrationType::IntegrationType1D2NGauss1Ip, NuTo::IpData::STATICDATA);

    myStructure.InterpolationTypeInfo(fibreInterpolationType);

    //**********************************************
    //          Fibre
    //**********************************************

    NuTo::FullVector<double, -1> nodeCoords(Parameters::mDimension);
    NuTo::FullVector<int, -1> nodeIndicesFibre(2);

    // Nodes
    nodeCoords[0] = 1.0;
    nodeCoords[1] = 1.0;
    int node0 = myStructure.NodeCreate(nodeCoords);

    nodeCoords[0] = 2.0;
    nodeCoords[1] = 1.0;
    int node1 = myStructure.NodeCreate(nodeCoords);

    nodeCoords[0] = 2.0;
    nodeCoords[1] = 2.0;
    int node2 = myStructure.NodeCreate(nodeCoords);

    nodeIndicesFibre[0] = node0;
    nodeIndicesFibre[1] = node1;
    myStructure.ElementCreate(fibreInterpolationType, nodeIndicesFibre,
                              NuTo::ElementData::eElementDataType::CONSTITUTIVELAWIP,
                              NuTo::IpData::eIpDataType::NOIPDATA);

    nodeIndicesFibre[0] = node1;
    nodeIndicesFibre[1] = node2;
    myStructure.ElementCreate(fibreInterpolationType, nodeIndicesFibre,
                              NuTo::ElementData::eElementDataType::CONSTITUTIVELAWIP,
                              NuTo::IpData::eIpDataType::NOIPDATA);

    nodeIndicesFibre[0] = node2;
    nodeIndicesFibre[1] = node0;
    myStructure.ElementCreate(fibreInterpolationType, nodeIndicesFibre,
                              NuTo::ElementData::eElementDataType::CONSTITUTIVELAWIP,
                              NuTo::IpData::eIpDataType::NOIPDATA);

    myStructure.SetVerboseLevel(10);
    myStructure.Info();

    myStructure.ElementTotalConvertToInterpolationType(1e-6, 10);
    myStructure.ElementTotalSetSection(fibreSection);
    myStructure.ElementTotalSetConstitutiveLaw(fibreMaterial);


    //**********************************************
    //          Boundary Conditions
    //**********************************************
    myStructure.ConstraintLinearSetDisplacementNode(node0, Parameters::mDirectionX, 0);
    myStructure.ConstraintLinearSetDisplacementNode(node0, Parameters::mDirectionY, 0);

    myStructure.ConstraintLinearSetDisplacementNode(node1, Parameters::mDirectionX, 0);
    myStructure.ConstraintLinearSetDisplacementNode(node1, Parameters::mDirectionY, 0);

    //**********************************************
    //          Loads
    //**********************************************
    NuTo::FullMatrix<double, 2, 2> dispRHS;
    dispRHS(0, 0) = 0;
    dispRHS(1, 0) = Parameters::mSimulationTime;
    dispRHS(0, 1) = 0;
    dispRHS(1, 1) = Parameters::mLoadFinalDisplacement;


    // myStructure.ConstraintLinearSetDisplacementNode(node2, Parameters::mDirectionX, 0);
    // myStructure.ConstraintLinearSetDisplacementNode(node2, Parameters::mDirectionY, 3);

    int load = 0;
    bool enableDisplacementControl = true;
    if (enableDisplacementControl)
    {
        load = myStructure.ConstraintLinearSetDisplacementNode(node2, Parameters::mDirectionY, 1);
        myIntegrationScheme.SetTimeDependentConstraint(load, dispRHS);
    }
    else
    {
        load = myStructure.LoadCreateNodeForce(0, node2, Parameters::mDirectionY, 1);
        myIntegrationScheme.SetTimeDependentLoadCase(load, dispRHS);
    }


    //**********************************************
    //          Visualisation
    //**********************************************

    myStructure.AddVisualizationComponentDisplacements();

    //**********************************************
    //          Solver
    //**********************************************

    myStructure.NodeBuildGlobalDofs();

    myStructure.ElementCheckCoefficientMatrix_0(1e-6);

    myIntegrationScheme.Solve(Parameters::mSimulationTime);
    myStructure.Info();

    myStructure.ExportVtkDataFileElements(Parameters::mOutputPath.string() + "results.vtk");
    std::cout << "************************************************************************************" << std::endl;
    std::cout << "*                                    END                                           *" << std::endl;
    std::cout << "************************************************************************************" << std::endl;
}
