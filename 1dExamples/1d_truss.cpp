

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

const unsigned int dimension = 1;

class Parameters
{
public:

    static const int mDimension = 1;

    static const bool mPerformLineSearch = true;
    static const bool mAutomaticTimeStepping = false;

    static constexpr double mFibreYoungsModulus = 20;
    static constexpr double mFibrePoissonsRatio = 0.0;
    static constexpr double mFibreCrossSection = 0.1;

    static constexpr double mTimeStep = 1e-1;
    static constexpr double mMinTimeStep = 1e-5;
    static constexpr double mMaxTimeStep = 1e-0;
    static constexpr double mToleranceForce = 1e-6;
    static constexpr double mSimulationTime = 1.0;
    static constexpr double mLoadFinalDisplacement = 1.0;

    static const boost::filesystem::path mOutputPath;
    static const boost::filesystem::path mMeshFilePath;

    static const NuTo::FullVector<double, 1> mDirectionX;
    static const NuTo::FullVector<double, 1> mDirectionY;
};

const boost::filesystem::path Parameters::mOutputPath("/home/phuschke/1d_truss/");

int main(int argc, char* argv[])
{

    //**********************************************
    //          Structure
    //**********************************************

    NuTo::Structure myStructure(Parameters::mDimension);
    myStructure.SetVerboseLevel(10);
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
    myIntegrationScheme.SetResultDirectory(Parameters::mOutputPath.string(), true);

    //**********************************************
    //          Section
    //**********************************************

    int fibreSection = myStructure.SectionCreate(NuTo::Section::TRUSS);
    myStructure.SectionSetArea(fibreSection, Parameters::mFibreCrossSection);

    //**********************************************
    //          Material
    //**********************************************

    int fibreMaterial = myStructure.ConstitutiveLawCreate(NuTo::Constitutive::eConstitutiveType::LINEAR_ELASTIC_ENGINEERING_STRESS);
    myStructure.ConstitutiveLawSetParameterDouble(fibreMaterial,NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS ,Parameters::mFibreYoungsModulus);
    myStructure.ConstitutiveLawSetParameterDouble(fibreMaterial,NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO ,Parameters::mFibrePoissonsRatio);

    //**********************************************
    //          Interpolation
    //**********************************************

    int fibreInterpolationType = myStructure.InterpolationTypeCreate(NuTo::Interpolation::eShapeType::TRUSS1D);
    myStructure.InterpolationTypeAdd(fibreInterpolationType, NuTo::Node::COORDINATES, NuTo::Interpolation::EQUIDISTANT1);
    myStructure.InterpolationTypeAdd(fibreInterpolationType, NuTo::Node::DISPLACEMENTS, NuTo::Interpolation::EQUIDISTANT1);
    myStructure.InterpolationTypeSetIntegrationType(fibreInterpolationType, NuTo::IntegrationType::IntegrationType1D2NGauss1Ip, NuTo::IpData::STATICDATA);

    myStructure.InterpolationTypeInfo(fibreInterpolationType);

    //**********************************************
    //          Fibre
    //**********************************************

    NuTo::FullVector<double, -1> nodeCoords(Parameters::mDimension);
    NuTo::FullVector<int, -1> nodeIndicesFibre(2);

    // Nodes
    nodeCoords[0] = 0.0;
    int node0 = myStructure.NodeCreate(nodeCoords);

    nodeCoords[0] = 999;
    int node1 = myStructure.NodeCreate(nodeCoords);

    nodeIndicesFibre[0] = node0;
    nodeIndicesFibre[1] = node1;
    myStructure.ElementCreate(fibreInterpolationType, nodeIndicesFibre, NuTo::ElementData::eElementDataType::CONSTITUTIVELAWIP, NuTo::IpData::eIpDataType::NOIPDATA);

    myStructure.ElementTotalConvertToInterpolationType(1e-6, 10);
    myStructure.ElementTotalSetSection(fibreSection);
    myStructure.ElementTotalSetConstitutiveLaw(fibreMaterial);

    NuTo::FullVector<double, 1> directionX;
    directionX[0] = 1;

    std::cout << "directionX"<< directionX << std::endl;

    myStructure.ConstraintLinearSetDisplacementNode(node0, directionX,0);

    int loadCase = 0;
    myStructure.LoadCreateNodeForce(loadCase, node1, directionX, 1);

    myStructure.NodeBuildGlobalDofs();
    myStructure.Info();
    NuTo::FullMatrix<double, 2, 2> timeDependentLoad;
    timeDependentLoad(0, 0) = 0;
    timeDependentLoad(1, 0) = 1;
    timeDependentLoad(0, 1) = 0;
    timeDependentLoad(1, 1) = 1;

    myIntegrationScheme.SetTimeDependentLoadCase(loadCase, timeDependentLoad);
    myIntegrationScheme.Solve(Parameters::mSimulationTime);

    myStructure.Info();

    std::cout << "************************************************************************************" << std::endl;
    std::cout << "*                                    END                                           *" << std::endl;
    std::cout << "************************************************************************************" << std::endl;
}

