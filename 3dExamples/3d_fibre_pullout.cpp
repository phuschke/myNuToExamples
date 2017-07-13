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

//**********************************************
//          Parameter Initialization
//**********************************************

class ParametersMaterial
{
public:
    static constexpr double mYoungsModulus = 10;
    static constexpr double mPoissonsRatio = 0.0;
};

class ParametersTimeIntegration
{
public:
    static constexpr double mTimeStep = 1.;
    static constexpr double mSimulationTime = 1.0;
    static constexpr double mLoad = 1.0;
};

class ParametersGeometry3D
{
public:
    static constexpr int mDimension = 3;
    static constexpr double mCrossSection = 0.1;
    static const NuTo::FullVector<double, 3> mDirectionX;
    static const NuTo::FullVector<double, 3> mDirectionY;
    static const NuTo::FullVector<double, 3> mDirectionZ;
};

const NuTo::FullVector<double, 3> ParametersGeometry3D::mDirectionX = NuTo::FullVector<double, 3>::UnitX();
const NuTo::FullVector<double, 3> ParametersGeometry3D::mDirectionY = NuTo::FullVector<double, 3>::UnitY();
const NuTo::FullVector<double, 3> ParametersGeometry3D::mDirectionZ = NuTo::FullVector<double, 3>::UnitZ();

//**********************************************
//          Run3d
//**********************************************

void Run3d(NuTo::FullVector<double, -1> rNodeCoords0, NuTo::FullVector<double, -1> rNodeCoords1,
           NuTo::FullVector<double, -1> rNodeCoords2, NuTo::FullVector<double, -1> rDirectionAligned,
           NuTo::FullVector<double, -1> rDirectionOrthogonal0, NuTo::FullVector<double, -1> rDirectionOrthogonal1)
{
    //**********************************************
    //          Structure
    //**********************************************

    NuTo::Structure myStructure(ParametersGeometry3D::mDimension);
    myStructure.SetVerboseLevel(10);
    //**********************************************
    //         Integration Scheme
    //**********************************************

    boost::filesystem::path resultPath(boost::filesystem::initial_path().string() +
                                       "/resultElementsInHigherDimensions/");
    boost::filesystem::remove_all(resultPath);
    boost::filesystem::create_directory(resultPath);

    NuTo::NewmarkDirect myIntegrationScheme(&myStructure);
    myIntegrationScheme.SetTimeStep(ParametersTimeIntegration::mTimeStep);
    myIntegrationScheme.SetResultDirectory(resultPath.string(), false);

    //**********************************************
    //          Section
    //**********************************************

    int fibreSection = myStructure.SectionCreate(NuTo::Section::TRUSS);
    myStructure.SectionSetArea(fibreSection, ParametersGeometry3D::mCrossSection);

    //**********************************************
    //          Material
    //**********************************************

    int fibreMaterial =
            myStructure.ConstitutiveLawCreate(NuTo::Constitutive::eConstitutiveType::LINEAR_ELASTIC_ENGINEERING_STRESS);
    myStructure.ConstitutiveLawSetYoungsModulus(fibreMaterial, ParametersMaterial::mYoungsModulus);
    myStructure.ConstitutiveLawSetPoissonsRatio(fibreMaterial, ParametersMaterial::mPoissonsRatio);

    //**********************************************
    //          Interpolation
    //**********************************************

    int fibreInterpolationType = myStructure.InterpolationTypeCreate(NuTo::Interpolation::eShapeType::TRUSSXD);
    myStructure.InterpolationTypeAdd(fibreInterpolationType, NuTo::Node::COORDINATES,
                                     NuTo::Interpolation::EQUIDISTANT2);
    myStructure.InterpolationTypeAdd(fibreInterpolationType, NuTo::Node::DISPLACEMENTS,
                                     NuTo::Interpolation::EQUIDISTANT2);
    myStructure.InterpolationTypeSetIntegrationType(
            fibreInterpolationType, NuTo::IntegrationType::IntegrationType1D2NGauss2Ip, NuTo::IpData::STATICDATA);

    //**********************************************
    //          Fibre
    //**********************************************

    // Nodes
    int node0 = myStructure.NodeCreate(rNodeCoords0);
    int node1 = myStructure.NodeCreate(rNodeCoords1);
    int node2 = myStructure.NodeCreate(rNodeCoords2);

    NuTo::FullVector<int, -1> nodeIndices(3);
    nodeIndices[0] = node0;
    nodeIndices[1] = node1;
    nodeIndices[2] = node2;
    myStructure.ElementCreate(fibreInterpolationType, nodeIndices,
                              NuTo::ElementData::eElementDataType::CONSTITUTIVELAWIP,
                              NuTo::IpData::eIpDataType::NOIPDATA);

    myStructure.ElementTotalConvertToInterpolationType(1e-6, 10);
    myStructure.ElementTotalSetSection(fibreSection);
    myStructure.ElementTotalSetConstitutiveLaw(fibreMaterial);

    //**********************************************
    //          Boundary Conditions
    //**********************************************
    myStructure.ConstraintLinearSetDisplacementNode(node0, ParametersGeometry3D::mDirectionX, 0);
    myStructure.ConstraintLinearSetDisplacementNode(node0, ParametersGeometry3D::mDirectionY, 0);
    myStructure.ConstraintLinearSetDisplacementNode(node0, ParametersGeometry3D::mDirectionZ, 0);

    myStructure.ConstraintLinearSetDisplacementNode(node1, rDirectionOrthogonal0, 0);
    myStructure.ConstraintLinearSetDisplacementNode(node2, rDirectionOrthogonal0, 0);

    myStructure.ConstraintLinearSetDisplacementNode(node1, rDirectionOrthogonal1, 0);
    myStructure.ConstraintLinearSetDisplacementNode(node2, rDirectionOrthogonal1, 0);

    //**********************************************
    //          Loads
    //**********************************************

    int load = myStructure.LoadCreateNodeForce(0, node2, rDirectionAligned, 1);

    NuTo::FullMatrix<double, 2, 2> timeDependentLoad;
    timeDependentLoad(0, 0) = 0;
    timeDependentLoad(1, 0) = ParametersTimeIntegration::mSimulationTime;
    timeDependentLoad(0, 1) = 0;
    timeDependentLoad(1, 1) = ParametersTimeIntegration::mLoad;

    myIntegrationScheme.SetTimeDependentLoadCase(load, timeDependentLoad);

    //**********************************************
    //          Solver
    //**********************************************

    myStructure.NodeBuildGlobalDofs();
    myStructure.CalculateMaximumIndependentSets();

    myStructure.CheckCoefficientMatrix_0(1e-6, true);
    myStructure.Info();

    myIntegrationScheme.Solve(ParametersTimeIntegration::mSimulationTime);

    myStructure.Info();

    NuTo::FullVector<double, -1> displacements;
    boost::filesystem::remove_all(resultPath);

    myStructure.NodeGetDisplacements(node2, displacements);

    if (std::abs(displacements.norm() - std::sqrt(3)) > 1e-8)
        throw NuTo::Exception("The calculated displacements do not agree with the analytical solution!");
}

int main(int argc, char* argv[])
{
    // node coordinates
    NuTo::FullVector<double, ParametersGeometry3D::mDimension> nodeCoords0;
    NuTo::FullVector<double, ParametersGeometry3D::mDimension> nodeCoords1;
    NuTo::FullVector<double, ParametersGeometry3D::mDimension> nodeCoords2;

    // directions
    NuTo::FullVector<double, ParametersGeometry3D::mDimension> directionAligned;
    NuTo::FullVector<double, ParametersGeometry3D::mDimension> directionOrthogonal0;
    NuTo::FullVector<double, ParametersGeometry3D::mDimension> directionOrthogonal1;

    try
    {
        //
        //  example 1: quadratic truss element, 0° inclined
        //
        std::cout << "\n \n *************** Start Example 1 ************************* \n \n " << std::endl;
        nodeCoords0[0] = 1.0;
        nodeCoords0[1] = 5.0;
        nodeCoords0[2] = 3.0;

        nodeCoords1[0] = 1.0 + 0.5 * std::sqrt(3);
        nodeCoords1[1] = 5.0;
        nodeCoords1[2] = 3.0;

        nodeCoords2[0] = 1.0 + 1.0 * std::sqrt(3);
        nodeCoords2[1] = 5.0;
        nodeCoords2[2] = 3.0;

        directionAligned[0] = 1.0;
        directionAligned[1] = 0.0;
        directionAligned[2] = 0.0;

        directionOrthogonal0[0] = 0.0;
        directionOrthogonal0[1] = 1.0;
        directionOrthogonal0[2] = 0.0;

        directionOrthogonal1[0] = 0.0;
        directionOrthogonal1[1] = 0.0;
        directionOrthogonal1[2] = 1.0;

        Run3d(nodeCoords0, nodeCoords1, nodeCoords2, directionAligned, directionOrthogonal0, directionOrthogonal1);

        //
        //  example 2: quadratic truss element, 90° inclined
        //

        std::cout << "\n \n *************** Start Example 2 ************************* \n \n " << std::endl;
        nodeCoords0[0] = 1.0;
        nodeCoords0[1] = 5.0;
        nodeCoords0[2] = 3.0;

        nodeCoords1[0] = 1.0;
        nodeCoords1[1] = 5.0;
        nodeCoords1[2] = 3.0 + 0.5 * std::sqrt(3);

        nodeCoords2[0] = 1.0;
        nodeCoords2[1] = 5.0;
        nodeCoords2[2] = 3.0 + 1.0 * std::sqrt(3);

        directionAligned[0] = 0.0;
        directionAligned[1] = 0.0;
        directionAligned[2] = 1.0;

        directionOrthogonal0[0] = 1.0;
        directionOrthogonal0[1] = 1.0;
        directionOrthogonal0[2] = 0.0;

        directionOrthogonal1[0] = -1.0;
        directionOrthogonal1[1] = 1.0;
        directionOrthogonal1[2] = 0.0;

        Run3d(nodeCoords0, nodeCoords1, nodeCoords2, directionAligned, directionOrthogonal0, directionOrthogonal1);

        //
        //  example 3: quadratic truss element, 45° inclined
        //

        std::cout << "\n \n *************** Start Example 3 ************************* \n \n " << std::endl;
        nodeCoords0[0] = 1.0;
        nodeCoords0[1] = 5.0;
        nodeCoords0[2] = 3.0;

        nodeCoords1[0] = 1.0 + 0.5;
        nodeCoords1[1] = 5.0 + 0.5;
        nodeCoords1[2] = 3.0 + 0.5;

        nodeCoords2[0] = 1.0 + 1.0;
        nodeCoords2[1] = 5.0 + 1.0;
        nodeCoords2[2] = 3.0 + 1.0;

        directionAligned[0] = 1.0;
        directionAligned[1] = 1.0;
        directionAligned[2] = 1.0;
        directionAligned.normalize();

        directionOrthogonal0[0] = -1.0;
        directionOrthogonal0[1] = -1.0;
        directionOrthogonal0[2] = 2.0;
        directionOrthogonal0.normalize();

        directionOrthogonal1[0] = 2.0;
        directionOrthogonal1[1] = -1.0;
        directionOrthogonal1[2] = -1.0;
        directionOrthogonal1.normalize();


        Run3d(nodeCoords0, nodeCoords1, nodeCoords2, directionAligned, directionOrthogonal0, directionOrthogonal1);
    }
    catch (NuTo::Exception& e)
    {
        std::cout << "## Test failed ##" << std::endl;
        std::cout << e.ErrorMessage();
        return -1;
    }

    std::cout << "## Test successful ##" << std::endl;
}
