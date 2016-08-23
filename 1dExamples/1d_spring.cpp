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
#include <boost/tokenizer.hpp>
#include <sstream>

#include "nuto/mechanics/structures/unstructured/Structure.h"
#include "nuto/math/FullMatrix.h"
#include "nuto/math/SparseMatrixCSRGeneral.h"
#include "nuto/mechanics/groups/Group.h"
#include "nuto/mechanics/nodes/NodeDof.h"
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
    static const bool mAutomaticTimeStepping = true;

    static constexpr double mSpringStiffness = 10;

    static constexpr double mTimeStep = 1e-0;
    static constexpr double mMinTimeStep = 1e-5;
    static constexpr double mMaxTimeStep = 1e-0;
    static constexpr double mToleranceForce = 1e-6;
    static constexpr double mSimulationTime = 1.0;
    static constexpr double mLoadFinalDisplacement = 1.0;

    static const boost::filesystem::path mOutputPath;
    static const boost::filesystem::path mMeshFilePath;

    static const NuTo::FullVector<double, 1> mDirectionX;

};

const boost::filesystem::path Parameters::mOutputPath("/home/phuschke/1d_spring/");

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

    //**********************************************
    //          Section
    //**********************************************

    int fibreSection = myStructure.SectionCreate(NuTo::Section::SPRING);

    //**********************************************
    //          Material
    //**********************************************

    int springMaterial = myStructure.ConstitutiveLawCreate(NuTo::Constitutive::eConstitutiveType::LINEAR_SPRING);
    myStructure.ConstitutiveLawSetParameterDouble(springMaterial, NuTo::Constitutive::eConstitutiveParameter::SPRING_STIFFNESS, Parameters::mSpringStiffness);

    //**********************************************
    //          Interpolation
    //**********************************************

    int fibreInterpolationType = myStructure.InterpolationTypeCreate(NuTo::Interpolation::eShapeType::SPRING);
    myStructure.InterpolationTypeAdd(fibreInterpolationType, NuTo::Node::COORDINATES, NuTo::Interpolation::EQUIDISTANT1);
    myStructure.InterpolationTypeAdd(fibreInterpolationType, NuTo::Node::DISPLACEMENTS, NuTo::Interpolation::EQUIDISTANT1);

    myStructure.InterpolationTypeInfo(fibreInterpolationType);

    //**********************************************
    //          Spring
    //**********************************************

    std::set<NuTo::Node::eAttributes> dofs;
    dofs.insert(NuTo::Node::COORDINATES);
    dofs.insert(NuTo::Node::DISPLACEMENTS);

    NuTo::FullVector<double, -1> nodeCoords(Parameters::mDimension);
    NuTo::FullVector<int, -1> nodeIndicesSpring(2);

    // Nodes
    nodeCoords[0] = 0.0;
    int node0 = myStructure.NodeCreate(nodeCoords, dofs);

    nodeCoords[0] = 0.0;
    int node1 = myStructure.NodeCreate(nodeCoords, dofs);




    myStructure.Info();

    nodeIndicesSpring[0] = node0;
    nodeIndicesSpring[1] = node1;
    myStructure.ElementCreate(fibreInterpolationType, nodeIndicesSpring, NuTo::ElementData::eElementDataType::CONSTITUTIVELAWIP, NuTo::IpData::eIpDataType::NOIPDATA);

    //myStructure.ElementTotalConvertToInterpolationType();
    myStructure.ElementTotalSetSection(fibreSection);
    myStructure.ElementTotalSetConstitutiveLaw(springMaterial);

    myStructure.SetVerboseLevel(10);
    NuTo::FullVector<double, 1> directionX;
    directionX[0] = 1;

    std::cout << "directionX" << directionX << std::endl;

    //**********************************************
    //          Boundary Conditions
    //**********************************************

    myStructure.ConstraintLinearSetDisplacementNode(node0, directionX, 0);

    //**********************************************
    //          Loads
    //**********************************************

    myStructure.LoadCreateNodeForce(0, node1, directionX, 1);

    //**********************************************
    //          Solver
    //**********************************************

    myStructure.Info();
    myStructure.NodeBuildGlobalDofs();
    myStructure.Info();

    // build global stiffness matrix and equivalent load vector which correspond to prescribed boundary values
    NuTo::SparseMatrixCSRVector2General<double> stiffnessMatrixVec;
    NuTo::FullVector<double, Eigen::Dynamic> dispForceVector;
    myStructure.CalculateMaximumIndependentSets();
    myStructure.BuildGlobalCoefficientMatrix0(stiffnessMatrixVec, dispForceVector);
    NuTo::SparseMatrixCSRGeneral<double> stiffnessMatrix(stiffnessMatrixVec);

    stiffnessMatrix.Info();
    // build global external load vector
    NuTo::FullVector<double, Eigen::Dynamic> extForceVector;
    myStructure.BuildGlobalExternalLoadVector(0, extForceVector);
    std::cout << "extForceVector \n" << extForceVector << std::endl;
    std::cout << "dispForceVector \n" << dispForceVector << std::endl;
    // calculate right hand side
    NuTo::FullVector<double, Eigen::Dynamic> rhsVector = dispForceVector + extForceVector;

    std::cout << "rhsVector \n" << rhsVector << std::endl;
    std::cout << "?????????????????????????????????????????????????????????????????" << std::endl;
    std::cout << "stiffnessMatrix.Info() \n" << std::endl;
    stiffnessMatrix.Info();
    // solve
    NuTo::SparseDirectSolverMUMPS mySolver;
    NuTo::FullVector<double, Eigen::Dynamic> displacementVector;
    stiffnessMatrix.SetOneBasedIndexing();
    mySolver.Solve(stiffnessMatrix, rhsVector, displacementVector);

    // write displacements to node
    myStructure.NodeMergeActiveDofValues(displacementVector);

    // calculate residual
    NuTo::FullVector<double, Eigen::Dynamic> intForceVector;
    myStructure.BuildGlobalGradientInternalPotentialVector(intForceVector);
    NuTo::FullVector<double, Eigen::Dynamic> residualVector = extForceVector - intForceVector;
    std::cout << "residual: " << residualVector.Norm() << std::endl;

    // visualize results
    myStructure.AddVisualizationComponentDisplacements();
    myStructure.AddVisualizationComponentEngineeringStrain();
    myStructure.AddVisualizationComponentEngineeringStress();

    myStructure.Info();

    std::cout << "************************************************************************************" << std::endl;
    std::cout << "*                                    END                                           *" << std::endl;
    std::cout << "************************************************************************************" << std::endl;
}

