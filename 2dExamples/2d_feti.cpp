#include "nuto/mechanics/structures/unstructured/Structure.h"
#include "nuto/mechanics/timeIntegration/NewmarkDirect.h"
#include "nuto/math/SparseDirectSolverMUMPS.h"
#include <boost/filesystem.hpp>
#include <iostream>
#include <fstream>
#include <string>

constexpr unsigned int dimension = 2;

class Parameters
{
public:

    static const int mDimension = dimension;

    static const bool mPerformLineSearch = true;
    static const bool mAutomaticTimeStepping = true;

    static constexpr double mMatrixYoungsModulus = 1.0e4;   // concrete
    static constexpr double mMatrixPoissonsRatio = 0.0;
    static constexpr double mMatrixThickness = 10;
    static constexpr double mMatrixNonlocalRadius = 2;
    static constexpr double mMatrixTensileStrength = 3;
    static constexpr double mMatrixCompressiveStrength = 30;
    static constexpr double mMatrixFractureEnergy = 1;

    static constexpr double mTimeStep = 1e-1;
    static constexpr double mMinTimeStep = 1e-3;
    static constexpr double mMaxTimeStep = 1e-1;
    static constexpr double mToleranceForce = 1e-6;
    static constexpr double mSimulationTime = 1.0;
    static constexpr double mLoad = 1.0;

    static const boost::filesystem::path mOutputPath;
    static const boost::filesystem::path mMeshFilePathMatrix;

    static const NuTo::FullVector<double, dimension> mDirectionX;
    static const NuTo::FullVector<double, dimension> mDirectionY;

};

const boost::filesystem::path Parameters::mOutputPath("/home/phuschke/2d_feti/");
const boost::filesystem::path Parameters::mMeshFilePathMatrix("/home/phuschke/develop/nuto/myNutoExamples/MeshFiles/feti00.msh");


const NuTo::FullVector<double, dimension> Parameters::mDirectionX = NuTo::FullVector<double, dimension>::UnitX();
const NuTo::FullVector<double, dimension> Parameters::mDirectionY = NuTo::FullVector<double, dimension>::UnitY();

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void GetSubdomainStiffnessMatrix(const boost::filesystem::path& rMeshPath, NuTo::FullMatrix<double, -1, -1>& rK, NuTo::FullVector<double, -1>& rF)
{

    std::cout << "***********************************" << std::endl;
    std::cout << "**      Structure                **" << std::endl;
    std::cout << "***********************************" << std::endl;

    NuTo::Structure structure00(dimension);
    structure00.SetVerboseLevel(10);
    structure00.SetShowTime(false);


    std::cout << "***********************************" << std::endl;
    std::cout << "**      Integration Scheme       **" << std::endl;
    std::cout << "***********************************" << std::endl;

    NuTo::NewmarkDirect myIntegrationScheme(&structure00);
    myIntegrationScheme.SetTimeStep(0.1);
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
    structure00.InterpolationTypeAdd(interpolationType00, NuTo::Node::DISPLACEMENTS, NuTo::Interpolation::EQUIDISTANT1);


    std::cout << "***********************************" << std::endl;
    std::cout << "**      Matrix                   **" << std::endl;
    std::cout << "***********************************" << std::endl;

    NuTo::FullMatrix<int, Eigen::Dynamic, Eigen::Dynamic> createdGroupIdMatrix = structure00.ImportFromGmsh(rMeshPath.string(), NuTo::ElementData::CONSTITUTIVELAWIP, NuTo::IpData::eIpDataType::STATICDATA);
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
    std::cout << "**      Loads                    **" << std::endl;
    std::cout << "***********************************" << std::endl;

    int groupNodeBCRight = structure00.GroupCreate(NuTo::Groups::eGroupId::Nodes);
    structure00.GroupAddNodeCoordinateRange(groupNodeBCRight, 0, 20.0-1e-6, 20.0+1e-6);

    //int timeDependentConstraint = structure00.ConstraintLinearSetDisplacementNodeGroup(groupNodeBCRight, Parameters::mDirectionX, 1);



    structure00.NodeBuildGlobalDofs();

    structure00.CalculateMaximumIndependentSets();

    NuTo::SparseMatrixCSRVector2General<double> stiffnessMatrixCSRVector2;
    NuTo::FullVector<double,Eigen::Dynamic> dispForceVector, displacementVector, residualVector;

    structure00.BuildGlobalCoefficientMatrix0(stiffnessMatrixCSRVector2, dispForceVector);
    NuTo::SparseMatrixCSRGeneral<double> stiffnessMatrix (stiffnessMatrixCSRVector2);

//
//    // solve
//    NuTo::SparseDirectSolverMUMPS mySolver;
//    stiffnessMatrix.SetOneBasedIndexing();
//    mySolver.Solve(stiffnessMatrix, dispForceVector, displacementVector);

    std::cout << "displacementVec \n" << displacementVector << std::endl;


    rK = stiffnessMatrix;
//    rF = dispForceVector;
    structure00.Info();




}








int main()
{


    boost::filesystem::path meshPath00("/home/phuschke/develop/nuto/myNutoExamples/MeshFiles/feti00.msh");
    boost::filesystem::path meshPath01("/home/phuschke/develop/nuto/myNutoExamples/MeshFiles/feti01.msh");
    try
    {
        NuTo::FullMatrix<double, -1, -1> K00;
        NuTo::FullVector<double, -1> f00(4);
        GetSubdomainStiffnessMatrix(meshPath00, K00, f00);
        NuTo::FullMatrix<double, -1, -1> B00(2,K00.rows());

        B00(0,2) = 1.0;
        B00(1,3) = 1.0;

        f00(0,0) = 1;
        f00(2,0) = 1;

        std::cout << "K00: \n" << K00 << std::endl;
        std::cout << "B00: \n" << B00 << std::endl;
        std::cout << "f00: \n" << f00 << std::endl;

        Eigen::VectorXd disp = K00.colPivHouseholderQr().solve(f00);
            std::cout << "disp \n" << disp << std::endl;

        NuTo::FullMatrix<double, -1, -1> K01;
        NuTo::FullVector<double, -1> f01(4);
        GetSubdomainStiffnessMatrix(meshPath01, K01, f01);
        NuTo::FullMatrix<double, -1, -1> B01(2,K01.rows());

        B01(0,0) = 1.0;
        B01(1,1) = 1.0;

        f01(0,0) = 1;
        f01(2,0) = 1;


//        B01(0,4) = 1.0;
//        B01(1,1) = 1.0;
//        B01(2,5) = 1.0;




        std::cout << "K01: \n" << K01 << std::endl;
        std::cout << "B01: \n" << B01 << std::endl;
        std::cout << "f01: \n" << f01 << std::endl;




        auto invK00 = K00.Inverse();
        auto invK01 = K01.Inverse();



        std::cout << "****************************" << std::endl;
        std::cout << "B00.transpose: \n" << B00.transpose() << std::endl;
        std::cout << "B00 \n" << B00 << std::endl;

        auto tmpMatrix = (B00*invK00*B00.transpose() + B01*invK01*B01.transpose());
        Eigen::MatrixXd lambda = tmpMatrix.colPivHouseholderQr().solve(B01*invK01*f01 - B00*invK00*f00);
        std::cout << "tmpMatrix" << tmpMatrix << std::endl;


                std::cout << "lambda \n" << lambda << std::endl;





        Eigen::MatrixXd u0 = invK00*(f00 + B00.transpose() * lambda);
        Eigen::MatrixXd u1 = invK01*(f01 - B01.transpose() * lambda);

        std::cout << "u0 \n" << u0 << std::endl;
        std::cout << "u1 \n" << u1<< std::endl ;

        std::cout << "B00 u0 \n" << B00 * u0<< std::endl;
        std::cout << "B01 u1 \n" << B01 * u1<< std::endl;

        std::cout << "***********************************" << std::endl;
        std::cout << "**      END                      **" << std::endl;
        std::cout << "***********************************" << std::endl;

    } catch (NuTo::MechanicsException& e)
    {
        std::cout << e.ErrorMessage();

    } catch (NuTo::MathException& e)
    {
        std::cout << e.ErrorMessage();

    }

}


