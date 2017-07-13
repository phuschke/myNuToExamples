#include "nuto/mechanics/structures/unstructured/Structure.h"
#include "nuto/mechanics/timeIntegration/NewmarkDirect.h"

#include <boost/filesystem.hpp>
#include <iostream>
#include <fstream>
#include <string>
#include <time.h>

// Stahlfasern: l = 13 mm, d = 0.2 mm, A = 0.0314 mm^2, U = 0.628 mm, V = 0.4082 mm^3

constexpr unsigned int dimension = 3;

class Parameters
{
public:
    static const int mDimension = dimension;

    static const bool mPerformLineSearch = true;
    static const bool mAutomaticTimeStepping = true;

    static constexpr double mMatrixYoungsModulus = 49083; // concrete
    static constexpr double mMatrixPoissonsRatio = 0.0;
    static constexpr double mMatrixNonlocalRadius = 0.5;
    static constexpr double mMatrixTensileStrength = 6.6;
    static constexpr double mMatrixCompressiveStrength = 115;
    static constexpr double mMatrixFractureEnergy = 0.01;
    static constexpr double mMatrixLengthX = 240.0;
    static constexpr double mMatrixLengthY = 20.0;
    static constexpr double mMatrixLengthZ = 40.0;

    static constexpr double mFibreYoungsModulus = 2.1e5; // steel
    static constexpr double mFibrePoissonsRatio = 0.2;
    static constexpr double mFibreCrossSection = 1.0; // 0.0314;
    static constexpr double mFibreCircumference = 0.628;
    static constexpr double mFibreLength = 13.0;
    static constexpr double mFibreVolume = 0.4082;

    static constexpr double mInterfaceNormalStiffness = 1e6;
    static constexpr double mAlpha = 1;
    static constexpr double mMaxBondStress = 3.e1;
    static constexpr double mResidualBondStress = 1.e1;
    static constexpr double mSlipAtMaxBondStress = 0.1;
    static constexpr double mSlipAtResidualBondStress = 1;

    static constexpr double mTimeStep = 1e-2;
    static constexpr double mMinTimeStep = 1e-4;
    static constexpr double mMaxTimeStep = 1e-1;
    static constexpr double mToleranceForce = 1e-4;
    static constexpr double mSimulationTime = 1.0;
    static constexpr double mLoad = 0.3;

    static const std::string mOutputPath;
    static const std::string mMeshFilePathMatrix;
    static const std::string mMeshFilePathFibre;

    static const NuTo::FullVector<double, dimension> mDirectionX;
    static const NuTo::FullVector<double, dimension> mDirectionY;
    static const NuTo::FullVector<double, dimension> mDirectionZ;
};

const std::string Parameters::mOutputPath("/home/phuschke/3d_gradient_uniaxial_fiber_constraint_dogbone_200/");
const std::string Parameters::mMeshFilePathMatrix(
        "/home/phuschke/develop/nuto/myNutoExamples/MeshFiles/3d/3d_uniaxial_matrix_ehlers.msh");
const std::string
        Parameters::mMeshFilePathFibre("/home/phuschke/develop/nuto/myNutoExamples/MeshFiles/3d/trusses200Dogbone.msh");


const NuTo::FullVector<double, dimension> Parameters::mDirectionX = NuTo::FullVector<double, dimension>::UnitX();
const NuTo::FullVector<double, dimension> Parameters::mDirectionY = NuTo::FullVector<double, dimension>::UnitY();
const NuTo::FullVector<double, dimension> Parameters::mDirectionZ = NuTo::FullVector<double, dimension>::UnitZ();

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
        myStructure.SetNumProcessors(4);

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
        myIntegrationScheme.SetResultDirectory(Parameters::mOutputPath, true);

        std::cout << "***********************************" << std::endl;
        std::cout << "**      Section                  **" << std::endl;
        std::cout << "***********************************" << std::endl;

        int matrixSection = myStructure.SectionCreate(NuTo::Section::VOLUME);

        int fibreSection = myStructure.SectionCreate(NuTo::Section::TRUSS);
        myStructure.SectionSetArea(fibreSection, Parameters::mFibreCrossSection);

        int bondSection = myStructure.SectionCreate(NuTo::Section::FIBRE_MATRIX_BOND);
        myStructure.SectionSetCircumference(bondSection, Parameters::mFibreCircumference);

        std::cout << "***********************************" << std::endl;
        std::cout << "**      Material                 **" << std::endl;
        std::cout << "***********************************" << std::endl;

        //        int matrixMaterial =
        //        myStructure.ConstitutiveLawCreate(NuTo::Constitutive::eConstitutiveType::LINEAR_ELASTIC_ENGINEERING_STRESS);
        //        myStructure.ConstitutiveLawSetParameterDouble(matrixMaterial,
        //        NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, Parameters::mMatrixYoungsModulus);
        //        myStructure.ConstitutiveLawSetParameterDouble(matrixMaterial,
        //        NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO, Parameters::mMatrixPoissonsRatio);

        NuTo::FullVector<double, Eigen::Dynamic> myDamageLaw(1);
        myDamageLaw(0) = NuTo::Constitutive::eDamageLawType::ISOTROPIC_EXPONENTIAL_SOFTENING;

        int matrixMaterial = myStructure.ConstitutiveLawCreate(
                NuTo::Constitutive::eConstitutiveType::GRADIENT_DAMAGE_ENGINEERING_STRESS);
        myStructure.ConstitutiveLawSetParameterDouble(matrixMaterial,
                                                      NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS,
                                                      Parameters::mMatrixYoungsModulus);
        myStructure.ConstitutiveLawSetParameterDouble(matrixMaterial,
                                                      NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO,
                                                      Parameters::mMatrixPoissonsRatio);
        myStructure.ConstitutiveLawSetParameterDouble(matrixMaterial,
                                                      NuTo::Constitutive::eConstitutiveParameter::NONLOCAL_RADIUS,
                                                      Parameters::mMatrixNonlocalRadius);
        myStructure.ConstitutiveLawSetParameterDouble(matrixMaterial,
                                                      NuTo::Constitutive::eConstitutiveParameter::TENSILE_STRENGTH,
                                                      Parameters::mMatrixTensileStrength);
        myStructure.ConstitutiveLawSetParameterDouble(matrixMaterial,
                                                      NuTo::Constitutive::eConstitutiveParameter::COMPRESSIVE_STRENGTH,
                                                      Parameters::mMatrixCompressiveStrength);
        myStructure.ConstitutiveLawSetParameterDouble(matrixMaterial,
                                                      NuTo::Constitutive::eConstitutiveParameter::FRACTURE_ENERGY,
                                                      Parameters::mMatrixFractureEnergy);
        myStructure.ConstitutiveLawSetParameterFullVectorDouble(
                matrixMaterial, NuTo::Constitutive::eConstitutiveParameter::DAMAGE_LAW, myDamageLaw);
        //


        int fibreMaterial = myStructure.ConstitutiveLawCreate(
                NuTo::Constitutive::eConstitutiveType::LINEAR_ELASTIC_ENGINEERING_STRESS);
        myStructure.ConstitutiveLawSetParameterDouble(fibreMaterial,
                                                      NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS,
                                                      Parameters::mFibreYoungsModulus);
        myStructure.ConstitutiveLawSetParameterDouble(fibreMaterial,
                                                      NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO,
                                                      Parameters::mFibrePoissonsRatio);

        int interfaceMaterial =
                myStructure.ConstitutiveLawCreate(NuTo::Constitutive::eConstitutiveType::FIBRE_MATRIX_BOND_STRESS_SLIP);
        myStructure.ConstitutiveLawSetParameterDouble(interfaceMaterial,
                                                      NuTo::Constitutive::eConstitutiveParameter::NORMAL_STIFFNESS,
                                                      Parameters::mInterfaceNormalStiffness);
        myStructure.ConstitutiveLawSetParameterDouble(
                interfaceMaterial, NuTo::Constitutive::eConstitutiveParameter::ALPHA, Parameters::mAlpha);
        myStructure.ConstitutiveLawSetParameterDouble(interfaceMaterial,
                                                      NuTo::Constitutive::eConstitutiveParameter::MAX_BOND_STRESS,
                                                      Parameters::mMaxBondStress);
        myStructure.ConstitutiveLawSetParameterDouble(interfaceMaterial,
                                                      NuTo::Constitutive::eConstitutiveParameter::RESIDUAL_BOND_STRESS,
                                                      Parameters::mResidualBondStress);
        myStructure.ConstitutiveLawSetParameterDouble(
                interfaceMaterial, NuTo::Constitutive::eConstitutiveParameter::SLIP_AT_MAX_BOND_STRESS,
                Parameters::mSlipAtMaxBondStress);
        myStructure.ConstitutiveLawSetParameterDouble(
                interfaceMaterial, NuTo::Constitutive::eConstitutiveParameter::SLIP_AT_RESIDUAL_BOND_STRESS,
                Parameters::mSlipAtResidualBondStress);

        std::cout << "***********************************" << std::endl;
        std::cout << "**      Interpolation Type       **" << std::endl;
        std::cout << "***********************************" << std::endl;

        int matrixInterpolationType =
                myStructure.InterpolationTypeCreate(NuTo::Interpolation::eShapeType::TETRAHEDRON3D);
        myStructure.InterpolationTypeAdd(matrixInterpolationType, NuTo::Node::COORDINATES,
                                         NuTo::Interpolation::EQUIDISTANT1);
        myStructure.InterpolationTypeAdd(matrixInterpolationType, NuTo::Node::DISPLACEMENTS,
                                         NuTo::Interpolation::EQUIDISTANT2);
        myStructure.InterpolationTypeAdd(matrixInterpolationType, NuTo::Node::NONLOCALEQSTRAIN,
                                         NuTo::Interpolation::EQUIDISTANT1);

        int fibreInterpolationType = myStructure.InterpolationTypeCreate(NuTo::Interpolation::eShapeType::TRUSSXD);
        myStructure.InterpolationTypeAdd(fibreInterpolationType, NuTo::Node::COORDINATES,
                                         NuTo::Interpolation::EQUIDISTANT1);
        myStructure.InterpolationTypeAdd(fibreInterpolationType, NuTo::Node::DISPLACEMENTS,
                                         NuTo::Interpolation::EQUIDISTANT1);

        int interfaceInterpolationType =
                myStructure.InterpolationTypeCreate(NuTo::Interpolation::eShapeType::INTERFACE);
        myStructure.InterpolationTypeAdd(interfaceInterpolationType, NuTo::Node::COORDINATES,
                                         NuTo::Interpolation::EQUIDISTANT1);
        myStructure.InterpolationTypeAdd(interfaceInterpolationType, NuTo::Node::DISPLACEMENTS,
                                         NuTo::Interpolation::EQUIDISTANT1);

        std::cout << "***********************************" << std::endl;
        std::cout << "**      Matrix                   **" << std::endl;
        std::cout << "***********************************" << std::endl;

        NuTo::FullMatrix<int, Eigen::Dynamic, Eigen::Dynamic> createdGroupIdMatrix =
                myStructure.ImportFromGmsh(Parameters::mMeshFilePathMatrix, NuTo::ElementData::CONSTITUTIVELAWIP,
                                           NuTo::IpData::eIpDataType::STATICDATA);
        int groupIdMatrix = createdGroupIdMatrix.GetValue(0, 0);

        myStructure.ElementGroupSetSection(groupIdMatrix, matrixSection);
        myStructure.ElementGroupSetConstitutiveLaw(groupIdMatrix, matrixMaterial);
        myStructure.ElementTotalConvertToInterpolationType(1.e-6, 10);
        myStructure.ElementGroupSetInterpolationType(groupIdMatrix, matrixInterpolationType);
        myStructure.InterpolationTypeSetIntegrationType(matrixInterpolationType,
                                                        NuTo::IntegrationType::IntegrationType3D4NGauss4Ip,
                                                        NuTo::IpData::eIpDataType::STATICDATA);
        myStructure.ElementTotalConvertToInterpolationType(1.e-6, 10);

        std::cout << "***********************************" << std::endl;
        std::cout << "**      Boundary Conditions      **" << std::endl;
        std::cout << "***********************************" << std::endl;

        NuTo::FullVector<double, Parameters::mDimension> nodeCoords;

        int groupNodeBCLeft = myStructure.GroupCreate(NuTo::Groups::eGroupId::Nodes);
        myStructure.GroupAddNodeCoordinateRange(groupNodeBCLeft, 0, 00.0 - 1e-6, 00.0 + 1e-6);

        myStructure.ConstraintLinearSetDisplacementNodeGroup(groupNodeBCLeft, Parameters::mDirectionX, 0);
        myStructure.ConstraintLinearSetDisplacementNodeGroup(groupNodeBCLeft, Parameters::mDirectionY, 0);
        myStructure.ConstraintLinearSetDisplacementNodeGroup(groupNodeBCLeft, Parameters::mDirectionZ, 0);

        int groupNodeBCRight = myStructure.GroupCreate(NuTo::Groups::eGroupId::Nodes);
        myStructure.GroupAddNodeCoordinateRange(groupNodeBCRight, 0, Parameters::mMatrixLengthX - 1e-6,
                                                Parameters::mMatrixLengthX + 1e-6);

        myStructure.ConstraintLinearSetDisplacementNodeGroup(groupNodeBCRight, Parameters::mDirectionY, 0);
        myStructure.ConstraintLinearSetDisplacementNodeGroup(groupNodeBCRight, Parameters::mDirectionZ, 0);


        //        int groupNodeBCBottom = myStructure.GroupCreate(NuTo::Groups::eGroupId::Nodes);
        //        nodeCoords[0] = 0.0;
        //        nodeCoords[1] = 0;
        //        nodeCoords[2] = 0;
        //        myStructure.GroupAddNodeRadiusRange(groupNodeBCBottom, nodeCoords, 0, 1e-2);
        //
        //        nodeCoords[0] = 0.0;
        //        nodeCoords[1] = 0;
        //        nodeCoords[2] = 40;
        //        myStructure.GroupAddNodeRadiusRange(groupNodeBCBottom, nodeCoords, 0, 1e-2);
        //
        ////        myStructure.GroupAddNodeCoordinateRange(groupNodeBCBottom, 1, 0.0 - 1e-6, 0.0 + 1e-6);
        //
        //        myStructure.ConstraintLinearSetDisplacementNodeGroup(groupNodeBCBottom, Parameters::mDirectionY, 0);
        //
        //
        //        int groupNodeBCFront = myStructure.GroupCreate(NuTo::Groups::eGroupId::Nodes);
        //
        //        nodeCoords[0] = 0.0;
        //        nodeCoords[1] = -20;
        //        nodeCoords[2] = 0;
        //        myStructure.GroupAddNodeRadiusRange(groupNodeBCFront, nodeCoords, 0, 1e-2);
        //
        //        nodeCoords[0] = 0.0;
        //        nodeCoords[1] = +20;
        //        nodeCoords[2] = 0;
        //        myStructure.GroupAddNodeRadiusRange(groupNodeBCFront, nodeCoords, 0, 1e-2);

        //        myStructure.GroupAddNodeCoordinateRange(groupNodeBCFront, 2, 0.0 - 1e-6, 0.0 + 1e-6);

        //        myStructure.ConstraintLinearSetDisplacementNodeGroup(groupNodeBCFront, Parameters::mDirectionZ, 0);


        std::cout << "***********************************" << std::endl;
        std::cout << "**      Fibre                    **" << std::endl;
        std::cout << "***********************************" << std::endl;

        NuTo::FullMatrix<int, Eigen::Dynamic, Eigen::Dynamic> createdGroupIdFibre =
                myStructure.ImportFromGmsh(Parameters::mMeshFilePathFibre, NuTo::ElementData::CONSTITUTIVELAWIP,
                                           NuTo::IpData::eIpDataType::STATICDATA);
        int groupIdFiber = createdGroupIdFibre.GetValue(0, 0);

        myStructure.ElementGroupSetSection(groupIdFiber, fibreSection);
        myStructure.ElementGroupSetConstitutiveLaw(groupIdFiber, fibreMaterial);
        myStructure.ElementGroupSetInterpolationType(groupIdFiber, fibreInterpolationType);
        myStructure.ElementConvertToInterpolationType(groupIdFiber);

        std::cout << "***********************************" << std::endl;
        std::cout << "**      Constraints              **" << std::endl;
        std::cout << "***********************************" << std::endl;

        int numGroups = 2;

        int groupMatrixNodes = myStructure.GroupCreate(NuTo::Groups::Nodes);
        myStructure.GroupAddNodeFromElementGroupCoordinateRange(groupMatrixNodes, groupIdMatrix, 0, 0, 10);

        int groupMatrixElements = myStructure.GroupCreate(NuTo::Groups::Elements);
        myStructure.GroupAddElementsFromNodes(groupMatrixElements, groupMatrixNodes, false);

        int groupConstraintNodes = myStructure.GroupCreate(NuTo::Groups::Nodes);
        myStructure.GroupAddNodeFromElementGroupCoordinateRange(groupConstraintNodes, groupIdFiber, 0, 0, 10);

        myStructure.GroupAddNodesFromElements(groupConstraintNodes, groupIdFiber);


        int numNearestNeighbours = 1;

        auto nodeIds = myStructure.GroupGetMemberIds(groupConstraintNodes);
        std::cout << "nodeIds.rows()" << nodeIds.rows() << std::endl;
        for (int iNode = 0; iNode < nodeIds.rows(); ++iNode)
        {
            myStructure.ConstraintLinearEquationNodeToElementCreate(
                    nodeIds.at(iNode, 0), groupIdMatrix, NuTo::Node::eAttributes::DISPLACEMENTS, numNearestNeighbours);
        }


        //        for (int i = 2; i < 3; ++i)
        //        {
        //
        //            constexpr int numNearestNeighbours = 1;
        //            int numSearchDomains = 5;
        //            constexpr double length = 110.0;
        //            double deltaLength = length / numSearchDomains;
        //
        //            clock_t t;
        //            t = clock();
        //
        //            for (int iDomain = 0; iDomain < numSearchDomains; ++iDomain)
        //            {
        //                int groupMatrixNodes = myStructure.GroupCreate(NuTo::Groups::Nodes);
        //                myStructure.GroupAddNodeFromElementGroupCoordinateRange(groupMatrixNodes, groupIdMatrix, 0,
        //                iDomain * deltaLength, (iDomain + 1) * deltaLength);
        //
        //                int groupMatrixElements = myStructure.GroupCreate(NuTo::Groups::Elements);
        //                myStructure.GroupAddElementsFromNodes(groupMatrixElements, groupMatrixNodes, false);
        //
        //                int groupConstraintNodes = myStructure.GroupCreate(NuTo::Groups::Nodes);
        //                myStructure.GroupAddNodeFromElementGroupCoordinateRange(groupConstraintNodes, groupIdFiber, 0,
        //                iDomain * deltaLength, (iDomain + 1) * deltaLength);
        //
        //                auto nodeIds = myStructure.GroupGetMemberIds(groupConstraintNodes);
        //
        //                for (int iNode = 0; iNode < nodeIds.rows(); ++iNode)
        //                {
        //                    myStructure.ConstraintLinearEquationNodeToElementCreate(nodeIds.at(iNode, 0),
        //                    groupMatrixElements, NuTo::Node::eAttributes::DISPLACEMENTS, numNearestNeighbours);
        //                }
        //            }
        //
        //            t = clock() - t;
        //            std::ofstream myfile;
        //            myfile.open("time.txt", std::ofstream::app);
        //            myfile << "numSearchDomain: \t" << numSearchDomains << std::endl;
        //            myfile << "time: \t" << ((float) t) / CLOCKS_PER_SEC << std::endl;
        //            myfile.close();
        //        }
        //
        //        return EXIT_SUCCESS;

        std::cout << "***********************************" << std::endl;
        std::cout << "**      Interface                **" << std::endl;
        std::cout << "***********************************" << std::endl;

        auto bondAndFibreElements =
                myStructure.InterfaceElementsCreate(groupIdFiber, interfaceInterpolationType, interfaceMaterial,
                                                    bondSection, fibreInterpolationType, fibreMaterial, fibreSection);

        int groupNewFibreElements = myStructure.GroupCreate(NuTo::Groups::Elements);
        for (int i = 0; i < bondAndFibreElements.rows(); ++i)
        {
            myStructure.GroupAddElement(groupNewFibreElements, bondAndFibreElements.at(i, 1));
        }


        std::cout << "***********************************" << std::endl;
        std::cout << "**      Loads                    **" << std::endl;
        std::cout << "***********************************" << std::endl;


        int timeDependentConstraint =
                myStructure.ConstraintLinearSetDisplacementNodeGroup(groupNodeBCRight, Parameters::mDirectionX, 1);


        std::cout << "***********************************" << std::endl;
        std::cout << "**      Visualization            **" << std::endl;
        std::cout << "***********************************" << std::endl;

        myStructure.AddVisualizationComponent(groupIdMatrix, NuTo::VisualizeBase::DISPLACEMENTS);
        //        myStructure.AddVisualizationComponent(groupIdMatrix, NuTo::VisualizeBase::ENGINEERING_STRAIN);
        //        myStructure.AddVisualizationComponent(groupIdMatrix, NuTo::VisualizeBase::ENGINEERING_STRESS);
        myStructure.AddVisualizationComponent(groupIdMatrix, NuTo::VisualizeBase::DAMAGE);
        myStructure.AddVisualizationComponent(groupIdMatrix, NuTo::VisualizeBase::LOCAL_EQ_STRAIN);

        myStructure.AddVisualizationComponent(groupNewFibreElements, NuTo::VisualizeBase::DISPLACEMENTS);


        std::cout << "***********************************" << std::endl;
        std::cout << "**      Solver                   **" << std::endl;
        std::cout << "***********************************" << std::endl;

        myStructure.NodeBuildGlobalDofs();


        myStructure.CalculateMaximumIndependentSets();


        int groupLoad = myStructure.GroupCreate(NuTo::Groups::eGroupId::Nodes);

        nodeCoords[0] = 240.0;
        nodeCoords[1] = 0;
        nodeCoords[2] = 0;
        myStructure.GroupAddNodeRadiusRange(groupLoad, nodeCoords, 0, 1e-2);

        myIntegrationScheme.AddResultGroupNodeForce("myforce", groupNodeBCRight);
        auto nodeID = myStructure.GroupGetMemberIds(groupLoad);
        myIntegrationScheme.AddResultNodeDisplacements("mydisplacements", nodeID.at(0, 0));


        NuTo::FullMatrix<double, 2, 2> timeDependentLoad;
        timeDependentLoad(0, 0) = 0;
        timeDependentLoad(1, 0) = 1.0 * Parameters::mSimulationTime;

        timeDependentLoad(0, 1) = 0;
        timeDependentLoad(1, 1) = 1.0 * Parameters::mLoad;

        myIntegrationScheme.SetTimeDependentConstraint(timeDependentConstraint, timeDependentLoad);

        myIntegrationScheme.Solve(Parameters::mSimulationTime);

        std::string command = "paste " + Parameters::mOutputPath + "myforce.dat " + Parameters::mOutputPath +
                              "mydisplacements.dat > " + Parameters::mOutputPath + "forceDisp.dat";
        system(command.c_str());

        std::cout << "***********************************" << std::endl;
        std::cout << "**      END                      **" << std::endl;
        std::cout << "***********************************" << std::endl;
    }
    catch (NuTo::MechanicsException& e)
    {
        std::cout << e.ErrorMessage();
    }
    catch (NuTo::MathException& e)
    {
        std::cout << e.ErrorMessage();
    }
}
