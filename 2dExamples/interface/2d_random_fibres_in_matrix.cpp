#include "nuto/mechanics/structures/unstructured/Structure.h"
#include "nuto/mechanics/timeIntegration/NewmarkDirect.h"

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

    static constexpr double mMatrixYoungsModulus = 3.0e4; // concrete
    static constexpr double mMatrixPoissonsRatio = 0.2;
    static constexpr double mMatrixThickness = 0.2;
    static constexpr double mMatrixNonlocalRadius = 2;
    static constexpr double mMatrixTensileStrength = 3;
    static constexpr double mMatrixCompressiveStrength = 50;
    static constexpr double mMatrixFractureEnergy = 0.01;

    static constexpr double mFibreYoungsModulus = 2.1e5; // steel
    static constexpr double mFibrePoissonsRatio = 0.2;
    static constexpr double mFibreCrossSection = 0.1;
    static constexpr double mFibreCircumference = 0.6;

    static constexpr double mInterfaceNormalStiffness = 1e6;
    static constexpr double mAlpha = 1;
    static constexpr double mMaxBondStress = 1e1;
    static constexpr double mResidualBondStress = 1e0;
    static constexpr double mSlipAtMaxBondStress = 0.1;
    static constexpr double mSlipAtResidualBondStress = 1;

    static constexpr double mTimeStep = 1e-2;
    static constexpr double mMinTimeStep = 1e-7;
    static constexpr double mMaxTimeStep = 1e-1;
    static constexpr double mToleranceForce = 1e-6;
    static constexpr double mSimulationTime = 1.0;
    static constexpr double mLoad = 0.035;

    static const boost::filesystem::path mOutputPath;
    static const boost::filesystem::path mMeshFilePath;

    static const NuTo::FullVector<double, dimension> mDirectionX;
    static const NuTo::FullVector<double, dimension> mDirectionY;
};

const boost::filesystem::path Parameters::mOutputPath("/home/phuschke/2d_random_fibres_in_matrix/");
const boost::filesystem::path Parameters::mMeshFilePath(
        "/home/phuschke/develop/nuto/myNutoExamples/MeshFiles/interface/2d_horizontal_fibre_double_notched.msh");

const NuTo::FullVector<double, dimension> Parameters::mDirectionX = NuTo::FullVector<double, dimension>::UnitX();
const NuTo::FullVector<double, dimension> Parameters::mDirectionY = NuTo::FullVector<double, dimension>::UnitY();

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

        std::cout << "***********************************" << std::endl;
        std::cout << "**      Section                  **" << std::endl;
        std::cout << "***********************************" << std::endl;

        int matrixSection = myStructure.SectionCreate(NuTo::Section::PLANE_STRESS);
        myStructure.SectionSetThickness(matrixSection, Parameters::mMatrixThickness);

        int fibreSection = myStructure.SectionCreate(NuTo::Section::TRUSS);
        myStructure.SectionSetArea(fibreSection, Parameters::mFibreCrossSection);

        int fibreMatrixBond = myStructure.SectionCreate(NuTo::Section::FIBRE_MATRIX_BOND);
        myStructure.SectionSetCircumference(fibreMatrixBond, Parameters::mFibreCircumference);

        std::cout << "***********************************" << std::endl;
        std::cout << "**      Material                 **" << std::endl;
        std::cout << "***********************************" << std::endl;

        //        int matrixMaterial =
        //        myStructure.ConstitutiveLawCreate(NuTo::Constitutive::eConstitutiveType::LINEAR_ELASTIC_ENGINEERING_STRESS);
        //        myStructure.ConstitutiveLawSetParameterDouble(matrixMaterial,
        //        NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, Parameters::mMatrixYoungsModulus);
        //        myStructure.ConstitutiveLawSetParameterDouble(matrixMaterial,
        //        NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO, Parameters::mMatrixPoissonsRatio);

        int matrixMaterial = myStructure.ConstitutiveLawCreate(
                NuTo::Constitutive::eConstitutiveType::GRADIENT_DAMAGE_ENGINEERING_STRESS);
        myStructure.ConstitutiveLawSetParameterDouble(matrixMaterial,
                                                      NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS,
                                                      Parameters::mMatrixYoungsModulus);
        myStructure.ConstitutiveLawSetParameterDouble(matrixMaterial,
                                                      NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO,
                                                      Parameters::mMatrixPoissonsRatio);
        myStructure.ConstitutiveLawSetParameterDouble(matrixMaterial,
                                                      NuTo::Constitutive::eConstitutiveParameter::TENSILE_STRENGTH,
                                                      Parameters::mMatrixTensileStrength);
        myStructure.ConstitutiveLawSetParameterDouble(matrixMaterial,
                                                      NuTo::Constitutive::eConstitutiveParameter::COMPRESSIVE_STRENGTH,
                                                      Parameters::mMatrixCompressiveStrength);
        myStructure.ConstitutiveLawSetParameterDouble(matrixMaterial,
                                                      NuTo::Constitutive::eConstitutiveParameter::NONLOCAL_RADIUS,
                                                      Parameters::mMatrixNonlocalRadius);
        myStructure.ConstitutiveLawSetParameterDouble(matrixMaterial,
                                                      NuTo::Constitutive::eConstitutiveParameter::FRACTURE_ENERGY,
                                                      Parameters::mMatrixFractureEnergy);

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

        NuTo::FullMatrix<int, Eigen::Dynamic, Eigen::Dynamic> gmshImportedGroupIds =
                myStructure.ImportFromGmsh(Parameters::mMeshFilePath.string(), NuTo::ElementData::CONSTITUTIVELAWIP,
                                           NuTo::IpData::eIpDataType::STATICDATA);

        int groupIdMatrix = gmshImportedGroupIds.GetValue(0, 0);
        int groupIdFibre = gmshImportedGroupIds.GetValue(1, 0);

        int matrixInterpolationType = myStructure.InterpolationTypeCreate(NuTo::Interpolation::eShapeType::TRIANGLE2D);
        myStructure.InterpolationTypeAdd(matrixInterpolationType, NuTo::Node::COORDINATES,
                                         NuTo::Interpolation::EQUIDISTANT2);
        myStructure.InterpolationTypeAdd(matrixInterpolationType, NuTo::Node::DISPLACEMENTS,
                                         NuTo::Interpolation::EQUIDISTANT2);
        myStructure.InterpolationTypeAdd(matrixInterpolationType, NuTo::Node::NONLOCALEQSTRAIN,
                                         NuTo::Interpolation::EQUIDISTANT1); // Gradient damage model

        int fibreInterpolationType = myStructure.InterpolationTypeCreate(NuTo::Interpolation::eShapeType::TRUSSXD);
        myStructure.InterpolationTypeAdd(fibreInterpolationType, NuTo::Node::COORDINATES,
                                         NuTo::Interpolation::EQUIDISTANT2);
        myStructure.InterpolationTypeAdd(fibreInterpolationType, NuTo::Node::DISPLACEMENTS,
                                         NuTo::Interpolation::EQUIDISTANT2);

        int interfaceInterpolationType =
                myStructure.InterpolationTypeCreate(NuTo::Interpolation::eShapeType::INTERFACE);
        myStructure.InterpolationTypeAdd(interfaceInterpolationType, NuTo::Node::COORDINATES,
                                         NuTo::Interpolation::EQUIDISTANT2);
        myStructure.InterpolationTypeAdd(interfaceInterpolationType, NuTo::Node::DISPLACEMENTS,
                                         NuTo::Interpolation::EQUIDISTANT2);
        myStructure.InterpolationTypeSetIntegrationType(interfaceInterpolationType,
                                                        NuTo::IntegrationType::IntegrationType1D2NLobatto4Ip,
                                                        NuTo::IpData::eIpDataType::STATICDATA);

        std::cout << "***********************************" << std::endl;
        std::cout << "**      Matrix                   **" << std::endl;
        std::cout << "***********************************" << std::endl;

        myStructure.ElementGroupSetInterpolationType(groupIdMatrix, matrixInterpolationType);
        myStructure.ElementGroupSetSection(groupIdMatrix, matrixSection);
        myStructure.ElementGroupSetConstitutiveLaw(groupIdMatrix, matrixMaterial);

        std::cout << "***********************************" << std::endl;
        std::cout << "**      Fibre                    **" << std::endl;
        std::cout << "***********************************" << std::endl;

        myStructure.ElementGroupSetInterpolationType(groupIdFibre, fibreInterpolationType);
        myStructure.ElementGroupSetSection(groupIdFibre, fibreSection);
        myStructure.ElementGroupSetConstitutiveLaw(groupIdFibre, fibreMaterial);

        myStructure.ElementTotalConvertToInterpolationType(1.e-6, 10);

        std::cout << "***********************************" << std::endl;
        std::cout << "**      Interface                **" << std::endl;
        std::cout << "***********************************" << std::endl;

        auto elementIds = myStructure.InterfaceElementsCreate(groupIdFibre, interfaceInterpolationType,
                                                              interfaceMaterial, fibreMatrixBond,
                                                              fibreInterpolationType, fibreMaterial, fibreSection);

        // myStructure.Info();

        std::cout << "***********************************" << std::endl;
        std::cout << "**      Boundary Conditions      **" << std::endl;
        std::cout << "***********************************" << std::endl;

        int groupNodeBCLeft = myStructure.GroupCreate(NuTo::Groups::eGroupId::Nodes);
        myStructure.GroupAddNodeCoordinateRange(groupNodeBCLeft, 0, -1e-6, +1e-6);

        myStructure.ConstraintLinearSetDisplacementNodeGroup(groupNodeBCLeft, Parameters::mDirectionX, 0);

        NuTo::FullVector<double, dimension> nodeCoords;

        nodeCoords[0] = 0.0;
        nodeCoords[1] = 20.0;
        int nodeLeft = myStructure.NodeGetIdAtCoordinate(nodeCoords, 1e-6);

        myStructure.ConstraintLinearSetDisplacementNode(nodeLeft, Parameters::mDirectionY, 0);

        std::cout << "***********************************" << std::endl;
        std::cout << "**      Loads                    **" << std::endl;
        std::cout << "***********************************" << std::endl;

        int groupNodeBCRight = myStructure.GroupCreate(NuTo::Groups::eGroupId::Nodes);
        myStructure.GroupAddNodeCoordinateRange(groupNodeBCRight, 0, 160 - 1e-6, 160 + 1e-6);

        int timeDependentConstraint =
                myStructure.ConstraintLinearSetDisplacementNodeGroup(groupNodeBCRight, Parameters::mDirectionX, 1);

        //        int loadCase = 0;
        //        myStructure.LoadCreateNodeForce(0, 1, Parameters::mDirectionX, 1);
        //        myStructure.LoadCreateNodeForce(0, 2, Parameters::mDirectionX, 1);
        //        myStructure.LoadCreateNodeForce(0, 9, Parameters::mDirectionX, 2);

        std::cout << "***********************************" << std::endl;
        std::cout << "**      Visualization            **" << std::endl;
        std::cout << "***********************************" << std::endl;

        myStructure.AddVisualizationComponentDisplacements();
        myStructure.AddVisualizationComponentConstitutive();
        // myStructure.AddVisualizationComponentBondStress();
        // myStructure.ElementGroupExportVtkDataFile(2, Parameters::mOutputPath.string() + "blabla.vtu", true);
        myStructure.AddVisualizationComponentDamage();

        std::cout << "***********************************" << std::endl;
        std::cout << "**      Solver                   **" << std::endl;
        std::cout << "***********************************" << std::endl;


        myStructure.NodeBuildGlobalDofs();


        myStructure.CalculateMaximumIndependentSets();
        // myStructure.Info();
        //        myStructure.CheckCoefficientMatrix_0(1e-6, true);
        //        myStructure.ElementCheckCoefficientMatrix_0(1e-6);


        int groupLoad = myStructure.GroupCreate(NuTo::Groups::eGroupId::Nodes);
        nodeCoords[0] = 160.0;
        nodeCoords[1] = 20.0;
        myStructure.GroupAddNodeRadiusRange(groupLoad, nodeCoords, 0, 1e-6);

        myIntegrationScheme.AddResultGroupNodeForce("myforce", groupLoad);
        auto nodeID = myStructure.GroupGetMemberIds(groupLoad);
        myIntegrationScheme.AddResultNodeDisplacements("mydisplacements", nodeID.at(0, 0));


        auto interfaceIds = myStructure.GroupGetMemberIds(2);
        for (int eleId = 0; eleId < interfaceIds.rows(); ++eleId)
        {
            char myString[200];
            sprintf(myString, "bondStress_%d", eleId);
            myIntegrationScheme.AddResultElementIpData(myString, interfaceIds.at(eleId, 0), NuTo::IpData::BOND_STRESS);

            sprintf(myString, "slip_%d", eleId);
            myIntegrationScheme.AddResultElementIpData(myString, interfaceIds.at(eleId, 0), NuTo::IpData::SLIP);
        }


        NuTo::FullMatrix<double, 2, 2> timeDependentLoad;
        timeDependentLoad(0, 0) = 0;
        timeDependentLoad(1, 0) = Parameters::mSimulationTime;
        timeDependentLoad(0, 1) = 0;
        timeDependentLoad(1, 1) = Parameters::mLoad;

        //        myIntegrationScheme.SetTimeDependentLoadCase(loadCase, timeDependentLoad);
        myIntegrationScheme.SetTimeDependentConstraint(timeDependentConstraint, timeDependentLoad);

        myIntegrationScheme.Solve(Parameters::mSimulationTime);


        std::cout << "\n\n\n Results written to " + Parameters::mOutputPath.string() << std::endl;
    }
    catch (NuTo::MechanicsException& e)
    {
        std::cout << e.ErrorMessage();
    }
    catch (NuTo::MathException& e)
    {
        std::cout << e.ErrorMessage();
    }
    catch (...)
    {
        std::cout << "Something else went wrong" << std::endl;
    }

    std::string command = "paste " + Parameters::mOutputPath.string() + "myforce.dat " +
                          Parameters::mOutputPath.string() + "mydisplacements.dat > " +
                          Parameters::mOutputPath.string() + "forceDisp.dat";
    system(command.c_str());
}
