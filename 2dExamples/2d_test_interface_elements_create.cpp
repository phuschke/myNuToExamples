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

    static const bool mPerformLineSearch = false;
    static const bool mAutomaticTimeStepping = false;

    static constexpr double mMatrixYoungsModulus = 4.0e4; // concrete
    static constexpr double mMatrixPoissonsRatio = 0.2;
    static constexpr double mMatrixThickness = 10;

    static constexpr double mFibreYoungsModulus = 2.1e11; // steel
    static constexpr double mFibrePoissonsRatio = 0.2;
    static constexpr double mFibreCrossSection = 0.01;

    static constexpr double mInterfaceNormalStiffness = 1e8;
    static constexpr double mAlpha = 1;
    static constexpr double mMaxBondStress = 2e5;
    static constexpr double mResidualBondStress = 1e5;
    static constexpr double mSlipAtMaxBondStress = 0.2;
    static constexpr double mSlipAtResidualBondStress = 0.5;

    static constexpr double mTimeStep = 1e-2;
    static constexpr double mMinTimeStep = 1e-5;
    static constexpr double mMaxTimeStep = 1e-1;
    static constexpr double mToleranceForce = 1e-6;
    static constexpr double mSimulationTime = 1.0;
    static constexpr double mLoad = 6;

    static const boost::filesystem::path mOutputPath;
    static const boost::filesystem::path mMeshFilePath;


    static const NuTo::FullVector<double, dimension> mDirectionX;
    static const NuTo::FullVector<double, dimension> mDirectionY;
};

const boost::filesystem::path Parameters::mOutputPath("/home/phuschke/2d_test_interface_elements_create/");
const boost::filesystem::path
        Parameters::mMeshFilePath("/home/phuschke/develop/nuto/myNutoExamples/MeshFiles/2d_interface_generation.msh");

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

        std::cout << "***********************************" << std::endl;
        std::cout << "**      Material                 **" << std::endl;
        std::cout << "***********************************" << std::endl;

        int matrixMaterial = myStructure.ConstitutiveLawCreate(
                NuTo::Constitutive::eConstitutiveType::LINEAR_ELASTIC_ENGINEERING_STRESS);
        myStructure.ConstitutiveLawSetParameterDouble(matrixMaterial,
                                                      NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS,
                                                      Parameters::mMatrixYoungsModulus);
        myStructure.ConstitutiveLawSetParameterDouble(matrixMaterial,
                                                      NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO,
                                                      Parameters::mMatrixPoissonsRatio);

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

        int matrixInterpolationType = myStructure.InterpolationTypeCreate(NuTo::Interpolation::eShapeType::QUAD2D);
        myStructure.InterpolationTypeAdd(matrixInterpolationType, NuTo::Node::COORDINATES,
                                         NuTo::Interpolation::EQUIDISTANT1);
        myStructure.InterpolationTypeAdd(matrixInterpolationType, NuTo::Node::DISPLACEMENTS,
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

        myStructure.ElementGroupSetInterpolationType(groupIdMatrix, matrixInterpolationType);
        myStructure.ElementGroupSetSection(groupIdMatrix, matrixSection);
        myStructure.ElementGroupSetConstitutiveLaw(groupIdMatrix, matrixMaterial);

        std::cout << "***********************************" << std::endl;
        std::cout << "**      Fibre                    **" << std::endl;
        std::cout << "***********************************" << std::endl;

        myStructure.ElementGroupSetInterpolationType(groupIdFibre, fibreInterpolationType);
        myStructure.ElementGroupSetSection(groupIdFibre, fibreSection);
        myStructure.ElementGroupSetConstitutiveLaw(groupIdFibre, fibreMaterial);

        myStructure.Info();
        myStructure.ElementTotalConvertToInterpolationType(1.e-6, 10);
        myStructure.NodeBuildGlobalDofs();
        myStructure.Info();
        myStructure.CheckCoefficientMatrix_0(1e-6, true);

        std::cout << "***********************************" << std::endl;
        std::cout << "**      Interface                **" << std::endl;
        std::cout << "***********************************" << std::endl;

        NuTo::FullVector<int, Eigen::Dynamic> elementIdsFibre;
        elementIdsFibre = myStructure.GroupGetMemberIds(groupIdFibre);

        int fibreElementGroup = myStructure.GroupCreate(NuTo::Groups::eGroupId::Elements);
        myStructure.GroupAddElement(fibreElementGroup, elementIdsFibre.at(0, 0));

        myStructure.InterfaceElementsCreate(fibreElementGroup, interfaceInterpolationType, interfaceMaterial,
                                            fibreInterpolationType, fibreMaterial, fibreSection);

        myStructure.Info();

        std::cout << "***********************************" << std::endl;
        std::cout << "**      Boundary Conditions      **" << std::endl;
        std::cout << "***********************************" << std::endl;

        int groupNodeBCLeft = myStructure.GroupCreate(NuTo::Groups::eGroupId::Nodes);
        myStructure.GroupAddNodeCoordinateRange(groupNodeBCLeft, 0, -1e-6, +1e-6);

        myStructure.ConstraintLinearSetDisplacementNodeGroup(groupNodeBCLeft, Parameters::mDirectionX, 0);


        NuTo::FullVector<double, dimension> nodeCoords;

        nodeCoords[0] = 0.0;
        nodeCoords[1] = 1;
        int nodeLeft = myStructure.NodeGetIdAtCoordinate(nodeCoords, 1e-6);


        myStructure.ConstraintLinearSetDisplacementNode(nodeLeft, Parameters::mDirectionY, 0);


        std::cout << "***********************************" << std::endl;
        std::cout << "**      Loads                **" << std::endl;
        std::cout << "***********************************" << std::endl;

        int groupNodeBCRight = myStructure.GroupCreate(NuTo::Groups::eGroupId::Nodes);
        myStructure.GroupAddNodeCoordinateRange(groupNodeBCRight, 0, 9 - 1e-6, 9 + 1e-6);

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

        std::cout << "***********************************" << std::endl;
        std::cout << "**      Solver                   **" << std::endl;
        std::cout << "***********************************" << std::endl;

        myStructure.NodeBuildGlobalDofs();
        myStructure.Info();
        myStructure.CheckCoefficientMatrix_0(1e-6, true);

        myStructure.ElementCheckCoefficientMatrix_0(1e-6);

        myStructure.NodeBuildGlobalDofs();
        myStructure.CalculateMaximumIndependentSets();


        int groupLoad = myStructure.GroupCreate(NuTo::Groups::eGroupId::Nodes);
        myStructure.GroupAddNode(groupLoad, 11);


        myIntegrationScheme.AddResultGroupNodeForce("myforce", groupLoad);
        myIntegrationScheme.AddResultNodeDisplacements("mydisplacements", 11);


        NuTo::FullMatrix<double, 2, 2> timeDependentLoad;
        timeDependentLoad(0, 0) = 0;
        timeDependentLoad(1, 0) = Parameters::mSimulationTime;
        timeDependentLoad(0, 1) = 0;
        timeDependentLoad(1, 1) = Parameters::mLoad;

        //        myIntegrationScheme.SetTimeDependentLoadCase(loadCase, timeDependentLoad);
        myIntegrationScheme.SetTimeDependentConstraint(timeDependentConstraint, timeDependentLoad);

        myIntegrationScheme.Solve(Parameters::mSimulationTime);


        std::string command = "paste " + Parameters::mOutputPath.string() + "myforce.dat " +
                              Parameters::mOutputPath.string() + "mydisplacements.dat > " +
                              Parameters::mOutputPath.string() + "forceDisp.dat";
        system(command.c_str());


        std::cout << "\n\n\n Results written to " + Parameters::mOutputPath.string() << std::endl;
    }
    catch (NuTo::MechanicsException& e)
    {
        std::cout << e.ErrorMessage();
        return -1;
    }
}
