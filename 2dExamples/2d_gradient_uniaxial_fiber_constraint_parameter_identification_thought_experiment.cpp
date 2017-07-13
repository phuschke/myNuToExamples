#include "nuto/mechanics/structures/unstructured/Structure.h"
#include "nuto/mechanics/timeIntegration/NewmarkDirect.h"

#include <boost/filesystem.hpp>
#include <iostream>
#include <fstream>
#include <string>
#include <time.h>
#include "nuto/mechanics/constitutive/laws/GradientDamageEngineeringStress.h"
#include "nuto/mechanics/constitutive/staticData/ConstitutiveStaticDataGradientDamage.h"
// Stahlfasern: l = 13 mm, d = 0.2 mm, A = 0.0314 mm^2, U = 0.628 mm, V = 0.4082 mm^3


int main(int argc, char* argv[])
{

    if (not argc == 5)
    {
        std::cout << "Invalid number of input arguments: matrixFractureEnergy, fibreCrossSection, maxBondStress, "
                     "slipAtMaxBondStress"
                  << std::endl;
        return 1;
    }

    // Parameters
    const int dimension = 2;
    const bool performLineSearch = true;
    const bool automaticTimeStepping = true;

    const double tol = 1.0e-4;

    const double matrixYoungsModulus = 40083; // concrete
    const double matrixTensileStrength = 6;
    const double matrixCompressiveStrength = 60;
    const double matrixPoissonsRatio = 0.0;

    const double matrixNonlocalRadius = std::stod(argv[2]);
    const double matrixFractureEnergy = std::stod(argv[1]); // 0.01;
    const double matrixLengthX = 40.0;
    // const double matrixLengthY = 20.0;
    const double matrixLengthZ = 1.0;

    const double fibreYoungsModulus = 2.1e5; // steel
    const double fibrePoissonsRatio = 0.2;
    const double fibreCrossSection = 0.0314;
    const double fibreCircumference = 0.628;
    const double fibreLength = 13.0;
    const double fibreVolume = 0.4082;

    const double mInterfaceNormalStiffness = 1e6;
    const double mAlpha = 1;
    const double maxBondStress = std::stod(argv[3]); // 3.0
    const double mResidualBondStress = 0.7 * maxBondStress;
    const double mSlipAtMaxBondStress = std::stod(argv[4]); // 0.01
    const double mSlipAtResidualBondStress = 10 * mSlipAtMaxBondStress;

    const double mTimeStep = 0.005;
    const double minTimeStep = 1e-5;
    const double maxTimeStep = 1;
    const double mToleranceDisp = 1e-3;
    const double mToleranceNlEqStrain = 1e-4;
    const double mSimulationTime = 1.0;
    const double mLoad = 1;

    const NuTo::FullVector<double, 2> directionX = NuTo::FullVector<double, 2>::UnitX();
    const NuTo::FullVector<double, 2> directionY = NuTo::FullVector<double, 2>::UnitY();

    char tmpChar[200];
    sprintf(tmpChar, "/home/phuschke/2d_parameter_identification/"
                     "2d_gradient_uniaxial_thought_experiment_fractureEnergy_%f_nonlocal_radius_%f_maxBondStress_%f_"
                     "mSlipAtMaxBondStress_%f/",
            matrixFractureEnergy, matrixNonlocalRadius, maxBondStress, mSlipAtMaxBondStress);
    const std::string resultDir = tmpChar;

    const std::string meshFileMatrix("/home/phuschke/meshFiles/2d/2d_uniaxial_matrix_markovic.msh");
    const std::string meshFileFibre("/home/phuschke/meshFiles/2d/trusses_markovic.msh");


    try
    {
        std::cout << "***********************************" << std::endl;
        std::cout << "**      Structure                **" << std::endl;
        std::cout << "***********************************" << std::endl;

        NuTo::Structure myStructure(dimension);
        myStructure.SetVerboseLevel(10);
        myStructure.SetShowTime(false);
        myStructure.SetNumProcessors(4);

        std::cout << "***********************************" << std::endl;
        std::cout << "**      Integration Scheme       **" << std::endl;
        std::cout << "***********************************" << std::endl;

        NuTo::NewmarkDirect myIntegrationScheme(&myStructure);
        myIntegrationScheme.SetTimeStep(mTimeStep);
        myIntegrationScheme.SetMinTimeStep(minTimeStep);
        myIntegrationScheme.SetMaxTimeStep(maxTimeStep);

        myIntegrationScheme.SetToleranceResidual(NuTo::Node::DISPLACEMENTS, mToleranceDisp);
        myIntegrationScheme.SetToleranceResidual(NuTo::Node::NONLOCALEQSTRAIN, mToleranceNlEqStrain);
        myIntegrationScheme.SetAutomaticTimeStepping(automaticTimeStepping);
        myIntegrationScheme.SetPerformLineSearch(performLineSearch);
        myIntegrationScheme.SetResultDirectory(resultDir, true);

        std::cout << "***********************************" << std::endl;
        std::cout << "**      Section                  **" << std::endl;
        std::cout << "***********************************" << std::endl;

        int matrixSection = myStructure.SectionCreate(NuTo::Section::PLANE_STRAIN);
        myStructure.SectionSetThickness(matrixSection, matrixLengthZ);

        int matrixSection_thin = myStructure.SectionCreate(NuTo::Section::PLANE_STRAIN);
        myStructure.SectionSetThickness(matrixSection, matrixLengthZ * 0.9);

        int fibreSection = myStructure.SectionCreate(NuTo::Section::TRUSS);
        myStructure.SectionSetArea(fibreSection, fibreCrossSection);

        int bondSection = myStructure.SectionCreate(NuTo::Section::FIBRE_MATRIX_BOND);
        myStructure.SectionSetCircumference(bondSection, fibreCircumference);

        std::cout << "***********************************" << std::endl;
        std::cout << "**      Material                 **" << std::endl;
        std::cout << "***********************************" << std::endl;

        int matrixMaterial_elastic = myStructure.ConstitutiveLawCreate(
                NuTo::Constitutive::eConstitutiveType::LINEAR_ELASTIC_ENGINEERING_STRESS);
        myStructure.ConstitutiveLawSetParameterDouble(matrixMaterial_elastic,
                                                      NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS,
                                                      matrixYoungsModulus);
        myStructure.ConstitutiveLawSetParameterDouble(matrixMaterial_elastic,
                                                      NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO,
                                                      matrixPoissonsRatio);

        int matrixMaterial = myStructure.ConstitutiveLawCreate(
                NuTo::Constitutive::eConstitutiveType::GRADIENT_DAMAGE_ENGINEERING_STRESS);
        myStructure.ConstitutiveLawSetParameterDouble(
                matrixMaterial, NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, matrixYoungsModulus);
        myStructure.ConstitutiveLawSetParameterDouble(
                matrixMaterial, NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO, matrixPoissonsRatio);
        myStructure.ConstitutiveLawSetParameterDouble(
                matrixMaterial, NuTo::Constitutive::eConstitutiveParameter::NONLOCAL_RADIUS, matrixNonlocalRadius);
        myStructure.ConstitutiveLawSetParameterDouble(
                matrixMaterial, NuTo::Constitutive::eConstitutiveParameter::TENSILE_STRENGTH, matrixTensileStrength);
        myStructure.ConstitutiveLawSetParameterDouble(matrixMaterial,
                                                      NuTo::Constitutive::eConstitutiveParameter::COMPRESSIVE_STRENGTH,
                                                      matrixCompressiveStrength);
        myStructure.ConstitutiveLawSetParameterDouble(
                matrixMaterial, NuTo::Constitutive::eConstitutiveParameter::FRACTURE_ENERGY, matrixFractureEnergy);
        myStructure.ConstitutiveLawSetParameterDouble(
                matrixMaterial, NuTo::Constitutive::eConstitutiveParameter::NONLOCAL_RADIUS_PARAMETER, 2);

        int fibreMaterial = myStructure.ConstitutiveLawCreate(
                NuTo::Constitutive::eConstitutiveType::LINEAR_ELASTIC_ENGINEERING_STRESS);
        myStructure.ConstitutiveLawSetParameterDouble(
                fibreMaterial, NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, fibreYoungsModulus);
        myStructure.ConstitutiveLawSetParameterDouble(
                fibreMaterial, NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO, fibrePoissonsRatio);

        int interfaceMaterial =
                myStructure.ConstitutiveLawCreate(NuTo::Constitutive::eConstitutiveType::FIBRE_MATRIX_BOND_STRESS_SLIP);
        myStructure.ConstitutiveLawSetParameterDouble(interfaceMaterial,
                                                      NuTo::Constitutive::eConstitutiveParameter::NORMAL_STIFFNESS,
                                                      mInterfaceNormalStiffness);
        myStructure.ConstitutiveLawSetParameterDouble(interfaceMaterial,
                                                      NuTo::Constitutive::eConstitutiveParameter::ALPHA, mAlpha);
        myStructure.ConstitutiveLawSetParameterDouble(
                interfaceMaterial, NuTo::Constitutive::eConstitutiveParameter::MAX_BOND_STRESS, maxBondStress);
        myStructure.ConstitutiveLawSetParameterDouble(interfaceMaterial,
                                                      NuTo::Constitutive::eConstitutiveParameter::RESIDUAL_BOND_STRESS,
                                                      mResidualBondStress);
        myStructure.ConstitutiveLawSetParameterDouble(
                interfaceMaterial, NuTo::Constitutive::eConstitutiveParameter::SLIP_AT_MAX_BOND_STRESS,
                mSlipAtMaxBondStress);
        myStructure.ConstitutiveLawSetParameterDouble(
                interfaceMaterial, NuTo::Constitutive::eConstitutiveParameter::SLIP_AT_RESIDUAL_BOND_STRESS,
                mSlipAtResidualBondStress);

        std::cout << "***********************************" << std::endl;
        std::cout << "**      Interpolation Type       **" << std::endl;
        std::cout << "***********************************" << std::endl;

        int matrixInterpolationType_elastic =
                myStructure.InterpolationTypeCreate(NuTo::Interpolation::eShapeType::TRIANGLE2D);
        myStructure.InterpolationTypeAdd(matrixInterpolationType_elastic, NuTo::Node::COORDINATES,
                                         NuTo::Interpolation::EQUIDISTANT1);
        myStructure.InterpolationTypeAdd(matrixInterpolationType_elastic, NuTo::Node::DISPLACEMENTS,
                                         NuTo::Interpolation::EQUIDISTANT2);

        int matrixInterpolationType_damage =
                myStructure.InterpolationTypeCreate(NuTo::Interpolation::eShapeType::TRIANGLE2D);
        myStructure.InterpolationTypeAdd(matrixInterpolationType_damage, NuTo::Node::COORDINATES,
                                         NuTo::Interpolation::EQUIDISTANT1);
        myStructure.InterpolationTypeAdd(matrixInterpolationType_damage, NuTo::Node::DISPLACEMENTS,
                                         NuTo::Interpolation::EQUIDISTANT2);
        myStructure.InterpolationTypeAdd(matrixInterpolationType_damage, NuTo::Node::NONLOCALEQSTRAIN,
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

        NuTo::FullMatrix<int, Eigen::Dynamic, Eigen::Dynamic> createdGroupIdMatrix = myStructure.ImportFromGmsh(
                meshFileMatrix, NuTo::ElementData::CONSTITUTIVELAWIP, NuTo::IpData::eIpDataType::STATICDATA);
        int groupIdMatrix = createdGroupIdMatrix.GetValue(0, 0);

        myStructure.ElementGroupSetSection(groupIdMatrix, matrixSection);
        myStructure.ElementGroupSetConstitutiveLaw(groupIdMatrix, matrixMaterial);
        myStructure.ElementGroupSetInterpolationType(groupIdMatrix, matrixInterpolationType_damage);
        myStructure.InterpolationTypeSetIntegrationType(matrixInterpolationType_damage,
                                                        NuTo::IntegrationType::IntegrationType2D3NGauss3Ip,
                                                        NuTo::IpData::eIpDataType::STATICDATA);
        myStructure.ElementTotalConvertToInterpolationType(tol, 10);


        std::cout << "***********************************" << std::endl;
        std::cout << "**      Predamage                **" << std::endl;
        std::cout << "***********************************" << std::endl;


        //        int groupNodePredamaged = myStructure.GroupCreate(NuTo::Groups::eGroupId::Nodes);
        //        myStructure.GroupAddNodeCoordinateRange(groupNodePredamaged, 0, matrixLengthX/2 - 1 - 1e-6,
        //        matrixLengthX/2 + 1 + 1e-6);
        //
        //        int groupElePredamaged = myStructure.GroupCreate(NuTo::Groups::eGroupId::Elements);
        //        myStructure.GroupAddElementsFromNodes(groupElePredamaged, groupNodePredamaged, true);
        //
        //        myStructure.ElementGroupSetSection(groupElePredamaged, matrixSection_thin);
        //
        //
        //
        //        auto eleIds = myStructure.GroupGetMemberIds(groupElePredamaged);
        //        for (int iEle = 0; iEle < eleIds.rows(); ++iEle)
        //        {
        //            auto elePtr = myStructure.ElementGetElementPtr(eleIds.at(iEle, 0));
        //
        //            if (elePtr->GetEnumType() == NuTo::Element::eElementType::CONTINUUMELEMENT)
        //            {
        //                for (int iIp = 0; iIp < elePtr->GetNumIntegrationPoints(); ++iIp)
        //                {
        //                    elePtr->GetStaticData(iIp)->AsGradientDamage()->SetKappa(2 * matrixTensileStrength /
        //                    matrixYoungsModulus);
        //                }
        //            }
        //        }


        std::cout << "***********************************" << std::endl;
        std::cout << "**      Boundary Conditions      **" << std::endl;
        std::cout << "***********************************" << std::endl;

        NuTo::FullVector<double, dimension> nodeCoords;

        int groupNodeBCLeft = myStructure.GroupCreate(NuTo::Groups::eGroupId::Nodes);
        myStructure.GroupAddNodeCoordinateRange(groupNodeBCLeft, 0, 0.0 - tol, 0.0 + tol);

        myStructure.ConstraintLinearSetDisplacementNodeGroup(groupNodeBCLeft, directionX, 0);


        int groupNodeBCRight = myStructure.GroupCreate(NuTo::Groups::eGroupId::Nodes);
        myStructure.GroupAddNodeCoordinateRange(groupNodeBCRight, 0, matrixLengthX - tol, matrixLengthX + tol);


        int groupNodeBCBottom = myStructure.GroupCreate(NuTo::Groups::eGroupId::Nodes);
        //        myStructure.GroupAddNodeCoordinateRange(groupNodeBCBottom, 1, 0 - tol, 0 + tol);


        //
        nodeCoords[0] = 0;
        nodeCoords[1] = 42.359;
        myStructure.GroupAddNodeRadiusRange(groupNodeBCBottom, nodeCoords, 0, 1e-2);
        myStructure.ConstraintLinearSetDisplacementNodeGroup(groupNodeBCBottom, directionY, 0);

        std::cout << "***********************************" << std::endl;
        std::cout << "**      Fibre                    **" << std::endl;
        std::cout << "***********************************" << std::endl;

        NuTo::FullMatrix<int, Eigen::Dynamic, Eigen::Dynamic> createdGroupIdFibre = myStructure.ImportFromGmsh(
                meshFileFibre, NuTo::ElementData::CONSTITUTIVELAWIP, NuTo::IpData::eIpDataType::STATICDATA);
        int groupIdFiber = createdGroupIdFibre.GetValue(0, 0);

        myStructure.ElementGroupSetSection(groupIdFiber, fibreSection);
        myStructure.ElementGroupSetConstitutiveLaw(groupIdFiber, fibreMaterial);
        myStructure.ElementGroupSetInterpolationType(groupIdFiber, fibreInterpolationType);
        myStructure.ElementConvertToInterpolationType(groupIdFiber);

        std::cout << "***********************************" << std::endl;
        std::cout << "**      Constraints              **" << std::endl;
        std::cout << "***********************************" << std::endl;


        int groupConstraintNodes = myStructure.GroupCreate(NuTo::Groups::Nodes);
        myStructure.GroupAddNodesFromElements(groupConstraintNodes, groupIdFiber);


        int numNearestNeighbours = 1;

        auto nodeIds = myStructure.GroupGetMemberIds(groupConstraintNodes);
        std::cout << "nodeIds.rows()" << nodeIds.rows() << std::endl;
        for (int iNode = 0; iNode < nodeIds.rows(); ++iNode)
        {
            myStructure.ConstraintLinearEquationNodeToElementCreate(
                    nodeIds.at(iNode, 0), groupIdMatrix, NuTo::Node::eDof::DISPLACEMENTS, numNearestNeighbours);
        }


        std::cout << "***********************************" << std::endl;
        std::cout << "**      Interface                **" << std::endl;
        std::cout << "***********************************" << std::endl;

        auto pairGroupFiberGroupBond =
                myStructure.InterfaceElementsCreate(groupIdFiber, interfaceInterpolationType, fibreInterpolationType);

        int groupEleFiber = pairGroupFiberGroupBond.first;
        std::cout << "groupEleFiber" << groupEleFiber << std::endl;
        int groupEleBond = pairGroupFiberGroupBond.second;
        std::cout << "groupEleBond" << groupEleBond << std::endl;

        myStructure.ElementGroupSetSection(groupEleFiber, fibreSection);
        myStructure.ElementGroupSetSection(groupEleBond, bondSection);

        myStructure.ElementGroupSetConstitutiveLaw(groupEleFiber, fibreMaterial);
        myStructure.ElementGroupSetConstitutiveLaw(groupEleBond, interfaceMaterial);

        std::cout << "***********************************" << std::endl;
        std::cout << "**      Loads                    **" << std::endl;
        std::cout << "***********************************" << std::endl;

        int timeDependentConstraint =
                myStructure.ConstraintLinearSetDisplacementNodeGroup(groupNodeBCRight, directionX, 1);

        std::cout << "***********************************" << std::endl;
        std::cout << "**      Visualization            **" << std::endl;
        std::cout << "***********************************" << std::endl;


        myStructure.AddVisualizationComponent(groupIdMatrix, NuTo::VisualizeBase::DAMAGE);
        myStructure.AddVisualizationComponent(groupIdMatrix, NuTo::VisualizeBase::DISPLACEMENTS);
        myStructure.AddVisualizationComponent(groupIdMatrix, NuTo::VisualizeBase::ENGINEERING_STRAIN);
        myStructure.AddVisualizationComponent(groupIdMatrix, NuTo::VisualizeBase::ENGINEERING_STRESS);
        //        myStructure.AddVisualizationComponent(groupIdMatrix, NuTo::VisualizeBase::NONLOCAL_EQ_STRAIN);


        //        myStructure.AddVisualizationComponent(groupIdMatrix, NuTo::VisualizeBase::SECTION);
        //        myStructure.AddVisualizationComponent(groupIdMatrix, NuTo::VisualizeBase::DAMAGE);
        //        myStructure.AddVisualizationComponent(groupIdMatrix, NuTo::VisualizeBase::LOCAL_EQ_STRAIN);

        myStructure.AddVisualizationComponent(groupEleFiber, NuTo::VisualizeBase::DISPLACEMENTS);
        myStructure.AddVisualizationComponent(groupEleFiber, NuTo::VisualizeBase::ENGINEERING_STRAIN);
        myStructure.AddVisualizationComponent(groupEleFiber, NuTo::VisualizeBase::ENGINEERING_STRESS);

        //        myStructure.AddVisualizationComponent(groupEleBond, NuTo::VisualizeBase::SLIP);
        //        myStructure.ElementGroupExportVtkDataFile(group_matrix_damage, std::string("bla1.vtu"), true);
        //        myStructure.ElementGroupExportVtkDataFile(group_matrix_elastic, std::string("bla2.vtu"), true);


        std::cout << "***********************************" << std::endl;
        std::cout << "**      Solver                   **" << std::endl;
        std::cout << "***********************************" << std::endl;

        myStructure.NodeBuildGlobalDofs();


        myStructure.CalculateMaximumIndependentSets();


        int groupLoad = myStructure.GroupCreate(NuTo::Groups::eGroupId::Nodes);

        nodeCoords[0] = matrixLengthX;
        nodeCoords[1] = 80;
        myStructure.GroupAddNodeRadiusRange(groupLoad, nodeCoords, 0, 1e-2);

        myIntegrationScheme.AddResultGroupNodeForce("myforce", groupNodeBCRight);
        auto nodeID = myStructure.GroupGetMemberIds(groupLoad);
        myIntegrationScheme.AddResultNodeDisplacements("mydisplacements", nodeID.at(0, 0));
        ////
        //        auto interfaceIds = myStructure.GroupGetMemberIds(groupEleBond);
        //        for (int eleId = 0; eleId < interfaceIds.rows(); ++eleId)
        //        {
        //            char myString[200];
        //            sprintf(myString, "bondStress_%d", eleId);
        //            myIntegrationScheme.AddResultElementIpData(myString, interfaceIds.at(eleId,0),
        //            NuTo::IpData::BOND_STRESS);
        //
        //            sprintf(myString, "slip_%d", eleId);
        //            myIntegrationScheme.AddResultElementIpData(myString, interfaceIds.at(eleId,0),
        //            NuTo::IpData::SLIP);
        //        }

        NuTo::FullMatrix<double, 2, 2> timeDependentLoad;
        timeDependentLoad(0, 0) = 0;
        timeDependentLoad(1, 0) = 1.0 * mSimulationTime;

        timeDependentLoad(0, 1) = 0;
        timeDependentLoad(1, 1) = 1.0 * mLoad;

        myIntegrationScheme.AddTimeDependentConstraint(timeDependentConstraint, timeDependentLoad);

        myIntegrationScheme.Solve(mSimulationTime);

        std::string command = "paste " + resultDir + "myforce.dat " + resultDir + "mydisplacements.dat > " + resultDir +
                              "forceDisp.dat";
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
