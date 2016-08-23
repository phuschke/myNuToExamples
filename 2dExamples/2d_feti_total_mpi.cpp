
#include <mpi.h>
#include <iostream>
#include <string.h>
#include <vector>
#include "nuto/mechanics/structures/unstructured/Structure.h"
#include "/opt/nuto/include/eigen3/Eigen/Core"
#include "/opt/nuto/include/eigen3/Eigen/Dense"
#include "ConjugateProjectedGradient.h"
#include <ctime>
class Parameters
{
public:

    static const int mDimension = 2;

    static constexpr double mMatrixYoungsModulus = 200;
    static constexpr double mMatrixPoissonsRatio = 0.2;
    static constexpr double mMatrixNonlocalRadius = 0.5  ;
    static constexpr double mMatrixTensileStrength = 3;
    static constexpr double mMatrixCompressiveStrength = 30;
    static constexpr double mMatrixFractureEnergy = 0.01;
    static constexpr double mMatrixThickness = 1;

    static const NuTo::FullVector<double, 2> mDirectionX;
    static const NuTo::FullVector<double, 2> mDirectionY;

};

const NuTo::FullVector<double, 2> Parameters::mDirectionX = NuTo::FullVector<double, 2>::UnitX();
const NuTo::FullVector<double, 2> Parameters::mDirectionY = NuTo::FullVector<double, 2>::UnitY();


void AssignSection(NuTo::Structure& rStructure, const int rMPIrank)
{
    std::cout << "***********************************" << std::endl;
    std::cout << "**      Section rank = " << rMPIrank << std::endl;
    std::cout << "***********************************" << std::endl;

    int section00 = rStructure.SectionCreate(NuTo::Section::PLANE_STRESS);
    rStructure.SectionSetThickness(section00, Parameters::mMatrixThickness);

    rStructure.ElementTotalSetSection(section00);
}

void AssignMaterial(NuTo::Structure& rStructure, const int rMPIrank)
{
    std::cout << "***********************************" << std::endl;
    std::cout << "**      Material rank = " << rMPIrank << std::endl;
    std::cout << "***********************************" << std::endl;

    int material00 = rStructure.ConstitutiveLawCreate(NuTo::Constitutive::eConstitutiveType::LINEAR_ELASTIC_ENGINEERING_STRESS);

    rStructure.ConstitutiveLawSetParameterDouble(material00, NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, Parameters::mMatrixYoungsModulus);
    rStructure.ConstitutiveLawSetParameterDouble(material00, NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO, Parameters::mMatrixPoissonsRatio);

    rStructure.ElementTotalSetConstitutiveLaw(material00);
//
//    NuTo::FullVector<double, Eigen::Dynamic> myDamageLaw(1);
//    myDamageLaw(0) = NuTo::Constitutive::eDamageLawType::ISOTROPIC_EXPONENTIAL_SOFTENING;
//
//    int matrixMaterial = rStructure.ConstitutiveLawCreate(NuTo::Constitutive::eConstitutiveType::GRADIENT_DAMAGE_ENGINEERING_STRESS);
//    rStructure.ConstitutiveLawSetParameterDouble(matrixMaterial, NuTo::Constitutive::eConstitutiveParameter::YOUNGS_MODULUS, Parameters::mMatrixYoungsModulus);
//    rStructure.ConstitutiveLawSetParameterDouble(matrixMaterial, NuTo::Constitutive::eConstitutiveParameter::POISSONS_RATIO, Parameters::mMatrixPoissonsRatio);
//    rStructure.ConstitutiveLawSetParameterDouble(matrixMaterial, NuTo::Constitutive::eConstitutiveParameter::NONLOCAL_RADIUS, Parameters::mMatrixNonlocalRadius);
//    rStructure.ConstitutiveLawSetParameterDouble(matrixMaterial, NuTo::Constitutive::eConstitutiveParameter::TENSILE_STRENGTH, Parameters::mMatrixTensileStrength);
//    rStructure.ConstitutiveLawSetParameterDouble(matrixMaterial, NuTo::Constitutive::eConstitutiveParameter::COMPRESSIVE_STRENGTH, Parameters::mMatrixCompressiveStrength);
//    rStructure.ConstitutiveLawSetParameterDouble(matrixMaterial, NuTo::Constitutive::eConstitutiveParameter::FRACTURE_ENERGY, Parameters::mMatrixFractureEnergy);
//    rStructure.ConstitutiveLawSetParameterFullVectorDouble(matrixMaterial, NuTo::Constitutive::eConstitutiveParameter::DAMAGE_LAW, myDamageLaw);
//
//    rStructure.ElementTotalSetConstitutiveLaw(matrixMaterial);
}

void ImportMesh(NuTo::Structure& rStructure, const std::string& rMeshFile, const int rMPIrank)
{
    std::cout << "***********************************" << std::endl;
    std::cout << "**      Import Mesh rank = " << rMPIrank << std::endl;
    std::cout << "***********************************" << std::endl;

    NuTo::FullMatrix<int, Eigen::Dynamic, Eigen::Dynamic> createdGroupIdMatrix = rStructure.ImportFromGmsh(rMeshFile, NuTo::ElementData::CONSTITUTIVELAWIP, NuTo::IpData::eIpDataType::STATICDATA);
    int subdomain00 = createdGroupIdMatrix.GetValue(0, 0);

    int interpolationType00 = rStructure.InterpolationTypeCreate(NuTo::Interpolation::eShapeType::QUAD2D);
    rStructure.InterpolationTypeAdd(interpolationType00, NuTo::Node::COORDINATES, NuTo::Interpolation::EQUIDISTANT1);
    rStructure.InterpolationTypeAdd(interpolationType00, NuTo::Node::DISPLACEMENTS, NuTo::Interpolation::EQUIDISTANT2);
//    rStructure.InterpolationTypeAdd(interpolationType00, NuTo::Node::NONLOCALEQSTRAIN, NuTo::Interpolation::EQUIDISTANT1);

    rStructure.ElementGroupSetInterpolationType(subdomain00, interpolationType00);
    rStructure.ElementTotalConvertToInterpolationType(1.e-6, 10);


}



////////////////////////////////
//          MAIN
////////////////////////////////
int main(int argc, char* argv[])
{
    MPI_Init(&argc, &argv);
    int rank;
    int size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);


    std::cout << "Number of tasks=" << size << " My rank=" << rank << std::endl;

    std::string meshFile;

    switch (rank)
    {
    case 0:
        meshFile = "/home/phuschke/develop/nuto/myNutoExamples/MeshFiles/feti00.msh";
        break;
    case 1:
        meshFile = "/home/phuschke/develop/nuto/myNutoExamples/MeshFiles/feti01.msh";
        break;
    default:
        throw NuTo::MechanicsException("Number of tasks > number of subdomains.");
        break;
    }

    NuTo::Structure structure(Parameters::mDimension);
    structure.SetVerboseLevel(10);
    structure.SetShowTime(false);

    ImportMesh(structure, meshFile, rank);
    AssignMaterial(structure, rank);
    AssignSection(structure, rank);

    const int numNodes = structure.GetNumNodes();

    // assemble stiffness matrix
    NuTo::FullVector<double,Eigen::Dynamic>     externalForce;
    NuTo::SparseMatrixCSRGeneral<double>        stiffnessMatrixCSR;
    structure.BuildGlobalCoefficientMatrix0(stiffnessMatrixCSR, externalForce);
    externalForce.setZero();
    NuTo::FullMatrix<double, Eigen::Dynamic, Eigen::Dynamic> stiffnessMatrix = stiffnessMatrixCSR;


    const double searchTol = 1.0e-6;
    // load for rank 1
    if (rank == 1)
    {
        NuTo::FullVector<double,2> nodeCoords;
        nodeCoords[0] = 40;
        nodeCoords[1] = 0;

        int loadNodeGroup = structure.GroupCreate(NuTo::Groups::Nodes);
//        structure.GroupAddNodeRadiusRange(loadNodeGroup, nodeCoords, 0, 1.0e-6);
        structure.GroupAddNodeCoordinateRange(loadNodeGroup, 0, 40 - searchTol, 40 + searchTol);
        auto loadNodes = structure.GroupGetMemberIds(loadNodeGroup);


        for (int iNode = 0; iNode < loadNodes.rows(); ++iNode)
        {
            NuTo::FullVector<int, -1> displacementDofs;
            structure.NodeGetDisplacementDofs(loadNodes.at(iNode, 0), displacementDofs);

            externalForce(displacementDofs.at(0, 0)) = 0;
            externalForce(displacementDofs.at(1, 0)) = 1;

        }

//        std::cout << "externalForce" << externalForce << std::endl;
    }




    // get all nodes on the interface
    int groupNodeInterface = structure.GroupCreate(NuTo::Groups::Nodes);
    structure.GroupAddNodeCoordinateRange(groupNodeInterface, 0, 20.0 - searchTol, 20.0 + searchTol);

    // get all nodes on the boundary
    int groupNodeBoundary = structure.GroupCreate(NuTo::Groups::Nodes);
    structure.GroupAddNodeCoordinateRange(groupNodeBoundary, 0, 0.0  - searchTol, 0.0  + searchTol);

    int numInterfaceAndBoundaryNodes = structure.GroupGetNumMembers(groupNodeInterface) + structure.GroupGetNumMembers(groupNodeBoundary);
    MPI_Bcast(&numInterfaceAndBoundaryNodes, 1,MPI_INT,0, MPI_COMM_WORLD);

    auto nodesInterface = structure.GroupGetMemberIds(groupNodeInterface);
    auto nodesBoundary  = structure.GroupGetMemberIds(groupNodeBoundary);

    Eigen::MatrixXd connectivityMatrix = Eigen::MatrixXd::Zero(Parameters::mDimension*numInterfaceAndBoundaryNodes, stiffnessMatrix.cols());

    for (int i = 0; i < nodesInterface.rows(); ++i)
    {
        NuTo::FullVector<double,2> nodeCoords;
        nodeCoords[0] = 20;
        nodeCoords[1] = i*0.125;

        int groupNodeTmp = structure.GroupCreate(NuTo::Groups::Nodes);
        structure.GroupAddNodeRadiusRange(groupNodeTmp, nodeCoords, 0, 1.0e-6);
        auto tmpNode = structure.GroupGetMemberIds(groupNodeTmp);


        NuTo::FullVector<int, -1> displacementDofs;
        structure.NodeGetDisplacementDofs(tmpNode.at(0, 0), displacementDofs);


        int index = 2*i;
        connectivityMatrix(index,   displacementDofs.at(0,0)) = ( rank==0 ? 1 : -1 );
        connectivityMatrix(index+1, displacementDofs.at(1,0)) = ( rank==0 ? 1 : -1 );

    }


    for (int i = 0; i < nodesBoundary.rows(); ++i)
    {
        NuTo::FullVector<int, -1> displacementDofs;
        structure.NodeGetDisplacementDofs(nodesBoundary.at(i, 0), displacementDofs);

        int index = 2*(i + nodesInterface.rows());
        connectivityMatrix(index,   displacementDofs.at(0,0)) = 1.0;
        connectivityMatrix(index+1, displacementDofs.at(1,0)) = 1.0;
    }


//    if(rank == 1)
//    {
//        std::cout << "connec" << std::endl << connectivityMatrix << std::endl;
//    }

    std::cout << "***********************************" << std::endl;
    std::cout << "**      FETI rank = " << rank        << std::endl;
    std::cout << "***********************************" << std::endl;

    const int numDofs           = structure.GetNumDofs();
    const int numInterfaceEqs   = connectivityMatrix.rows();



    // compute the total number of rigid body modes
    int numRigidBodyModesLocal    = 3;
    Eigen::MatrixXd rigidBodyModes = Eigen::MatrixXd::Zero(numDofs,numRigidBodyModesLocal);


    for (int iNode = 0; iNode < numNodes; ++iNode)
    {
        const int startRow = 2 * iNode;
        constexpr int startCol = 0;

        const Eigen::Matrix<double, 2, 1> coordinates = structure.NodeGetNodePtr(iNode)->GetCoordinates2D();

        rigidBodyModes.block(startRow, startCol, Parameters::mDimension, 3) <<  1.,     0.,    -coordinates.at(1, 0),
                                                                                0.,     1.,     coordinates.at(0, 0);
    }

    int numRigidBodyModesGlobal   = 0;


    MPI_Allreduce(&numRigidBodyModesLocal, &numRigidBodyModesGlobal, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);


    std::cout << "numRigidBodyModesGlobal" << std::endl << numRigidBodyModesGlobal << std::endl;




    Eigen::MatrixXd interfaceRigidBodyModesLocal    = connectivityMatrix * rigidBodyModes;
    Eigen::VectorXd rigidBodyForceVectorLocal       = rigidBodyModes.transpose() * externalForce;


    std::clock_t begin_time = std::clock();
    Eigen::LDLT<Eigen::MatrixXd> solverLDLT(stiffnessMatrix);
    if (rank==0)
        std::cout << float( std::clock () - begin_time ) /  CLOCKS_PER_SEC<< std::endl;

    begin_time = std::clock();
    Eigen::VectorXd displacementGapLocal    = connectivityMatrix * solverLDLT.solve(externalForce);
    if (rank==0)
            std::cout << float( std::clock () - begin_time ) /  CLOCKS_PER_SEC << std::endl;
    Eigen::VectorXd displacementGapSum      = Eigen::VectorXd::Zero(numInterfaceEqs);

    MPI_Reduce(displacementGapLocal.data(), displacementGapSum.data(), displacementGapLocal.size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);


    Eigen::VectorXd rigidBodyForceVectorGlobal      = Eigen::VectorXd::Zero(numRigidBodyModesGlobal);
    Eigen::MatrixXd interfaceRigidBodyModesGlobal   = Eigen::MatrixXd::Zero(numInterfaceEqs,numRigidBodyModesGlobal);

//    std::cout << "numRigidBodyModesGlobal before rank " <<  rank << std::endl << numRigidBodyModesGlobal << std::endl;
//    std::cout << "interfaceRigidBodyModesLocal before rank " <<  rank << std::endl << interfaceRigidBodyModesLocal << std::endl;
//    std::cout << "interfaceRigidBodyModesGlobal before rank " <<  rank  << std::endl << interfaceRigidBodyModesGlobal << std::endl;





    MPI_Gather(interfaceRigidBodyModesLocal.data(), interfaceRigidBodyModesLocal.size(), MPI_DOUBLE, interfaceRigidBodyModesGlobal.data(), interfaceRigidBodyModesLocal.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(rigidBodyForceVectorLocal.data(), rigidBodyForceVectorLocal.size(), MPI_DOUBLE, rigidBodyForceVectorGlobal.data(), rigidBodyForceVectorLocal.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);




//    if (rank == 0)
//    {
//        std::cout << "interfaceRigidBodyModesGlobal after rank " <<  rank  << std::endl << interfaceRigidBodyModesGlobal << std::endl;
//        std::cout << "rigidBodyForceVectorGlobal after rank " <<  rank  << std::endl << rigidBodyForceVectorGlobal << std::endl;
//    }

//    Eigen::MatrixXd interfaceFlexibilityLocal    = connectivityMatrix * PseudoInverse * connectivityMatrix.transpose();
    begin_time = std::clock();
    Eigen::MatrixXd interfaceFlexibilityLocal    = connectivityMatrix * solverLDLT.solve(connectivityMatrix.transpose());
    if (rank==0)
            std::cout << float( std::clock () - begin_time ) /  CLOCKS_PER_SEC<< std::endl;
    Eigen::MatrixXd interfaceFlexibilitySum      = Eigen::MatrixXd::Zero(numInterfaceEqs, numInterfaceEqs);

    MPI_Reduce(interfaceFlexibilityLocal.data(), interfaceFlexibilitySum.data(), interfaceFlexibilityLocal.size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    Eigen::VectorXd lambda      = Eigen::VectorXd::Zero(numInterfaceEqs);
    Eigen::VectorXd alphaGlobal = Eigen::VectorXd::Zero(numRigidBodyModesGlobal);

    if (rank == 0)
    {
        std::cout << "***********************************" << std::endl;
        std::cout << "**   Conjugate Projected Gradient**" << std::endl;
        std::cout << "***********************************" << std::endl;

//        std::cout << "interfaceRigidBodyModesLocal" << std::endl << interfaceRigidBodyModesLocal << std::endl;
//        std::cout << "interfaceRigidBodyModesGlobal" << std::endl << interfaceRigidBodyModesGlobal << std::endl;
//        std::cout << "connectivityMatrix" << std::endl << connectivityMatrix << std::endl;

        Eigen::MatrixXd projection =        Eigen::MatrixXd::Identity(numInterfaceEqs,numInterfaceEqs)
                                        -   interfaceRigidBodyModesGlobal * (interfaceRigidBodyModesGlobal.transpose() * interfaceRigidBodyModesGlobal).inverse() * interfaceRigidBodyModesGlobal.transpose();

        lambda  = interfaceRigidBodyModesGlobal * (interfaceRigidBodyModesGlobal.transpose() * interfaceRigidBodyModesGlobal).inverse() * rigidBodyForceVectorGlobal;

        Eigen::ConjugateProjectedGradient<Eigen::MatrixXd,Eigen::Lower,Eigen::IdentityPreconditioner> CPGSolver(interfaceFlexibilitySum);
        CPGSolver.setTolerance(1.0e-14);

        std::cout << "egienvalues" << interfaceFlexibilitySum.eigenvalues() << std::endl;

        CPGSolver._solveWithGuess(displacementGapSum, lambda, projection);

//        std::cout << "projection" << std::endl << projection<< std::endl;
//        std::cout << "lambda" << std::endl <<lambda<< std::endl;
//        std::cout << "displacementGapSum" << std::endl <<displacementGapSum<< std::endl;



//
//        Eigen::VectorXd r       = displacementGapSum - interfaceFlexibilitySum*lambda;
//        Eigen::VectorXd w       = projection * r;
//        Eigen::VectorXd w_old   = w;
//        Eigen::VectorXd p       = w;
//
//        for (int iter = 0; iter < 100; ++iter)
//        {
//            double gamma = w_old.dot(w_old) / p.dot(interfaceFlexibilitySum*p);
//            lambda += gamma * p;
//            r +=   - gamma * interfaceFlexibilitySum * p;
//            w = projection * r;
//
//            double beta = w.dot(w) / w_old.dot(w_old);
//            p = w + beta*p;
//            w_old = w;
//        }

        // Conjugate projected gradient method with re-orthogonalization
//        Eigen::VectorXd r = displacementGapSum - interfaceFlexibilitySum * lambda;
//        Eigen::VectorXd w = projection * r;
//        Eigen::VectorXd w_old = w;
//        const int maxIter = interfaceFlexibilitySum.rows();
//        Eigen::MatrixXd p = Eigen::MatrixXd::Zero(w.rows(), maxIter+1);
//
//
//        p.col(0) = w;
//
//
//        for (int iter = 0; iter < maxIter; ++iter)
//        {
//
//            double gamma = w_old.dot(w_old) / p.col(iter).dot(interfaceFlexibilitySum * p.col(iter));
//
//            lambda += gamma * p.col(iter);
//
//            r += -gamma * interfaceFlexibilitySum * p.col(iter);
//
//            w = projection * r;
//
//            p.col(iter+1) = w;
//            for (int i = 0; i < iter+1; ++i)
//            {
//                p.col(iter+1) -= w.dot(interfaceFlexibilitySum * p.col(i)) / p.col(i).dot(interfaceFlexibilitySum * p.col(i)) * p.col(i);
//            }
//            w_old = w;
//        }





        if (rank==0)
           {
               std::cout << "cpg max iterations" << CPGSolver.maxIterations() << std::endl;
               std::cout << "cpg tolerance" << CPGSolver.tolerance() << std::endl;
               std::cout << "cpg iterations" << CPGSolver.iterations() << std::endl;
               std::cout << "cpg error" << CPGSolver.error() << std::endl;
//               std::cout << "lambda" << lambda << std::endl;
           }

        alphaGlobal = (interfaceRigidBodyModesGlobal.transpose() * interfaceRigidBodyModesGlobal).inverse() * interfaceRigidBodyModesGlobal.transpose() * (displacementGapSum - interfaceFlexibilitySum*lambda);

        std::cout << "alphaGlobal" << alphaGlobal << std::endl;
    }

    MPI_Bcast(lambda.data(), lambda.size(),MPI_DOUBLE,0, MPI_COMM_WORLD);

    Eigen::VectorXd alphaLocal = Eigen::VectorXd::Zero(numRigidBodyModesLocal);
    MPI_Scatter(alphaGlobal.data(), alphaLocal.size(), MPI_DOUBLE, alphaLocal.data(), alphaLocal.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);

    Eigen::VectorXd displacement = Eigen::VectorXd::Zero(numDofs);

    displacement.head(structure.GetNumActiveDofs()) = solverLDLT.solve(externalForce - connectivityMatrix.transpose() * lambda) - rigidBodyModes * alphaLocal;

//    std::cout << "displacement" << std::endl << displacement << std::endl;


    std::cout << "***********************************" << std::endl;
    std::cout << "**      Visualization rank = " << rank        << std::endl;
    std::cout << "***********************************" << std::endl;


    int groupAllElements = structure.GroupCreate(NuTo::Groups::Elements);
    structure.GroupAddElementsTotal(groupAllElements);
    structure.AddVisualizationComponent(groupAllElements, NuTo::VisualizeBase::DISPLACEMENTS);
//    structure.AddVisualizationComponent(groupAllElements, NuTo::VisualizeBase::DAMAGE);
//    structure.AddVisualizationComponent(groupAllElements, NuTo::VisualizeBase::NONLOCAL_EQ_STRAIN);
    structure.AddVisualizationComponent(groupAllElements, NuTo::VisualizeBase::ENGINEERING_STRESS);
    structure.AddVisualizationComponent(groupAllElements, NuTo::VisualizeBase::ENGINEERING_STRAIN);

    char outputFile[200];
    sprintf(outputFile, "/home/phuschke/structure_rank_%03d_0.vtk", rank);
    structure.ExportVtkDataFileElements(outputFile);


//    if(rank==0)
//        {
//            std::cout << "connec" << std::endl << connectivityMatrix*displacement.head(structure.GetNumActiveDofs()) << std::endl;
//
//        }


    for (int iNode = 0; iNode < numNodes; ++iNode)
    {

        NuTo::FullVector<int, -1> displacementDofs;
        structure.NodeGetDisplacementDofs(iNode, displacementDofs);

        NuTo::FullVector<double, Eigen::Dynamic> nodeDisplacements(Parameters::mDimension);
        nodeDisplacements(0) = displacement(displacementDofs(0));
        nodeDisplacements(1) = displacement(displacementDofs(1));

        structure.NodeSetDisplacements(iNode, nodeDisplacements);

    }
//
//    for (int iNode = 0; iNode < numNodes; ++iNode)
//    {
//
//        NuTo::FullVector<double, Eigen::Dynamic> nodeDisplacements(Parameters::mDimension);
//        nodeDisplacements(0) = displacement(2*iNode);
//        nodeDisplacements(1) = displacement(2*iNode +1);
//
//        structure.NodeSetDisplacements(iNode, displacement.at(iNode,0));
//
//    }


    sprintf(outputFile, "/home/phuschke/structure_rank_%03d_1.vtk", rank);
    structure.ExportVtkDataFileElements(outputFile);




//    if(rank==1)
//    {
//        std::cout << "connec" << std::endl << connectivityMatrix*displacement.head(structure.GetNumActiveDofs()) << std::endl;
//
//        structure.Info();
//    }
//
//    if(rank==0)
//    {
//        sleep(2);
//        std::cout << connectivityMatrix << std::endl;
//        structure.Info();
//
//
//        for (int i = 0; i < numDofs; ++i)
//        {
//            std::cout << "dof " << i << " disp " << displacement.at(i,0) << std::endl;
//        }
//    }

    // Tell the MPI library to release all resources it is using:

    if (rank==0)
    {
        std::cout << "size local problem: " << stiffnessMatrix.rows() << std::endl;
        std::cout << "size global problem: " << interfaceFlexibilitySum.rows() << std::endl;
    }
    MPI_Finalize();
}
