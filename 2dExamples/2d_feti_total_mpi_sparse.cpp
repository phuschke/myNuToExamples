
#include <mpi.h>
#include <iostream>
#include <string.h>
#include <vector>
#include "nuto/mechanics/structures/unstructured/Structure.h"
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/SparseQR>
#include "ConjugateProjectedGradient.h"
#include <ctime>
#include <cstdlib>
#include <chrono>
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



void AssembleStiffnesMatrix(NuTo::Structure& structure, Eigen::SparseMatrix<double>& stiffnessMatrixSparse)
{
    // assemble stiffness matrix
    NuTo::StructureOutputBlockMatrix stiffnessMatrix = structure.BuildGlobalHessian0();
    NuTo::SparseMatrixCSRGeneral<double> stiffnessMatrixCSR(stiffnessMatrix.JJ(NuTo::Node::DISPLACEMENTS, NuTo::Node::DISPLACEMENTS));

    std::cout << "stiffnessMatrixCSR.GetNumEntries()" << stiffnessMatrixCSR.GetNumEntries() << std::endl;
    std::cout << "stiffnessMatrixCSR.GetNumColumns()" << stiffnessMatrixCSR.GetColumns().size() << std::endl;
    std::cout << "stiffnessMatrixCSR.GetNumRows()" << stiffnessMatrixCSR.GetRowIndex().size() << std::endl;

    std::vector<Eigen::Triplet<double>> tripletList;

    std::vector<double> val = stiffnessMatrixCSR.GetValues();
    std::vector<int> colInd = stiffnessMatrixCSR.GetColumns();
    std::vector<int> rowInd = stiffnessMatrixCSR.GetRowIndex();

    for (unsigned i = 0; i < rowInd.size() - 1; ++i)
    {
        for (int k = rowInd[i]; k < rowInd[i + 1]; ++k)
            tripletList.push_back(Eigen::Triplet<double>(i, colInd[k], val[k]));
    }


    stiffnessMatrixSparse.setFromTriplets(tripletList.begin(), tripletList.end());

}


void my_conjugate_projected_gradient(const Eigen::MatrixXd& mat, const Eigen::VectorXd& rhs, Eigen::VectorXd& x, const Eigen::MatrixXd& projectionMatrix,
                        const Eigen::MatrixXd& precond, int& iters,
                        double& tol_error)
{

//  using std::sqrt;
//  using std::abs;

//  double tol = tol_error;
//  int maxIters = iters;

//  int n = mat.cols();

//  Eigen::VectorXd residual = rhs - mat * x; //initial residual

//  double rhsNorm2 = rhs.squaredNorm();
//  if(rhsNorm2 == 0)
//  {
//    x.setZero();
//    iters = 0;
//    tol_error = 0;
//    return;
//  }

//  double threshold = tol*tol*rhsNorm2;

//  double residualNorm2 = residual.squaredNorm();
//  if (residualNorm2 < threshold)
//  {
//    iters = 0;
//    tol_error = sqrt(residualNorm2 / rhsNorm2);
//    return;
//  }


//  Eigen::VectorXd w = projectionMatrix * precond.ldlt().solve(residual);      //initial search direction
//  Eigen::VectorXd p = w;

////  MatrixXd pInitial(p.rows(), maxIters+1);  // test for reorthogonlization
////  pInitial.col(0) = w;



//  Eigen::VectorXd z(n), tmp(n);
//  double absNew = w.dot(w);  // the square of the absolute value of r scaled by invM
//  double wNorm2;
//  int i = 0;

//  while(i < maxIters)
//  {
//    tmp.noalias() = mat * p;              // the bottleneck of the algorithm

//    double alpha = absNew / p.dot(tmp);   // the amount we travel on dir
//    x += alpha * p;                       // update solution
//    residual -= alpha * tmp;              // update residue

//    residualNorm2 = residual.squaredNorm();
//    if(residualNorm2 < threshold)
//      break;

//    w = projectionMatrix * residual;

////   wNorm2 = w.squaredNorm();
////    if (wNorm2 < threshold)
////        break;

//    z = precond.ldlt().solve(w);          // approximately solve for "A z = residual"

//    double absOld = absNew;
//    absNew = w.dot(z);     // update the absolute value of r
//    double beta = absNew / absOld;            // calculate the Gram-Schmidt value used to create the new search direction
//    p = z + beta * p;                             // update search direction

////    pInitial.col(i+1) = z;
////    std::cout << "bla" << std::endl;
////    for(int k = 0; k < i; ++k)
////    {
////        pInitial.col(i+1) -= z.dot(mat*pInitial.col(k))/pInitial.col(k).dot(mat*pInitial.col(k))*pInitial.col(k);
////        std::cout << "bla" << std::endl;
////    }
////    p = pInitial.col(i+1);



//    i++;
//  }
//  tol_error = sqrt(residualNorm2 / rhsNorm2);
//  iters = i;
}

void AssenbleExternalForceVector(const double searchTol, NuTo::FullVector<double, Eigen::Dynamic>& externalForce, NuTo::Structure& structure, double xCoord)
{
    NuTo::FullVector<double, 2> nodeCoords;
    nodeCoords[0] = xCoord;
    nodeCoords[1] = 0;
    int loadNodeGroup = structure.GroupCreate(NuTo::Groups::Nodes);
    //        structure.GroupAddNodeRadiusRange(loadNodeGroup, nodeCoords, 0, 1.0e-6);
    structure.GroupAddNodeCoordinateRange(loadNodeGroup, 0, xCoord - searchTol, xCoord + searchTol);
    auto loadNodes = structure.GroupGetMemberIds(loadNodeGroup);
    for (int iNode = 0; iNode < loadNodes.rows(); ++iNode)
    {
        NuTo::FullVector<int, -1> displacementDofs;
        structure.NodeGetDisplacementDofs(loadNodes.at(iNode, 0), displacementDofs);
        externalForce(displacementDofs.at(0, 0)) = 0;
        externalForce(displacementDofs.at(1, 0)) = 0.1;
    }
}

void AssenbleExternalForceVectorThreePointBending(const double searchTol, NuTo::FullVector<double, Eigen::Dynamic>& externalForce, NuTo::Structure& structure, double xCoord)
{
    NuTo::FullVector<double, 2> nodeCoords;
    nodeCoords[0] = 30;
    nodeCoords[1] = 10;
    int loadNodeGroup = structure.GroupCreate(NuTo::Groups::Nodes);
    structure.GroupAddNodeRadiusRange(loadNodeGroup, nodeCoords, 0, 1.0e-6);

    auto loadNodes = structure.GroupGetMemberIds(loadNodeGroup);
    for (int iNode = 0; iNode < loadNodes.rows(); ++iNode)
    {
        NuTo::FullVector<int, -1> displacementDofs;
        structure.NodeGetDisplacementDofs(loadNodes.at(iNode, 0), displacementDofs);
        externalForce(displacementDofs.at(0, 0)) = 0;
        externalForce(displacementDofs.at(1, 0)) = 1;
    }
}


Eigen::SparseMatrix<double> AssembleConnectivityMatrixThreePointBending(const double searchTol, int rank, NuTo::Structure& structure, int nodesOnEdge)
{
    // get all nodes on the interface 01
    int groupNodeInterface01 = structure.GroupCreate(NuTo::Groups::Nodes);
    structure.GroupAddNodeCoordinateRange(groupNodeInterface01, 0, 20.0 - searchTol, 20.0 + searchTol);

    // get all nodes on the interface 12
    int groupNodeInterface12 = structure.GroupCreate(NuTo::Groups::Nodes);
    structure.GroupAddNodeCoordinateRange(groupNodeInterface12, 0, 40.0 - searchTol, 40.0 + searchTol);

    int numInterface01 = structure.GroupGetNumMembers(groupNodeInterface01);
    int numInterface12 = structure.GroupGetNumMembers(groupNodeInterface12);

    int numInterface01Sum = 0;
    MPI_Allreduce(&numInterface01, &numInterface01Sum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    int numInterface12Sum = 0;
    MPI_Allreduce(&numInterface12, &numInterface12Sum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    // the domain is divided in such a way that each node on an interface is shared by exactly 2 subdomains

    numInterface01Sum /= 2;
    numInterface12Sum /= 2;

    // get all nodes on the boundary
    NuTo::FullVector<double, 2> nodeCoords;
    nodeCoords[0] = 0.0;
    nodeCoords[1] = 0.0;
    int groupNodeBoundaryLeft = structure.GroupCreate(NuTo::Groups::Nodes);
    structure.GroupAddNodeRadiusRange(groupNodeBoundaryLeft, nodeCoords, 0, 1);
    int numBoundaryNodesLeft = structure.GroupGetNumMembers(groupNodeBoundaryLeft);

    nodeCoords[0] = 60.0;
    nodeCoords[1] = 0.0;
    int groupNodeBoundaryRight = structure.GroupCreate(NuTo::Groups::Nodes);
    structure.GroupAddNodeRadiusRange(groupNodeBoundaryRight, nodeCoords, 0, 1);
    int numBoundaryNodesRight = structure.GroupGetNumMembers(groupNodeBoundaryRight);

    int numBoundaryNodesLeftSum = 0;
    MPI_Allreduce(&numBoundaryNodesLeft, &numBoundaryNodesLeftSum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    int numBoundaryNodesRightSum = 0;
    MPI_Allreduce(&numBoundaryNodesRight, &numBoundaryNodesRightSum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    auto nodesBoundaryLeft      = structure.GroupGetMemberIds(groupNodeBoundaryLeft);
    auto nodesBoundaryRight     = structure.GroupGetMemberIds(groupNodeBoundaryRight);

    const int numLagrangeMultiplier = Parameters::mDimension * (numBoundaryNodesRightSum + numBoundaryNodesLeftSum + numInterface01Sum + numInterface12Sum);

    Eigen::SparseMatrix<double> connectivityMatrix(numLagrangeMultiplier, structure.GetNumDofs(NuTo::Node::DISPLACEMENTS));

    if (rank==0)
    {
        std::cout << "numLagrangeMultiplier \n" << numLagrangeMultiplier << std::endl;
        std::cout << "numBoundaryNodesLeftSum  \n" << numBoundaryNodesLeftSum << std::endl;
        std::cout << "numBoundaryNodesRightSum  \n" << numBoundaryNodesRightSum << std::endl;
        std::cout << "numInterface01Sum \n" << numInterface01Sum << std::endl;
        std::cout << "numInterface12Sum \n" << numInterface12Sum << std::endl;
    }

    for (int i = 0; i < numInterface01; ++i)
    {
        NuTo::FullVector<double, 2> nodeCoords;
        nodeCoords[0] = 20;
        nodeCoords[1] = i * 10.0/(2*(nodesOnEdge-1));
        int groupNodeTmp = structure.GroupCreate(NuTo::Groups::Nodes);
        structure.GroupAddNodeRadiusRange(groupNodeTmp, nodeCoords, 0, 1.0e-6);
        auto tmpNode = structure.GroupGetMemberIds(groupNodeTmp);
        NuTo::FullVector<int, -1> displacementDofs;
        structure.NodeGetDisplacementDofs(tmpNode.at(0, 0), displacementDofs);
        int index = 2 * i;
        connectivityMatrix.insert(index, displacementDofs.at(0, 0)) = (rank == 0 ? 1 : -1);
        connectivityMatrix.insert(index + 1, displacementDofs.at(1, 0)) = (rank == 0 ? 1 : -1);
    }



    for (int i = 0; i < numInterface12; ++i)
    {
        NuTo::FullVector<double, 2> nodeCoords;
        nodeCoords[0] = 40;
        nodeCoords[1] = i * 10.0 / (2 * (nodesOnEdge - 1));
        int groupNodeTmp = structure.GroupCreate(NuTo::Groups::Nodes);
        structure.GroupAddNodeRadiusRange(groupNodeTmp, nodeCoords, 0, 1.0e-6);
        auto tmpNode = structure.GroupGetMemberIds(groupNodeTmp);
        NuTo::FullVector<int, -1> displacementDofs;
        structure.NodeGetDisplacementDofs(tmpNode.at(0, 0), displacementDofs);
        int index = 2 * (i + numInterface01Sum);
        connectivityMatrix.insert(index, displacementDofs.at(0, 0)) = (rank == 1 ? 1 : -1);
        connectivityMatrix.insert(index + 1, displacementDofs.at(1, 0)) = (rank == 1 ? 1 : -1);
    }



    for (int i = 0; i < numBoundaryNodesLeft; ++i)
    {
        NuTo::FullVector<int, -1> displacementDofs;
        structure.NodeGetDisplacementDofs(nodesBoundaryLeft.at(i, 0), displacementDofs);
        int index = 2 * (i + numInterface01Sum + numInterface12Sum);
        connectivityMatrix.insert(index, displacementDofs.at(0, 0)) = 1.0;
        connectivityMatrix.insert(index + 1, displacementDofs.at(1, 0)) = 1.0;
    }

    for (int i = 0; i < numBoundaryNodesRight; ++i)
    {
        NuTo::FullVector<int, -1> displacementDofs;
        structure.NodeGetDisplacementDofs(nodesBoundaryRight.at(i, 0), displacementDofs);
        int index = 2 * (i + numInterface01Sum + numInterface12Sum + numBoundaryNodesLeft);
        connectivityMatrix.insert(index, displacementDofs.at(0, 0)) = 1.0;
        connectivityMatrix.insert(index + 1, displacementDofs.at(1, 0)) = 1.0;
    }




    return connectivityMatrix;
}


Eigen::SparseMatrix<double> AssembleConnectivityMatrix(const double searchTol, int rank, NuTo::Structure& structure, int nodesOnEdge)
{
    // get all nodes on the interface 01
    int groupNodeInterface01 = structure.GroupCreate(NuTo::Groups::Nodes);
    structure.GroupAddNodeCoordinateRange(groupNodeInterface01, 0, 20.0 - searchTol, 20.0 + searchTol);

    // get all nodes on the interface 12
    int groupNodeInterface12 = structure.GroupCreate(NuTo::Groups::Nodes);
    structure.GroupAddNodeCoordinateRange(groupNodeInterface12, 0, 40.0 - searchTol, 40.0 + searchTol);

    int numInterface01 = structure.GroupGetNumMembers(groupNodeInterface01);
    int numInterface12 = structure.GroupGetNumMembers(groupNodeInterface12);

    int numInterface01Sum = 0;
    MPI_Allreduce(&numInterface01, &numInterface01Sum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    int numInterface12Sum = 0;
    MPI_Allreduce(&numInterface12, &numInterface12Sum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    // the domain is divided in such a way that each node on an interface is shared by exactly 2 subdomains

    numInterface01Sum /= 2;
    numInterface12Sum /= 2;


    // get all nodes on the boundary
    int groupNodeBoundary = structure.GroupCreate(NuTo::Groups::Nodes);
    structure.GroupAddNodeCoordinateRange(groupNodeBoundary, 0, 0.0 - searchTol, 0.0 + searchTol);
    int numBoundaryNodes = structure.GroupGetNumMembers(groupNodeBoundary);

    int numBoundaryNodesSum = 0;
    MPI_Allreduce(&numBoundaryNodes, &numBoundaryNodesSum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    auto nodesBoundary      = structure.GroupGetMemberIds(groupNodeBoundary);

    const int numLagrangeMultiplier = Parameters::mDimension * (numBoundaryNodesSum + numInterface01Sum + numInterface12Sum);

    Eigen::SparseMatrix<double> connectivityMatrix(numLagrangeMultiplier, structure.GetNumDofs(NuTo::Node::DISPLACEMENTS));


    if (rank==0)
    {
        std::cout << "numLagrangeMultiplier \n" << numLagrangeMultiplier << std::endl;
        std::cout << "numBoundaryNodes  \n" << numBoundaryNodes << std::endl;
        std::cout << "numInterface01Sum \n" << numInterface01Sum << std::endl;
        std::cout << "numInterface12Sum \n" << numInterface12Sum << std::endl;
    }

    for (int i = 0; i < numInterface01; ++i)
    {
        NuTo::FullVector<double, 2> nodeCoords;
        nodeCoords[0] = 20;
        nodeCoords[1] = i * 10.0/(2*(nodesOnEdge-1));
        int groupNodeTmp = structure.GroupCreate(NuTo::Groups::Nodes);
        structure.GroupAddNodeRadiusRange(groupNodeTmp, nodeCoords, 0, 1.0e-6);
        auto tmpNode = structure.GroupGetMemberIds(groupNodeTmp);
        NuTo::FullVector<int, -1> displacementDofs;
        structure.NodeGetDisplacementDofs(tmpNode.at(0, 0), displacementDofs);
        int index = 2 * i;
        connectivityMatrix.insert(index, displacementDofs.at(0, 0)) = (rank == 0 ? 1 : -1);
        connectivityMatrix.insert(index + 1, displacementDofs.at(1, 0)) = (rank == 0 ? 1 : -1);
    }

    for (int i = 0; i < numInterface12; ++i)
    {
        NuTo::FullVector<double, 2> nodeCoords;
        nodeCoords[0] = 40;
        nodeCoords[1] = i * 10.0 / (2 * (nodesOnEdge - 1));
        int groupNodeTmp = structure.GroupCreate(NuTo::Groups::Nodes);
        structure.GroupAddNodeRadiusRange(groupNodeTmp, nodeCoords, 0, 1.0e-6);
        auto tmpNode = structure.GroupGetMemberIds(groupNodeTmp);
        NuTo::FullVector<int, -1> displacementDofs;
        structure.NodeGetDisplacementDofs(tmpNode.at(0, 0), displacementDofs);
        int index = 2 * (i + numInterface01Sum);
        connectivityMatrix.insert(index, displacementDofs.at(0, 0)) = (rank == 1 ? 1 : -1);
        connectivityMatrix.insert(index + 1, displacementDofs.at(1, 0)) = (rank == 1 ? 1 : -1);
    }

    for (int i = 0; i < numBoundaryNodes; ++i)
    {
        NuTo::FullVector<int, -1> displacementDofs;
        structure.NodeGetDisplacementDofs(nodesBoundary.at(i, 0), displacementDofs);
        int index = 2 * (i + numInterface01Sum + numInterface12Sum);
        connectivityMatrix.insert(index, displacementDofs.at(0, 0)) = 1.0;
        connectivityMatrix.insert(index + 1, displacementDofs.at(1, 0)) = 1.0;
    }

    return connectivityMatrix;
}

void AssembleRigidBodyModes(const int numNodes, NuTo::Structure& structure, Eigen::MatrixXd& rigidBodyModes)
{
    for (int iNode = 0; iNode < numNodes; ++iNode)
    {
        const int startRow = 2 * iNode;
        constexpr int startCol = 0;
        const Eigen::Matrix<double, 2, 1> coordinates = structure.NodeGetNodePtr(iNode)->Get(NuTo::Node::COORDINATES);
        rigidBodyModes.block(startRow, startCol, Parameters::mDimension, 3) << 1., 0., -coordinates.at(1, 0), 0., 1., coordinates.at(0, 0);
    }
}


void CreateGmshFile(const std::string& rGmshFile, double x_start, double y_start, double x_end, double y_end, int nodesOnEdge)
{
    std::ofstream GmshFile;

    GmshFile.open(rGmshFile + ".geo");

    GmshFile << "Mesh.RecombineAll = 1;" << std::endl;
    GmshFile << "x_start = " << x_start << ";" << std::endl;
    GmshFile << "y_start = " << y_start << ";" << std::endl;
    GmshFile << "x_end   = " << x_end   << ";" << std::endl;
    GmshFile << "y_end   = " << y_end   << ";" << std::endl;
    GmshFile << "meshSize       = 10;" << std::endl;

    GmshFile << "Point(1)  = {x_start   ,   y_start     ,   0,  meshSize};" << std::endl;
    GmshFile << "Point(2)  = {x_end     ,   y_start     ,   0,  meshSize};" << std::endl;
    GmshFile << "Point(3)  = {x_end     ,   y_end       ,   0,  meshSize};" << std::endl;
    GmshFile << "Point(4)  = {x_start   ,   y_end       ,   0,  meshSize};" << std::endl;

    GmshFile << "l1  = newreg; Line(l1)  = {1,2};" << std::endl;
    GmshFile << "l2  = newreg; Line(l2)  = {2,3};" << std::endl;
    GmshFile << "l3  = newreg; Line(l3)  = {3,4};" << std::endl;
    GmshFile << "l4  = newreg; Line(l4)  = {4,1};" << std::endl;

    GmshFile << "Transfinite Line {l1:l4} = " << nodesOnEdge << ";" << std::endl;

    GmshFile << "loop1 = newll; Line Loop(loop1) = {l1,l2,l3,l4};" << std::endl;

    GmshFile << "plane1 = news; Plane Surface(plane1) = {loop1};" << std::endl;

    GmshFile << "Transfinite Surface {plane1};" << std::endl;
    GmshFile << "Physical Surface(111) = {plane1};" << std::endl;

    GmshFile.close();

    std::string command("gmsh -2 -order 1 " + rGmshFile + ".geo");
    std::cout << command << std::endl;
    std::system(command.c_str());
}

void ImportMesh(NuTo::Structure& rStructure, const int nodesOnEdge, const int rMPIrank)
{
    std::cout << "***********************************" << std::endl;
    std::cout << "**      Import Mesh rank = " << rMPIrank << std::endl;
    std::cout << "***********************************" << std::endl;



    std::string meshFile;
    std::string geoFile;

    switch (rMPIrank)
    {
    case 0:
        geoFile = "/home/phuschke/meshFiles/feti00";
//        CreateGmshFile(geoFile, 0, 0, 20, 10, nodesOnEdge);

        meshFile = "/home/phuschke/meshFiles/feti00.msh";
        break;
    case 1:
        geoFile = "/home/phuschke/meshFiles/feti01";
//        CreateGmshFile(geoFile, 20, 0, 40, 10, nodesOnEdge);

        meshFile = "/home/phuschke/meshFiles/feti01.msh";
        break;
    case 2:
        geoFile = "/home/phuschke/meshFiles/feti02";
//        CreateGmshFile(geoFile, 40, 0, 60, 10, nodesOnEdge);

        meshFile = "/home/phuschke/meshFiles/feti02.msh";
        break;
    default:
        throw NuTo::MechanicsException("Number of tasks > number of subdomains.");
        break;
    }


    NuTo::FullMatrix<int, Eigen::Dynamic, Eigen::Dynamic> createdGroupIdMatrix = rStructure.ImportFromGmsh(meshFile, NuTo::ElementData::CONSTITUTIVELAWIP, NuTo::IpData::eIpDataType::STATICDATA);
    int subdomain00 = createdGroupIdMatrix.GetValue(0, 0);
    int interpolationType00 = createdGroupIdMatrix.GetValue(0, 1);

    rStructure.InterpolationTypeAdd(interpolationType00, NuTo::Node::DISPLACEMENTS, NuTo::Interpolation::EQUIDISTANT2);
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
    int numProcesses;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcesses);

    using Eigen::VectorXd;
    using Eigen::MatrixXd;
    using std::cout;
    using std::endl;

    constexpr int nodesOnEdge = 11;

    NuTo::Structure structure(Parameters::mDimension);
    structure.SetVerboseLevel(10);
    structure.SetShowTime(false);

    ImportMesh(structure, nodesOnEdge, rank);
    AssignMaterial(structure, rank);
    AssignSection(structure, rank);

    structure.NodeBuildGlobalDofs();
    const int numDofs  = structure.GetNumDofs(NuTo::Node::DISPLACEMENTS);
    const int numNodes = structure.GetNumNodes();
    structure.CalculateMaximumIndependentSets();

    Eigen::SparseMatrix<double,Eigen::ColMajor> stiffnessMatrixSparse(numDofs, numDofs);
    AssembleStiffnesMatrix(structure, stiffnessMatrixSparse);
    stiffnessMatrixSparse.makeCompressed();

    NuTo::FullVector<double,Eigen::Dynamic> externalForce(numDofs);
    constexpr double searchTol = 1.0e-6;

    structure.Info();

    // assemble load vector for rank 2
//    if (rank == 2)
//    {
//        constexpr double xCoord = 60.0;
//        AssenbleExternalForceVector(searchTol, externalForce, structure, xCoord);
//    }

    // assemble load vector for rank 1
    if (rank == 1)
    {
        constexpr double xCoord = 60.0;
        AssenbleExternalForceVectorThreePointBending(searchTol, externalForce, structure, xCoord);
    }

    // get all nodes on the interface
//    Eigen::SparseMatrix<double> connectivityMatrix = AssembleConnectivityMatrix(searchTol, rank, structure, nodesOnEdge);
    Eigen::SparseMatrix<double> connectivityMatrix = AssembleConnectivityMatrixThreePointBending(searchTol, rank, structure, nodesOnEdge);

    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();

    std::cout << "***********************************" << std::endl;
    std::cout << "**      FETI rank = " << rank        << std::endl;
    std::cout << "***********************************" << std::endl;

    const int numInterfaceEqs   = connectivityMatrix.rows();

    // compute the total number of rigid body modes
    constexpr int numRigidBodyModesLocal    = 3;
    const int numRigidBodyModesGlobal       = numProcesses * numRigidBodyModesLocal;

    Eigen::MatrixXd rigidBodyModes = Eigen::MatrixXd::Zero(numDofs,numRigidBodyModesLocal);

    AssembleRigidBodyModes(numNodes, structure, rigidBodyModes);

    Eigen::MatrixXd interfaceRigidBodyModesLocal    = connectivityMatrix * rigidBodyModes;
    Eigen::VectorXd rigidBodyForceVectorLocal       = rigidBodyModes.transpose() * externalForce;

    // Robust Cholesky decomposition of a matrix with pivoting.
//    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double,Eigen::ColMajor>> solver;
    Eigen::SparseQR<Eigen::SparseMatrix<double>,Eigen::COLAMDOrdering<int>> solver;
//    Eigen::SparseLU<Eigen::SparseMatrix<double,Eigen::ColMajor>> solver;

    solver.compute(stiffnessMatrixSparse);

    Eigen::VectorXd displacementGapLocal    = connectivityMatrix * solver.solve(externalForce);
    Eigen::VectorXd displacementGapSum      = Eigen::VectorXd::Zero(numInterfaceEqs);

    MPI_Allreduce(displacementGapLocal.data(), displacementGapSum.data(), displacementGapLocal.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    Eigen::MatrixXd interfaceRigidBodyModesGlobal   = Eigen::MatrixXd::Zero(numInterfaceEqs,numRigidBodyModesGlobal);
    MPI_Allgather(interfaceRigidBodyModesLocal.data(), interfaceRigidBodyModesLocal.size(), MPI_DOUBLE, interfaceRigidBodyModesGlobal.data(), interfaceRigidBodyModesLocal.size(), MPI_DOUBLE, MPI_COMM_WORLD);

    Eigen::VectorXd rigidBodyForceVectorGlobal      = Eigen::VectorXd::Zero(numRigidBodyModesGlobal);
    MPI_Allgather(rigidBodyForceVectorLocal.data(), rigidBodyForceVectorLocal.size(), MPI_DOUBLE, rigidBodyForceVectorGlobal.data(), rigidBodyForceVectorLocal.size(), MPI_DOUBLE, MPI_COMM_WORLD);

    Eigen::VectorXd alphaGlobal = Eigen::VectorXd::Zero(numRigidBodyModesGlobal);

    std::cout << "***********************************" << std::endl;
    std::cout << "**   Conjugate Projected Gradient**" << std::endl;
    std::cout << "***********************************" << std::endl;

    // this matrix is very small, i.e. 3 * number of subdomains in 2d. The inversion is therefore not expensive
    Eigen::MatrixXd GtransGinv = (interfaceRigidBodyModesGlobal.transpose() * interfaceRigidBodyModesGlobal).inverse();

    // this projection matrix guarantees that the solution for lambda satisfies the constraint at every iteration, if the following initial start vector for lambda is chosen
    Eigen::MatrixXd projection = Eigen::MatrixXd::Identity(numInterfaceEqs, numInterfaceEqs) - interfaceRigidBodyModesGlobal * GtransGinv * interfaceRigidBodyModesGlobal.transpose();
    Eigen::VectorXd lambda = interfaceRigidBodyModesGlobal * GtransGinv * rigidBodyForceVectorGlobal;

    constexpr double tol    = 1.0e-6;
    const int maxIters      = 1000;

    VectorXd rhs    = displacementGapSum;
    VectorXd x      = lambda;

    // it is cheaper to compute B_s * K_s^+ B_s^T * x than to assemble the F_i matrix directly
    VectorXd tmpLocal = connectivityMatrix * solver.solve(connectivityMatrix.transpose() * x);
    VectorXd tmpSum;
    tmpSum.setZero(tmpLocal.rows());
    MPI_Allreduce(tmpLocal.data(), tmpSum.data(), tmpLocal.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    //initial residual
    VectorXd residual = rhs - tmpSum;

//    double rhsNorm2 = rhs.squaredNorm();
//    double threshold = tol*tol*rhsNorm2;

    //initial projected search direction
    VectorXd w = projection * residual;
    VectorXd p = w;


//    double residualNorm2 = residual.squaredNorm();
    double threshold = tol * w.squaredNorm();


    double absNew = w.dot(w);
    double wNorm2;
    int i = 0;

    while(i < maxIters)
    {
      // at every iteration i the residual has to be recomputd which is quite expensive
      tmpLocal.noalias() = connectivityMatrix * solver.solve(connectivityMatrix.transpose() * p);
      tmpSum.setZero();
      MPI_Allreduce(tmpLocal.data(), tmpSum.data(), tmpLocal.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

      // step size
      double alpha = absNew / p.dot(tmpSum);

      // update solution
      x += alpha * p;

      // update residue
      residual -= alpha * tmpSum;

      // project the residual to guarantee the constraints
      w = projection * residual;

      wNorm2 = w.squaredNorm();
      if (wNorm2 < threshold)
          break;

      double absOld = absNew;
      absNew = w.dot(w);                        // update the absolute value of r
      double beta = absNew / absOld;            // calculate the Gram-Schmidt value used to create the new search direction
      p = w + beta * p;                         // update search direction
      i++;
    }

    cout << "numInterations" << i << endl;
    lambda = x;

    tmpLocal.noalias() = connectivityMatrix * solver.solve(connectivityMatrix.transpose() * x);
    tmpSum.setZero();
    MPI_Allreduce(tmpLocal.data(), tmpSum.data(), tmpLocal.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    alphaGlobal = GtransGinv * interfaceRigidBodyModesGlobal.transpose() * (displacementGapSum - tmpSum);

    Eigen::VectorXd alphaLocal = alphaGlobal.segment(numRigidBodyModesLocal * rank, numRigidBodyModesLocal);

    VectorXd tmp = externalForce - connectivityMatrix.transpose() * lambda;
    Eigen::VectorXd displacement = solver.solve(tmp);

    displacement -= rigidBodyModes * alphaLocal;

    if (rank == 0)
    {
        end = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds = end - start;
        std::cout << "elapsed time: " << elapsed_seconds.count() << "s\n";
    }

    std::cout << "***********************************" << std::endl;
    std::cout << "**      Visualization rank = " << rank        << std::endl;
    std::cout << "***********************************" << std::endl;


    int groupAllElements = structure.GroupCreate(NuTo::Groups::Elements);
    structure.GroupAddElementsTotal(groupAllElements);
    structure.AddVisualizationComponent(groupAllElements, NuTo::VisualizeBase::DISPLACEMENTS);
    structure.AddVisualizationComponent(groupAllElements, NuTo::VisualizeBase::ENGINEERING_STRESS);
    structure.AddVisualizationComponent(groupAllElements, NuTo::VisualizeBase::ENGINEERING_STRAIN);

    char outputFile[200];
    sprintf(outputFile, "/home/phuschke/structure_rank_%03d_0.vtk", rank);
    structure.ExportVtkDataFileElements(outputFile);

    for (int iNode = 0; iNode < numNodes; ++iNode)
    {

        NuTo::FullVector<int, -1> displacementDofs;
        structure.NodeGetDisplacementDofs(iNode, displacementDofs);

        NuTo::FullVector<double, Eigen::Dynamic> nodeDisplacements(Parameters::mDimension);
        nodeDisplacements(0) = displacement(displacementDofs(0));
        nodeDisplacements(1) = displacement(displacementDofs(1));

        structure.NodeSetDisplacements(iNode, nodeDisplacements);

    }

    sprintf(outputFile, "/home/phuschke/structure_rank_%03d_1.vtk", rank);
    structure.ExportVtkDataFileElements(outputFile);

    if (rank == 0)
    {
        cout << "size local: " << stiffnessMatrixSparse.cols() << endl;
        cout << "size global: " << connectivityMatrix.rows() << endl;

    }

    MPI_Finalize();
}
