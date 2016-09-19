
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
#include "ImportMesh.h"
#include "FetiSolver.h"

#include "../myNutoExamples/EnumsAndTypedefs.h"

#include "nuto/mechanics/nodes/NodeBase.h"

#include "nuto/math/SparseMatrixCSR.h"
#include "nuto/math/SparseMatrixCSRGeneral.h"
#include "nuto/mechanics/dofSubMatrixStorage/BlockSparseMatrix.h"
#include "nuto/mechanics/structures/StructureOutputBlockMatrix.h"


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

    int section00 = rStructure.SectionCreate(NuTo::eSectionType::PLANE_STRESS);
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
    NuTo::SparseMatrixCSRGeneral<double> stiffnessMatrixCSR(stiffnessMatrix.JJ(NuTo::Node::eDof::DISPLACEMENTS, NuTo::Node::eDof::DISPLACEMENTS));

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

void AssembleExternalForceVectorThreePointBending(const double searchTol, NuTo::FullVector<double, Eigen::Dynamic>& externalForce, NuTo::Structure& structure, double xCoord)
{
    NuTo::FullVector<double, 2> nodeCoords;
    nodeCoords[0] = 30;
    nodeCoords[1] = 0;
    int loadNodeGroup = structure.GroupCreate(NuTo::eGroupId::Nodes);
    structure.GroupAddNodeRadiusRange(loadNodeGroup, nodeCoords, 0, 1.0e+1);

    auto loadNodes = structure.GroupGetMemberIds(loadNodeGroup);
    for (int iNode = 0; iNode < loadNodes.rows(); ++iNode)
    {
        NuTo::FullVector<int, -1> displacementDofs;
        structure.NodeGetDisplacementDofs(loadNodes.at(iNode, 0), displacementDofs);
        externalForce(displacementDofs.at(0, 0)) = 0;
        externalForce(displacementDofs.at(1, 0)) = -1;
    }
}

void AssembleRigidBodyModes(const int numNodes, NuTo::Structure& structure, Eigen::MatrixXd& rigidBodyModes)
{

    const auto& nodeMap = structure.NodeGetNodeMap();

    for (const auto& node : nodeMap)
    {
        NuTo::FullVector<int, Eigen::Dynamic> displacementDofs;
        structure.NodeGetDisplacementDofs(node.first, displacementDofs);
        const Eigen::Matrix<double, 2, 1> coordinates = node.second->Get(NuTo::Node::eDof::COORDINATES);
        rigidBodyModes.block(displacementDofs[0], 0, 1, 3) << 1., 0., -coordinates.at(1, 0);
        rigidBodyModes.block(displacementDofs[1], 0, 1, 3) << 0., 1.,  coordinates.at(0, 0);
    }



}

ImportContainer ImportMesh(int rank)
{
    std::string meshFile;

    switch (rank)
    {
    case 0:
        meshFile = "/home/phuschke/nuto_project/nuto/myNutoExamples/meshFiles/2d/feti/fine/feti.msh_000001";
        break;
    case 1:
        meshFile = "/home/phuschke/nuto_project/nuto/myNutoExamples/meshFiles/2d/feti/fine/feti.msh_000002";
        break;
    case 2:
        meshFile = "/home/phuschke/nuto_project/nuto/myNutoExamples/meshFiles/2d/feti/fine/feti.msh_000003";
        break;
    default:
        throw NuTo::MechanicsException("Number of tasks > number of subdomains.");
        break;
    }
    ImportContainer importContainer = ImportMeshFile(meshFile);

    return importContainer;
}

int main(int argc, char* argv[])
{
    MPI_Init(&argc, &argv);
    int rank            = MPI::COMM_WORLD.Get_rank();
    int numProcesses    = MPI::COMM_WORLD.Get_size();

    try
    {
        using Eigen::VectorXd;
        using Eigen::MatrixXd;
        using std::cout;
        using std::endl;

        constexpr int dim = 2;
        NuTo::Structure structure(dim);
        structure.SetVerboseLevel(10);
        structure.SetShowTime(false);

        ImportContainer importContainer = ImportMesh(rank);

        for (const auto& node : importContainer.mNodeList)
            structure.NodeCreate(node.mId, node.mCoordinates.head(dim));


        int interpolationType = structure.InterpolationTypeCreate(NuTo::Interpolation::eShapeType::QUAD2D);
        structure.InterpolationTypeAdd(interpolationType, NuTo::Node::eDof::COORDINATES, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);
        structure.InterpolationTypeAdd(interpolationType, NuTo::Node::eDof::DISPLACEMENTS, NuTo::Interpolation::eTypeOrder::EQUIDISTANT1);


        for (const auto& element : importContainer.mElementList)
            structure.ElementCreate(element.mId,interpolationType, element.mNodeIds,NuTo::ElementData::eElementDataType::CONSTITUTIVELAWIP,NuTo::IpData::eIpDataType::STATICDATA);


        structure.ElementTotalConvertToInterpolationType();

        AssignMaterial(structure, rank);
        AssignSection(structure, rank);

        structure.NodeBuildGlobalDofs();

        const int numDofs  = structure.GetNumDofs(NuTo::Node::eDof::DISPLACEMENTS);
        const int numNodes = structure.GetNumNodes();
        structure.CalculateMaximumIndependentSets();

        Eigen::SparseMatrix<double,Eigen::ColMajor> stiffnessMatrixSparse(numDofs, numDofs);
        AssembleStiffnesMatrix(structure, stiffnessMatrixSparse);
        stiffnessMatrixSparse.makeCompressed();

        NuTo::FullVector<double,Eigen::Dynamic> externalForce;
        externalForce.setZero(numDofs);
        constexpr double searchTol = 1.0e-6;

        constexpr double xCoord = 60.0;
        AssembleExternalForceVectorThreePointBending(searchTol, externalForce, structure, xCoord);


        const int num_interface_nodes_global = 42;
        const int num_boundary_nodes_global = 42;
        const int num_lagrange_multipliers  = dim * (num_interface_nodes_global + num_boundary_nodes_global);
        Eigen::SparseMatrix<double> connectivityMatrix(num_lagrange_multipliers, numDofs);

        for (const auto& interface : importContainer.mInterfaceList)
            for (const auto& nodeId : interface.mNodeIdsMap)
            {
                int globalIndex = dim * nodeId.first;
                NuTo::FullVector<int, Eigen::Dynamic> displacementDofs;
                structure.NodeGetDisplacementDofs(nodeId.second, displacementDofs);
                connectivityMatrix.insert(globalIndex   , displacementDofs[0]) = interface.mValue;
                connectivityMatrix.insert(globalIndex +1, displacementDofs[1]) = interface.mValue;
            }

        for (const auto& boundary : importContainer.mBoundaryList)
            for (const auto& nodeId : boundary.mNodeIdsMap)
            {
                int globalIndex = dim * (nodeId.first + num_interface_nodes_global);
                NuTo::FullVector<int, Eigen::Dynamic> displacementDofs;
                structure.NodeGetDisplacementDofs(nodeId.second, displacementDofs);
                connectivityMatrix.insert(globalIndex   , displacementDofs[0]) = 1.0;
                connectivityMatrix.insert(globalIndex +1, displacementDofs[1]) = 1.0;
            }





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


        //    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double,Eigen::ColMajor>> solver;
        //    Eigen::SparseLU<Eigen::SparseMatrix<double,Eigen::ColMajor>> solver;
        Eigen::SparseQR<Eigen::SparseMatrix<double>,Eigen::COLAMDOrdering<int>> solver;
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
        const int maxIters      = 100;

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

//        auto displacement = solver.solve(externalForce);

        int groupAllElements = structure.GroupCreate(NuTo::eGroupId::Elements);
        structure.GroupAddElementsTotal(groupAllElements);
        structure.AddVisualizationComponent(groupAllElements, NuTo::eVisualizeWhat::DISPLACEMENTS);
        structure.AddVisualizationComponent(groupAllElements, NuTo::eVisualizeWhat::ENGINEERING_STRESS);
        structure.AddVisualizationComponent(groupAllElements, NuTo::eVisualizeWhat::ENGINEERING_STRAIN);

        char outputFile[200];
        sprintf(outputFile, "/home/phuschke/structure_rank_%03d_0.vtk", rank);
        structure.ExportVtkDataFileElements(outputFile);


        const auto& nodeMap = structure.NodeGetNodeMap();

        for (const auto& node : nodeMap)
        {
            NuTo::FullVector<int, Eigen::Dynamic> displacementDofs;
            structure.NodeGetDisplacementDofs(node.first, displacementDofs);

            NuTo::FullVector<double, Eigen::Dynamic> nodeDisplacements(Parameters::mDimension);
            nodeDisplacements(0) = displacement(displacementDofs(0));
            nodeDisplacements(1) = displacement(displacementDofs(1));

            structure.NodeSetDisplacements(node.first, nodeDisplacements);

        }

        sprintf(outputFile, "/home/phuschke/structure_rank_%03d_1.vtk", rank);
        structure.ExportVtkDataFileElements(outputFile);

        if (rank == 0)
        {
            cout << "size local: " << stiffnessMatrixSparse.cols() << endl;
            cout << "size global: " << connectivityMatrix.rows() << endl;

        }
    }
    catch(std::exception e)
    {
        std::cout << e.what() << std::endl;
    }
    catch(...)
    {
        std::cout << "Something went wrong" << std::endl;
    }

    MPI_Finalize();
}
