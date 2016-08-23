
#include <iostream>
#include <string.h>

#include "nuto/mechanics/structures/unstructured/Structure.h"
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/SparseLU>
#include <memory>

using FullMatrixXi      = NuTo::FullMatrix<int, Eigen::Dynamic, Eigen::Dynamic>;
using FullVectorXi      = NuTo::FullVector<int, Eigen::Dynamic>;
using StructureList     = std::vector<std::unique_ptr<NuTo::Structure>>;
using SparseMatrixD     = Eigen::SparseMatrix<double>;
using SparseMatrixList  = std::vector<std::unique_ptr<SparseMatrixD>>;

using MatrixList        = std::vector<std::unique_ptr<Eigen::MatrixXd>>;
using VectorList        = std::vector<std::unique_ptr<Eigen::VectorXd>>;
using Solver            = Eigen::SparseLU<Eigen::SparseMatrix<double,Eigen::ColMajor>>;
using SolverList        = std::vector<std::unique_ptr<Solver>>;
////////////////////////////////
//          MAIN
////////////////////////////////
int main(int argc, char* argv[])
{

//    try
//    {
//        constexpr int dimension             = 2;
//        constexpr int num_subdomains        = 3;
//        constexpr int num_rigid_body_modes  = 3;
//        constexpr double tol                = 1e-8;
//        StructureList       structureList;
//        SparseMatrixList    connectivityMatrixList;
//        SparseMatrixList    stiffnessMatrixList;
//        MatrixList          rigidBodyModesList;
//        MatrixList          interfaceRigidBodyModesList;

//        VectorList          externalForceVectorList;
//        VectorList          rigidBodyExternalForceVectorList;
//        SolverList          solverList;

//        std::set<int>       nodeIDsInterfaceGlobal;
//        std::set<int>       nodeIDsBoundaryGlobal;
//        std::vector<std::vector<int>>    nodeIDsInterfaceMasterList;
//        std::vector<std::vector<int>>    nodeIDsInterfaceSlaveList;
//        std::vector<std::vector<int>>    nodeIDsBoundaryList;
//        std::vector<int> gmshPartition = {2,1,3};

//        for (int i = 0; i<num_subdomains; ++i)
//        {
//            structureList.push_back(std::make_unique<NuTo::Structure>(dimension));
//            connectivityMatrixList.push_back(std::make_unique<SparseMatrixD>());
//            stiffnessMatrixList.push_back(std::make_unique<SparseMatrixD>());

//            rigidBodyModesList.push_back(std::make_unique<Eigen::MatrixXd>());
//            interfaceRigidBodyModesList.push_back(std::make_unique<Eigen::MatrixXd>());
//            externalForceVectorList.push_back(std::make_unique<Eigen::VectorXd>());
//            rigidBodyExternalForceVectorList.push_back(std::make_unique<Eigen::VectorXd>());
//            solverList.push_back(std::make_unique<Solver>());
//        }


//        // assemble interface and boundary node lists
//        for (int i_subdomain = 0; i_subdomain < num_subdomains; ++i_subdomain)
//        {
//            auto& structure = structureList[i_subdomain];
//            structure->SetVerboseLevel(10);
//            structure->SetShowTime(false);
//            char meshFile[200];
//            sprintf(meshFile, "/home/phuschke/meshFiles/2d/feti_beam.msh_00000%d", i_subdomain +1);
//            std::cout << "subdomain" << std::endl << meshFile << std::endl;


//            NuTo::FullMatrix<int, Eigen::Dynamic, Eigen::Dynamic> ids = structure->ImportFromGmsh(meshFile, NuTo::ElementData::CONSTITUTIVELAWIP, NuTo::IpData::eIpDataType::STATICDATA);
////            std::vector<int> physicalGroups = structure->ImportGmshMeshFile(meshFile, NuTo::ElementData::CONSTITUTIVELAWIP, NuTo::IpData::eIpDataType::STATICDATA);
//            std::vector<int> physicalGroups;
//            for (int i = 0; i < ids.rows(); ++i)
//                physicalGroups.push_back(ids.at(i,0));

//            std::cout << "finished importing" << std::endl;

//            int groupNodesBoundary          = structure->GroupCreate(NuTo::Groups::Nodes);
//            int groupNodesInterfaceMaster   = structure->GroupCreate(NuTo::Groups::Nodes);
//            int groupNodesInterfaceSlave    = structure->GroupCreate(NuTo::Groups::Nodes);
//            int groupElementsSubdomain      = 0;

//            for (auto& physicalGroupId : physicalGroups)
//            {
//                if      (physicalGroupId == 1000 + gmshPartition[i_subdomain])
//                {
//                    structure->GroupAddNodesFromElements(groupNodesBoundary, physicalGroupId);
////                    structure->GroupDeleteElements(physicalGroupId);
//                }
//                else if (physicalGroupId == 2000 + gmshPartition[i_subdomain])
//                {
//                    structure->GroupAddNodesFromElements(groupNodesInterfaceMaster, physicalGroupId);
////                    structure->GroupDeleteElements(physicalGroupId);
//                }
//                else if (physicalGroupId == 3000 + gmshPartition[i_subdomain])
//                {
//                    structure->GroupAddNodesFromElements(groupNodesInterfaceSlave, physicalGroupId);
////                    structure->GroupDeleteElements(physicalGroupId);
//                }
//                else if (physicalGroupId == 4000 + gmshPartition[i_subdomain])
//                {
//                    groupElementsSubdomain = physicalGroupId;
//                }
//                else
//                {
//                    std::cout << "Check mesh file" << std::endl;
//                    return -1;
//                }
//            }

//            int sectionID  = structure->SectionCreate(NuTo::Section::PLANE_STRAIN);
//            structure->ElementGroupSetSection(groupElementsSubdomain,sectionID);

//            int materialID = structure->ConstitutiveLawCreate(NuTo::Constitutive::LINEAR_ELASTIC_ENGINEERING_STRESS);
//            structure->ElementGroupSetConstitutiveLaw(groupElementsSubdomain,materialID);

//            int interpolationType00 = structure->InterpolationTypeCreate(NuTo::Interpolation::eShapeType::QUAD2D);
//            structure->InterpolationTypeAdd(interpolationType00, NuTo::Node::COORDINATES, NuTo::Interpolation::EQUIDISTANT1);
//            structure->InterpolationTypeAdd(interpolationType00, NuTo::Node::DISPLACEMENTS, NuTo::Interpolation::EQUIDISTANT1);
//            structure->ElementGroupSetInterpolationType(groupElementsSubdomain, interpolationType00);
//            structure->Info();

//            structure->ElementConvertToInterpolationType(groupElementsSubdomain, 1.e-8, 10);
//            structure->NodeBuildGlobalDofs();





//            Eigen::VectorXi     nodeIDsInterfaceMasterEigen          = structure->GroupGetMemberIds(groupNodesInterfaceMaster);
//            nodeIDsInterfaceMasterList.push_back(std::vector<int>(nodeIDsInterfaceMasterEigen.data(), nodeIDsInterfaceMasterEigen.data() + nodeIDsInterfaceMasterEigen.size()));

//            Eigen::VectorXi     nodeIDsInterfaceSlaveEigen          = structure->GroupGetMemberIds(groupNodesInterfaceMaster);
//            nodeIDsInterfaceSlaveList.push_back(std::vector<int>(nodeIDsInterfaceSlaveEigen.data(), nodeIDsInterfaceSlaveEigen.data() + nodeIDsInterfaceSlaveEigen.size()));

//            Eigen::VectorXi     nodeIDsBoundaryEigen          = structure->GroupGetMemberIds(groupNodesBoundary);
//            nodeIDsBoundaryList.push_back(std::vector<int>(nodeIDsBoundaryEigen.data(), nodeIDsBoundaryEigen.data() + nodeIDsBoundaryEigen.size()));

//            for (const auto& nodeID : nodeIDsInterfaceMasterList[i_subdomain])
//                nodeIDsInterfaceGlobal.insert(nodeID);

//            for (const auto& nodeID : nodeIDsBoundaryList[i_subdomain])
//                nodeIDsBoundaryGlobal.insert(nodeID);

//            structure->AddVisualizationComponent(groupElementsSubdomain, NuTo::VisualizeBase::DISPLACEMENTS);


//        }

//        int subdomainID = 0;
//        for (const auto& subdomain : nodeIDsInterfaceMasterList)
//        {
//            std::cout << "Master interface subdomain " << subdomainID << std::endl;
//            for (const auto& nodeID : subdomain)
//            {
//                std::cout << nodeID << std::endl;
//            }
//            ++subdomainID;
//        }

//        subdomainID = 0;
//        for (const auto& subdomain : nodeIDsBoundaryList)
//        {
//            std::cout << "Boundary subdomain " << subdomainID << std::endl;
//            for (const auto& nodeID : subdomain)
//            {
//                std::cout << nodeID << std::endl;
//            }
//            ++subdomainID;
//        }

//        subdomainID = 0;
//        for (const auto& subdomain : nodeIDsInterfaceSlaveList)
//        {
//            std::cout << "Slave interface subdomain " << subdomainID << std::endl;
//            for (const auto& nodeID : subdomain)
//            {
//                std::cout << nodeID << std::endl;
//            }
//            ++subdomainID;
//        }


//        // assemble connectivity matrix
//        for (int i_subdomain = 0; i_subdomain < num_subdomains; ++i_subdomain)
//        {
//            auto& structure = structureList[i_subdomain];
//            connectivityMatrixList.push_back(std::make_unique<Eigen::SparseMatrix<double>>());

//            auto& connectivityMatrix = connectivityMatrixList[i_subdomain];
//            connectivityMatrix->resize(dimension * nodeIDsInterfaceGlobal.size(), structure->GetNumDofs(NuTo::Node::DISPLACEMENTS));

//            auto& nodeIDsInterfaceMaster = nodeIDsInterfaceMasterList[i_subdomain];
//            auto& nodeIDsInterfaceSlave = nodeIDsInterfaceSlaveList[i_subdomain];

//            int i_node = 0;
//            for (const auto& nodeID : nodeIDsInterfaceGlobal)
//            {
//                if (std::find(nodeIDsInterfaceMaster.begin(), nodeIDsInterfaceMaster.end(), nodeID) != nodeIDsInterfaceMaster.end())
//                {
//                    FullVectorXi displacementDofs;
//                    structure->NodeGetDisplacementDofs(nodeID, displacementDofs);
//                    int index = 2 * i_node;
//                    connectivityMatrix->insert(index     , displacementDofs[0]) = 1;
//                    connectivityMatrix->insert(index+1   , displacementDofs[1]) = 1;
//                }
//                else if (std::find(nodeIDsInterfaceSlave.begin(), nodeIDsInterfaceSlave.end(), nodeID) != nodeIDsInterfaceSlave.end())
//                {
//                    FullVectorXi displacementDofs;
//                    structure->NodeGetDisplacementDofs(nodeID, displacementDofs);
//                    int index = 2 * i_node;
//                    connectivityMatrix->insert(index     , displacementDofs[0]) = -1;
//                    connectivityMatrix->insert(index+1   , displacementDofs[1]) = -1;
//                }

//                ++i_node;
//            }


//        }

//        // assemble stiffness matrix
//        for (int i_subdomain = 0; i_subdomain < num_subdomains; ++i_subdomain)
//        {
//            auto& structure = structureList[i_subdomain];
//            NuTo::StructureOutputBlockMatrix stiffnessMatrix = structure->BuildGlobalHessian0();

//            NuTo::SparseMatrixCSRGeneral<double> stiffnessMatrixCSR(stiffnessMatrix.JJ(NuTo::Node::DISPLACEMENTS, NuTo::Node::DISPLACEMENTS));
//            std::vector<Eigen::Triplet<double>> tripletList;


//            std::vector<double> val = stiffnessMatrixCSR.GetValues();
//            std::vector<int> colInd = stiffnessMatrixCSR.GetColumns();
//            std::vector<int> rowInd = stiffnessMatrixCSR.GetRowIndex();

//            for (unsigned i = 0; i < rowInd.size() - 1; ++i)
//            {
//                for (int k = rowInd[i]; k < rowInd[i + 1]; ++k)
//                    tripletList.push_back(Eigen::Triplet<double>(i, colInd[k], val[k]));
//            }


//            stiffnessMatrixList[i_subdomain]->setFromTriplets(tripletList.begin(), tripletList.end());
//            std::cout << "stiffnessMatrixList[i_subdomain]->rows(); " << stiffnessMatrixList[i_subdomain]->rows() << std::endl;
//            std::cout << "stiffnessMatrixList[i_subdomain]->cols(); " << stiffnessMatrixList[i_subdomain]->cols() << std::endl;

//        }

//        // assemble rigid body modes
//        for (int i_subdomain = 0; i_subdomain < num_subdomains; ++i_subdomain)
//        {
//            auto& structure                 = structureList[i_subdomain];
//            auto& rigidBodyModes            = *rigidBodyModesList[i_subdomain];
//            auto& interfaceRigidBodyModes   = *interfaceRigidBodyModesList[i_subdomain];
//            auto& connectivityMatrix        = *connectivityMatrixList[i_subdomain];

//            rigidBodyModes.resize(structure->GetNumDofs(NuTo::Node::DISPLACEMENTS), num_rigid_body_modes);

//            for (const auto& node : structure->NodeGetNodeMap())
//            {
//                const int          nodeId      = node->first;
//                auto const&        nodePtr     = node->second;
//                std::vector<int>   nodeDofs    = structure->NodeGetDisplacementDofs(nodeId);

//                const Eigen::Vector2d coords = nodePtr->Get(NuTo::Node::COORDINATES);
//                rigidBodyModes.block(nodeDofs[0], 0, 1, 3) << 1., 0., -coords[1];
//                rigidBodyModes.block(nodeDofs[1], 0, 1, 3) << 0., 1.,  coords[0];
//            }

//            interfaceRigidBodyModes = connectivityMatrix * rigidBodyModes;

//        }

//        // assemble external force vector
//        for (int i_subdomain = 0; i_subdomain < num_subdomains; ++i_subdomain)
//        {
//            auto& structure                 = structureList[i_subdomain];
//            auto& externalForceVector       = *externalForceVectorList[i_subdomain];
//            externalForceVector.resize(structure->GetNumDofs(NuTo::Node::DISPLACEMENTS));
//            auto& rigidBodyExternalForceVector = *rigidBodyExternalForceVectorList[i_subdomain];
//            auto& rigidBodyModes            = *rigidBodyModesList[i_subdomain];


//            int loadNodeGroup = structure->GroupCreate(NuTo::Groups::Nodes);
//            structure->GroupAddNodeCoordinateRange(loadNodeGroup, 0, 60 - tol, 60 + tol);
//            auto loadNodes = structure->GroupGetMemberIds(loadNodeGroup);
//            for (int iNode = 0; iNode < loadNodes.rows(); ++iNode)
//            {
//                std::vector<int>   nodeDofs    = structure->NodeGetDisplacementDofs(loadNodes[iNode]);
//                externalForceVector(nodeDofs[0]) = 0;
//                externalForceVector(nodeDofs[1]) = 0.1;
//            }
//            rigidBodyExternalForceVector = rigidBodyModes.transpose() * externalForceVector;


//            structure->Info();


//            char outputFile[200];
//            sprintf(outputFile, "result_00000%d.vtu", i_subdomain);
//            structure->ExportVtkDataFileElements(outputFile,true);
//        }




//        Eigen::VectorXd displacementGapSum = Eigen::VectorXd::Zero(dimension * nodeIDsInterfaceGlobal.size());
//        // assemble PCG
//        for (int i_subdomain = 0; i_subdomain < num_subdomains; ++i_subdomain)
//        {
//            auto& structure                 = structureList[i_subdomain];
//            auto& solver                    = *solverList[i_subdomain];
//            auto& stiffnessMatrix           = *stiffnessMatrixList[i_subdomain];
//            auto& externalForceVector       = *externalForceVectorList[i_subdomain];
//            auto& connectivityMatrix        = *connectivityMatrixList[i_subdomain];
//            std::cout << "subdomain " << i_subdomain << std::endl;
//            std::cout << "stiffnessMatrix.rows(); " << stiffnessMatrix.rows() << std::endl;
//            std::cout << "stiffnessMatrix.cols(); " << stiffnessMatrix.cols() << std::endl;

//            solver.compute(stiffnessMatrix);
//            std::cout << "subdomain " << i_subdomain << std::endl;
//            displacementGapSum += connectivityMatrix * solver.solve(externalForceVector);


//        }


//    }
//    catch (NuTo::MechanicsException& e)
//    {
//        std::cout << e.ErrorMessage();
//        return EXIT_FAILURE;
//    }
//    catch(...)
//    { std::cout << "Damn it" << std::endl;}


}



























