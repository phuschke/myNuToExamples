#include "nuto/mechanics/timeIntegration/TimeIntegrationBase.h"
#include "nuto/mechanics/timeIntegration/NewmarkDirect.h"
#include "nuto/base/ErrorEnum.h"
#include "nuto/mechanics/structures/StructureBase.h"
#include "nuto/mechanics/MechanicsException.h"
#include "nuto/base/Timer.h"

#include "nuto/mechanics/structures/StructureBase.h"
#include "nuto/mechanics/structures/unstructured/StructureFETI.h"
#include "nuto/mechanics/structures/StructureBaseEnum.h"
#include "nuto/mechanics/structures/StructureOutputBlockMatrix.h"
#include "nuto/mechanics/structures/StructureOutputDummy.h"
#include <mpi.h>

#include "nuto/base/CallbackInterface.h"
#include "nuto/math/SparseMatrixCSRGeneral.h"
#include <eigen3/Eigen/Dense>
#include "nuto/mechanics/constitutive/ConstitutiveEnum.h"
#include "nuto/mechanics/nodes/NodeEnum.h"
#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveCalculateStaticData.h"
#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveIOMap.h"
#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveTimeStep.h"
#include "nuto/mechanics/constitutive/staticData/ConstitutiveStaticDataGradientDamage.h"
#include <cmath>

namespace NuTo
{
class FetiSolver : public NewmarkDirect
{
public:

    using VectorXd      = Eigen::VectorXd;
    using MatrixXd      = Eigen::MatrixXd;
    using SparseMatrix  = Eigen::SparseMatrix<double>;

    FetiSolver(StructureBase* rStructure) : NewmarkDirect (rStructure)
    {

    }
    std::string GetTypeId() const override
    {
        throw MechanicsException(__PRETTY_FUNCTION__, "Not implemented!");
    }

    bool HasCriticalTimeStep() const override
    {
        throw MechanicsException(__PRETTY_FUNCTION__, "Not implemented!");
    }

    double CalculateCriticalTimeStep() const override
    {
        throw MechanicsException(__PRETTY_FUNCTION__, "Not implemented!");
    }

//    Eigen::SparseMatrix<double> ConvertToEigenSparseMatrix(NuTo::SparseMatrixCSRVector2<double>& sparseMatrix)
//    {
//        NuTo::SparseMatrixCSRGeneral<double> stiffnessMatrixCSR(sparseMatrix);
//        std::vector<Eigen::Triplet<double>> tripletList;

//        std::vector<double> val = stiffnessMatrixCSR.GetValues();
//        std::vector<int> colInd = stiffnessMatrixCSR.GetColumns();
//        std::vector<int> rowInd = stiffnessMatrixCSR.GetRowIndex();

//        for (unsigned i = 0; i < rowInd.size() - 1; ++i)
//        {
//            for (int k = rowInd[i]; k < rowInd[i + 1]; ++k)
//                tripletList.push_back(Eigen::Triplet<double>(i, colInd[k], val[k]));
//        }

//        Eigen::SparseMatrix<double> stiffnessMatrixSparse(sparseMatrix.GetNumRows(), sparseMatrix.GetNumColumns());
//        stiffnessMatrixSparse.setFromTriplets(tripletList.begin(), tripletList.end());
//        stiffnessMatrixSparse.makeCompressed();

//        return stiffnessMatrixSparse;
//    }

//    Eigen::SparseMatrix<double> ConvertToEigenSparseMatrix(StructureOutputBlockMatrix rMatrix)
//    {
//        std::vector<Eigen::Triplet<double>> tripletList;
//        std::vector<double> val;
//        std::vector<int> colInd;
//        std::vector<int> rowInd;

//        // JJ starts at row = 0, col = 0
//        NuTo::SparseMatrixCSRGeneral<double> stiffnessMatrixCSRJJ(rMatrix.JJ(Node::eDof::DISPLACEMENTS, Node::eDof::DISPLACEMENTS));
//        int numColsJJ = stiffnessMatrixCSRJJ.GetNumColumns();
//        int numRowsJJ = stiffnessMatrixCSRJJ.GetNumRows();

//        val = stiffnessMatrixCSRJJ.GetValues();
//        colInd = stiffnessMatrixCSRJJ.GetColumns();
//        rowInd = stiffnessMatrixCSRJJ.GetRowIndex();

//        for (unsigned i = 0; i < rowInd.size() - 1; ++i)
//        {
//            for (int k = rowInd[i]; k < rowInd[i + 1]; ++k)
//                tripletList.push_back(Eigen::Triplet<double>(i, colInd[k], val[k]));
//        }

//        // JK
//        NuTo::SparseMatrixCSRGeneral<double> stiffnessMatrixCSRJK(rMatrix.JK(Node::eDof::DISPLACEMENTS, Node::eDof::DISPLACEMENTS));

//        val = stiffnessMatrixCSRJK.GetValues();
//        colInd = stiffnessMatrixCSRJK.GetColumns();
//        rowInd = stiffnessMatrixCSRJK.GetRowIndex();

//        for (unsigned i = 0; i < rowInd.size() - 1; ++i)
//        {
//            for (int k = rowInd[i]; k < rowInd[i + 1]; ++k)
//                tripletList.push_back(Eigen::Triplet<double>(i, colInd[k] + numColsJJ, val[k]));
//        }

//        // KJ
//        NuTo::SparseMatrixCSRGeneral<double> stiffnessMatrixCSRKJ(rMatrix.KJ(Node::eDof::DISPLACEMENTS, Node::eDof::DISPLACEMENTS));

//        val = stiffnessMatrixCSRKJ.GetValues();
//        colInd = stiffnessMatrixCSRKJ.GetColumns();
//        rowInd = stiffnessMatrixCSRKJ.GetRowIndex();

//        for (unsigned i = 0; i < rowInd.size() - 1; ++i)
//        {
//            for (int k = rowInd[i]; k < rowInd[i + 1]; ++k)
//                tripletList.push_back(Eigen::Triplet<double>(i + numRowsJJ, colInd[k], val[k]));
//        }

//        // KK
//        NuTo::SparseMatrixCSRGeneral<double> stiffnessMatrixCSRKK(rMatrix.KK(Node::eDof::DISPLACEMENTS, Node::eDof::DISPLACEMENTS));
//        val = stiffnessMatrixCSRKK.GetValues();
//        colInd = stiffnessMatrixCSRKK.GetColumns();
//        rowInd = stiffnessMatrixCSRKK.GetRowIndex();

//        for (unsigned i = 0; i < rowInd.size() - 1; ++i)
//        {
//            for (int k = rowInd[i]; k < rowInd[i + 1]; ++k)
//                tripletList.push_back(Eigen::Triplet<double>(i + numRowsJJ, colInd[k] + numRowsJJ, val[k]));
//        }


//        int numDofs = mStructure->GetNumDofs(Node::eDof::DISPLACEMENTS);
//        Eigen::SparseMatrix<double> stiffnessMatrixSparse(numDofs,numDofs);
//        stiffnessMatrixSparse.setFromTriplets(tripletList.begin(), tripletList.end());
//        stiffnessMatrixSparse.makeCompressed();

//        return stiffnessMatrixSparse;
//    }

    MatrixXd GatherInterfaceRigidBodyModes(const Eigen::MatrixXd& interfaceRigidBodyModes, const int numRigidBodyModesGlobal)
    {

        std::vector<int> recvCount;
        std::vector<int> displ;
        MpiGatherRecvCountAndDispls(recvCount, displ, interfaceRigidBodyModes.size());

        const int numInterfaceEqs               = interfaceRigidBodyModes.rows();
        MatrixXd interfaceRigidBodyModesGlobal   = MatrixXd::Zero(numInterfaceEqs,numRigidBodyModesGlobal);
        MPI_Allgatherv(interfaceRigidBodyModes.data(),
                       interfaceRigidBodyModes.size(),
                       MPI_DOUBLE,
                       interfaceRigidBodyModesGlobal.data(),
                       recvCount.data(),
                       displ.data(),
                       MPI_DOUBLE,
                       MPI_COMM_WORLD);

        return interfaceRigidBodyModesGlobal;
    }

    VectorXd GatherRigidBodyForceVector(const Eigen::VectorXd& rigidBodyForceVectorLocal, const int numRigidBodyModesGlobal)
    {

        const int numRigidBodyModesLocal        = rigidBodyForceVectorLocal.rows();
        std::vector<int> recvCount;
        std::vector<int> displ;
        MpiGatherRecvCountAndDispls(recvCount, displ, numRigidBodyModesLocal);

        VectorXd rigidBodyForceVectorGlobal      =VectorXd::Zero(numRigidBodyModesGlobal);
        MPI_Allgatherv(rigidBodyForceVectorLocal.data(),
                       rigidBodyForceVectorLocal.size(),
                       MPI_DOUBLE,
                       rigidBodyForceVectorGlobal.data(),
                       recvCount.data(),
                       displ.data(),
                       MPI_DOUBLE,
                       MPI_COMM_WORLD);

        return rigidBodyForceVectorGlobal;
    }

    void MpiGatherRecvCountAndDispls(std::vector<int>& recvCount, std::vector<int>& displs, const int numValues)
    {
        const int numProcesses = MPI::COMM_WORLD.Get_size();
        // recvCount:
        // Contais the number of elements that are received from each process.
        recvCount.clear();
        recvCount.resize(numProcesses, 0);


        MPI_Allgather(  &numValues,
                        1,
                        MPI_INT,
                        recvCount.data(),
                        1,
                        MPI_INT,
                        MPI_COMM_WORLD );

        // displs:
        // Entry i specifies the displacement (relative to recvbuf) at which to place the incoming data from process i.
        displs.clear();
        displs.resize(numProcesses, 0);
        for (int i = 1; i < numProcesses; ++i)
            displs[i] = displs[i-1] + recvCount[i-1];


    }

    BlockScalar CalculateNormResidual(BlockFullVector<double>& residual_mod, const std::set<Node::eDof>& activeDofSet)
    {
        StructureFETI* structure = static_cast<StructureFETI*>(mStructure);
        for (const auto& interface : structure->mInterfaces)
            for (const auto& nodePair: interface.mNodeIdsMap)
            {
                for (const auto& dof : activeDofSet)
                {
                    const std::vector<int> dofIds = structure->NodeGetDofIds(nodePair.second, dof);
                    for (const auto& id : dofIds)
                        residual_mod[dof][id] = 0;

                }

            }

        BlockScalar normResidual = residual_mod.CalculateInfNorm();

        for (const auto& dof : activeDofSet)
            MPI_Allreduce(MPI_IN_PLACE,  &normResidual[dof], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        return normResidual;

    }



    int BiCgStab(const MatrixXd& projection, VectorXd& x, const VectorXd& rhs)
    {

        StructureFETI* structure = static_cast<StructureFETI*>(mStructure);
        const SparseMatrix B      = structure->GetConnectivityMatrix();
        const SparseMatrix Btrans = B.transpose();


        // solve the linear system Ax = b
        VectorXd Ax = B * mSolver.solve(Btrans*x);
        MPI_Allreduce(MPI_IN_PLACE, Ax.data(), Ax.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        VectorXd& Ay = Ax;
        VectorXd& Az = Ax;

        const int n = rhs.rows();
        VectorXd r  = rhs - Ax;
        VectorXd r0 = r;

        VectorXd rProj   = projection * r;
        VectorXd rProj0  = rProj;

        double rProj0_sqnorm = rProj0.squaredNorm();
        double rhs_sqnorm = rhs.squaredNorm();

        double rho    = 1.;
        double alpha  = 1.;
        double w      = 1.;

        VectorXd v = VectorXd::Zero(n);
        VectorXd p = VectorXd::Zero(n);
        VectorXd y(n);
        VectorXd z(n);
        VectorXd s(n);
        VectorXd t(n);
        VectorXd tProj(n);
        VectorXd vProj(n);

        const double tol    = mCpgTolerance;
        const double tol2   = tol*tol*rhs_sqnorm;

        //const double eps    = std::numeric_limits<double>::epsilon();
        //const double eps2   = eps*eps;

        int i = 0;
        int restarts = 0;



        const int maxIters = mCpgMaxIterations;


        while (rProj.squaredNorm() > tol2 and i<maxIters )
        {
            double rho_old = rho;

            rho = rProj0.dot(rProj);



            if (std::fabs(rho) < 1.0e-20)
            {

                if (structure->mRank == 0)
                {
                    std::cout << "Shit is happending! (rho) " << (rho) <<  std::endl << std::endl;
                    std::cout << "Shit is happending! std::fabs(rho) " << std::fabs(rho) <<  std::endl << std::endl;
                    std::cout << "Shit is happending! (tol) " << tol <<  std::endl << std::endl;
                }

                // The new residual vector became too orthogonal to the arbitrarily chosen direction r0
                // Let's restart with a new r0:

                Ax.noalias() = B * mSolver.solve(Btrans*x);
                MPI_Allreduce(MPI_IN_PLACE, Ax.data(), Ax.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

                rProj  = projection*(rhs - Ax);
                rProj0 = rProj;
                rho = rProj0_sqnorm = rProj.squaredNorm();
                if(restarts++ == 0)
                    i = 0;
            }


            double beta = (rho/rho_old) * (alpha / w);

            p = rProj + beta * (p - w * v);


            //      y = precond.solve(p);
            y = p;

            Ay.noalias() = B * mSolver.solve(Btrans*y);
            MPI_Allreduce(MPI_IN_PLACE, Ay.data(), Ay.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

            v.noalias() = projection * Ay;

            alpha = rho / rProj0.dot(v);

            s = rProj - alpha * v;

            if (s.norm() < 1.0e-15)
            {
                x += alpha *p;

                if (structure->mRank == 0)
                    std::cout << "Early exit!"<< std::endl << std::endl;

                break;
            }

            //      z = precond.solve(s);
            z = s;

            Az.noalias() = B * mSolver.solve(Btrans*z);
            MPI_Allreduce(MPI_IN_PLACE, Az.data(), Az.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

            t.noalias() = projection * Az;

            double tmp = t.dot(t);

            if(tmp>0.0)
                w = t.dot(s) / tmp;
            else
                w = 0.0;

            x += alpha * y + w * z;
            rProj = z - w * t;

            if (structure->mRank == 0)
                std::cout << "BiCGStab rel. error = " << rProj.squaredNorm()/rhs_sqnorm << "\t at iteration = " << i << "/" << maxIters << std::endl;

            ++i;



        }
        //    tol_error = sqrt(r.squaredNorm()/rhs_sqnorm);
        return i;

    }


    int CPG(const MatrixXd& projection, VectorXd& x, const VectorXd& rhs)
    {

        VectorXd        Ap;

        StructureFETI* structure = static_cast<StructureFETI*>(mStructure);
        const SparseMatrix B      = structure->GetConnectivityMatrix();
        const SparseMatrix Btrans = B.transpose();



        VectorXd Ax = B * mSolver.solve(Btrans*x);
        MPI_Allreduce(MPI_IN_PLACE, Ax.data(), Ax.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        //initial residual
        VectorXd r = rhs - Ax;

        //initial projected search direction
        VectorXd w = projection * r;
        VectorXd p = w;

        //    double r = r.squaredNorm();
        double threshold = mCpgTolerance * w.squaredNorm();


        double wTw = w.dot(w);
        int iteration = 0;
        while(iteration < mCpgMaxIterations)
        {
            // at every iteration i the r has to be recomputd which is quite expensive
            Ap.noalias() = B * mSolver.solve(Btrans*p);
            MPI_Allreduce(MPI_IN_PLACE, Ap.data(), Ap.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

            // step size
            double alpha = wTw / p.dot(Ap);

            // update solution
            x += alpha * p;

            // update residual
            r -= alpha * Ap;

            // project the r to guarantee the constraints
            w = projection * r;

            if (w.squaredNorm() < threshold)
                break;

            double wTwOld = wTw;
            wTw = w.dot(w);                        // update the absolute value of r
            double beta = wTw / wTwOld;            // calculate the Gram-Schmidt value used to create the new search direction
            p = w + beta * p;                         // update search direction
            iteration++;
        }

        return iteration;
    }





    //! @brief Conjugate projected gradient method
    StructureOutputBlockVector ConjugateProjectedGradientMethod(BlockFullVector<double> residual_mod, const std::set<Node::eDof>& activeDofSet,VectorXd& lambda)
    {

        StructureFETI* structure = static_cast<StructureFETI*>(mStructure);

        const SparseMatrix& B                    = structure->GetConnectivityMatrix();
        const SparseMatrix& Btrans               = B.transpose();
        const int& numLagrangeMultipliers       = B.rows();

        structure->mNumRigidBodyModes           = mSolver.cols() -  mSolver.rank();
        const int& numRigidBodyModesLocal       = structure->mNumRigidBodyModes;
        int numRigidBodyModesGlobal             = 0;
        MPI_Allreduce(&numRigidBodyModesLocal,
                      &numRigidBodyModesGlobal,
                      1,
                      MPI_INT,
                      MPI_SUM,
                      MPI_COMM_WORLD);




        MPI_Barrier(MPI_COMM_WORLD);

        structure->mInterfaceRigidBodyModes.resize(numLagrangeMultipliers, numRigidBodyModesLocal);

        std::cout << "Rank: " << structure->mRank << " => number of rigid body modes " << numRigidBodyModesLocal << std::endl;
        MPI_Barrier(MPI_COMM_WORLD);

        if (structure->mNumRigidBodyModes > 0)
        {
            // extract the last columns of the Q matrix which represents the null space of K
            Eigen::MatrixXd rbHelper;
            rbHelper.setZero(mSolver.rows(), numRigidBodyModesLocal);
            rbHelper.bottomLeftCorner(numRigidBodyModesLocal,numRigidBodyModesLocal) = Eigen::MatrixXd::Identity(numRigidBodyModesLocal, numRigidBodyModesLocal);
            structure->mRigidBodyModes = mSolver.matrixQ() * rbHelper;
            std::cout << "Rank: " << structure->mRank << " => K*R.maxCoeff() " << mSolver.solve(structure->mRigidBodyModes).maxCoeff() << std::endl;

            structure->mInterfaceRigidBodyModes = B * structure->mRigidBodyModes;

        }



        const auto& rigidBodyModes              = structure->GetRigidBodyModes();
        const auto& interfaceRigidBodyModes     = structure->GetInterfaceRigidBodyModes();

        MPI_Barrier(MPI_COMM_WORLD);


        // G.size() = (number of Lagrange multipliers) x (total number of rigid body modes)
        MatrixXd G          =   GatherInterfaceRigidBodyModes(interfaceRigidBodyModes, numRigidBodyModesGlobal);

        // GtransGinv.size() = (total number of rigid body modes) x (total number of rigid body modes)
        MatrixXd GtransGinv =   (G.transpose() * G).inverse();

        // The projection matrix guarantees that the solution for lambda satisfies the constraint at every iteration
        MatrixXd projection =   MatrixXd::Identity(G.rows(), G.rows()) - G * GtransGinv * G.transpose();

        VectorXd rigidBodyForceVectorLocal(numRigidBodyModesLocal);
        if (structure->mNumRigidBodyModes > 0)
            rigidBodyForceVectorLocal = rigidBodyModes.transpose() * residual_mod.Export();

        VectorXd rigidBodyForceVectorGlobal = GatherRigidBodyForceVector(rigidBodyForceVectorLocal, numRigidBodyModesGlobal);

        MPI_Barrier(MPI_COMM_WORLD);

        VectorXd displacementGap    = B * mSolver.solve(residual_mod.Export());
        MPI_Allreduce(MPI_IN_PLACE,
                      displacementGap.data(),
                      displacementGap.size(),
                      MPI_DOUBLE,
                      MPI_SUM,
                      MPI_COMM_WORLD);

        lambda = G * GtransGinv * rigidBodyForceVectorGlobal;

        VectorXd rhs    = displacementGap;
        VectorXd x      = lambda;

        VectorXd alphaGlobal =VectorXd::Zero(numRigidBodyModesGlobal);


        if (structure->mRank == 0)
        {
            structure->GetLogger() << "\n*************************************\n";
            structure->GetLogger() << "       Start Interface Problem           ";
            structure->GetLogger() << "\n*************************************\n\n";
        }


        MPI_Barrier(MPI_COMM_WORLD);

//        int iteration = CPG(projection,x,rhs);
        int iteration = BiCgStab(projection,x,rhs);

        if (iteration >= mCpgMaxIterations)
            throw MechanicsException(__PRETTY_FUNCTION__,"Maximum number of iterations exceeded.");

        if (structure->mRank == 0)
        {
            structure->GetLogger() << "\n*************************************\n";
            structure->GetLogger() << "       End Interface Problem             ";
            structure->GetLogger() << "\n*************************************\n\n";
        }



        MPI_Barrier(MPI_COMM_WORLD);
        lambda = x;

        VectorXd tmp = B * mSolver.solve(Btrans*x);
        MPI_Allreduce(MPI_IN_PLACE, tmp.data(), tmp.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        alphaGlobal = GtransGinv * G.transpose() * (displacementGap - tmp);

        MPI_Barrier(MPI_COMM_WORLD);

        std::vector<int> recvCount;
        std::vector<int> displs;
        MpiGatherRecvCountAndDispls(recvCount, displs, numRigidBodyModesLocal);
        VectorXd alphaLocal = alphaGlobal.segment(displs[structure->mRank], numRigidBodyModesLocal);

        StructureOutputBlockVector  delta_dof_dt0(structure->GetDofStatus());
        VectorXd delta_dof_active = mSolver.solve(residual_mod.Export() - Btrans * lambda);


        if (structure->mNumRigidBodyModes > 0)
            delta_dof_active -= rigidBodyModes * alphaLocal;

    int offset = 0;
    for (const auto& dof : activeDofSet)
    {
        int numActiveDofs = structure->GetNumActiveDofs(dof);
        delta_dof_dt0.J[dof] = delta_dof_active.segment(offset, structure->GetNumActiveDofs(dof));
        offset += numActiveDofs;
    }
        return delta_dof_dt0;


    }

    //! @brief perform the time integration
    //! @param rTimeDelta ... length of the simulation
    NuTo::eError Solve(double rTimeDelta) override
    {


        try
        {
            StructureFETI* structure = static_cast<StructureFETI*>(mStructure);
            structure->NodeBuildGlobalDofs(__PRETTY_FUNCTION__);

            structure->AssembleConnectivityMatrix();
            const SparseMatrix& B       = structure->GetConnectivityMatrix();
            const SparseMatrix& Btrans  = B.transpose();
            VectorXd lambda             = VectorXd::Zero(B.rows());

            double curTime  = mTime;
            double timeStep = mTimeStep;
            mStructure->SetPrevTime(curTime);
            mStructure->SetTime(curTime);

            mStructure->NodeBuildGlobalDofs(__PRETTY_FUNCTION__);

            CalculateStaticAndTimeDependentExternalLoad();

            const DofStatus& dofStatus = mStructure->GetDofStatus();


            if(mStepActiveDofs.empty())
            {
                mStepActiveDofs.push_back(mStructure->DofTypesGetActive());
            }
            else
            {
                for(unsigned int i=0; i<mStepActiveDofs.size();++i)
                {
                    if(mStepActiveDofs[i].empty())
                    {
                        throw MechanicsException(__PRETTY_FUNCTION__, "Calculation step " +std::to_string(i)+ " has no active DOFs.");
                    }
                }
            }


            /*---------------------------------*\
            |        Allocate Variables         |
            \*---------------------------------*/

            StructureOutputBlockMatrix  hessian0(dofStatus, true);

            StructureOutputBlockVector  delta_dof_dt0(dofStatus, true);
            StructureOutputBlockVector  trial_dof_dt0(dofStatus, true); // e.g. disp

            StructureOutputBlockVector  dof_dt0(dofStatus, true); // e.g. disp

            StructureOutputBlockVector  lastConverged_dof_dt0(dofStatus, true); // e.g. disp
            StructureOutputBlockVector  lastConverged_dof_dt1(dofStatus, true); // e.g. velocity
            StructureOutputBlockVector  lastConverged_dof_dt2(dofStatus, true); // e.g. accelerations

            StructureOutputBlockVector  extForce(dofStatus, true);
            StructureOutputBlockVector  intForce(dofStatus, true);
            StructureOutputBlockVector  residual(dofStatus, true);


            StructureOutputBlockVector  prevIntForce(dofStatus, true);
            StructureOutputBlockVector  prevExtForce(dofStatus, true);
            StructureOutputBlockVector  prevResidual(dofStatus, true);


            // for constraints
            // ---------------

            BlockFullVector<double> residual_mod(dofStatus);
            const auto& cmat = mStructure->GetConstraintMatrix();

            /*---------------------------------*\
            |    Declare and fill Output Maps   |
            \*---------------------------------*/

            // Declare output maps

            std::map<NuTo::eStructureOutput, NuTo::StructureOutputBase*> evalInternalGradient;
            std::map<NuTo::eStructureOutput, NuTo::StructureOutputBase*> evalInternalGradientAndHessian0;
            std::map<NuTo::eStructureOutput, NuTo::StructureOutputBase*> evalHessian0;

            evalInternalGradient                [eStructureOutput::INTERNAL_GRADIENT]   = &intForce;

            evalInternalGradientAndHessian0     [eStructureOutput::INTERNAL_GRADIENT]   = &intForce;
            evalInternalGradientAndHessian0     [eStructureOutput::HESSIAN0]            = &hessian0;

            evalHessian0                        [eStructureOutput::HESSIAN0]            = &hessian0;

            /*---------------------------------*\
            |    Declare and fill Input map     |
            \*---------------------------------*/

            ConstitutiveInputMap inputMap;
            inputMap[Constitutive::eInput::CALCULATE_STATIC_DATA] = std::make_unique<ConstitutiveCalculateStaticData>(
                    eCalculateStaticData::EULER_BACKWARD);

            ExtractDofValues(lastConverged_dof_dt0, lastConverged_dof_dt1, lastConverged_dof_dt2);

            UpdateAndGetAndMergeConstraintRHS(curTime, lastConverged_dof_dt0);

            PostProcess(residual);

            while (curTime < rTimeDelta)
            {

                mStructure->DofTypeActivateAll();

                // calculate Delta_BRhs and Delta_ExtForce
                auto bRHS           = UpdateAndGetAndMergeConstraintRHS(mTime, lastConverged_dof_dt0);
                auto prevExtForce   = CalculateCurrentExternalLoad(mTime);

                curTime += timeStep;
                SetTimeAndTimeStep(curTime, timeStep, rTimeDelta);     //check whether harmonic excitation, check whether curTime is too close to the time data
                mStructure->SetTime(curTime);

                auto deltaBRHS = UpdateAndGetConstraintRHS(curTime) - bRHS;
                auto extForce = CalculateCurrentExternalLoad(curTime);

                for (const auto& activeDofSet : mStepActiveDofs)
                {
                    mStructure->DofTypeSetIsActive(activeDofSet);

                    // ******************************************************
                    mStructure->Evaluate(inputMap, evalInternalGradientAndHessian0);
                    // ******************************************************

                    residual =  prevExtForce - extForce;
                    delta_dof_dt0.J.SetZero();
                    delta_dof_dt0.K = deltaBRHS;
                    residual += hessian0 * delta_dof_dt0;

                    residual.ApplyCMatrix(residual_mod, mStructure->GetConstraintMatrix());
                    hessian0.ApplyCMatrix(mStructure->GetConstraintMatrix());

//                    residual_mod = BlockFullVector<double>(residual_mod.Export() + Btrans * lambda, dofStatus);
                    mSolver.compute(hessian0.JJ.ExportToEigenSparseMatrix());

                    // -1 * residual_mod because it actually represents the rhs of the linear system...
                    delta_dof_dt0 = ConjugateProjectedGradientMethod(-1*residual_mod, activeDofSet, lambda);



                    delta_dof_dt0.K = deltaBRHS - mStructure->GetConstraintMatrix()*delta_dof_dt0.J;
                    dof_dt0 = lastConverged_dof_dt0 + delta_dof_dt0;


                    mStructure->NodeMergeDofValues(dof_dt0);



                    // ******************************************************
                    mStructure->Evaluate(inputMap, evalInternalGradient);
                    // ******************************************************


                    residual = intForce - extForce;
                    residual.ApplyCMatrix(residual_mod, cmat);


                    // WARNING! This function modifies residual_mod. Need to find a better way to
                    // calculate the norm and calculate the "true" residual_mod. The problem is that
                    // the internal forces at the interfaces are always nonzero but balance each other out.
                    BlockScalar normResidual = CalculateNormResidual(residual_mod, activeDofSet);

//                    residual_mod = BlockFullVector<double>(residual_mod.Export() + Btrans * lambda, dofStatus);
//                    BlockScalar normResidual = residual_mod.CalculateInfNorm();
//                    for (const auto& dof : activeDofSet)
//                        MPI_Allreduce(MPI_IN_PLACE,  &normResidual[dof], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);


                    MPI_Barrier(MPI_COMM_WORLD);

                    int iteration = 0;
                    while( normResidual > mToleranceResidual and iteration < mMaxNumIterations)
                    {


                            std::cout << "rank " << structure->mRank << " entered while loop" << std::endl;
                            MPI_Barrier(MPI_COMM_WORLD);

                        // ******************************************************
                        mStructure->Evaluate(inputMap, evalHessian0);
                        // ******************************************************

                        // ******************************************************
//                        Eigen::SparseMatrix<double> stiffnessMatrixSparse = ConvertToEigenSparseMatrix(hessian0.JJ(Node::eDof::DISPLACEMENTS,Node::eDof::DISPLACEMENTS));
                        mSolver.compute(hessian0.JJ.ExportToEigenSparseMatrix());

                        // -1 * residual_mod because it actually represents the rhs of the linear system...
                        delta_dof_dt0 = ConjugateProjectedGradientMethod(-1*residual_mod, activeDofSet,lambda);
                        delta_dof_dt0.K = cmat*delta_dof_dt0.J*(-1.);
                        // ******************************************************


                        //perform a line search
                        double alpha = 1;
                        BlockScalar trialNormResidual(mStructure->GetDofStatus());
                        do
                        {
                            //calculate line search trial state
                            trial_dof_dt0 = dof_dt0 + delta_dof_dt0 * alpha;

                            mStructure->NodeMergeDofValues(trial_dof_dt0);

                            // ******************************************************
                            mStructure->Evaluate(inputMap, evalInternalGradient);
                            // ******************************************************

                            residual = intForce - prevExtForce;
                            residual.ApplyCMatrix(residual_mod, cmat);


                            // WARNING! This function modifies residual_mod. Need to find a better way to
                            // calculate the norm and calculate the "true" residual_mod. The problem is that
                            // the internal forces at the interfaces are always nonzero but balance each other out.
                            trialNormResidual = CalculateNormResidual(residual_mod, activeDofSet);

//                            residual_mod = BlockFullVector<double>(residual_mod.Export() + Btrans * lambda, dofStatus);
//                            trialNormResidual = residual_mod.CalculateInfNorm();
//                            for (const auto& dof : activeDofSet)
//                                MPI_Allreduce(MPI_IN_PLACE,  &normResidual[dof], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

                            if (structure->mRank == 0)
                                mStructure->GetLogger() << "[Linesearch a=" << std::to_string(alpha).substr(0, 6) << "] Trial residual: " << trialNormResidual <<  "\n";

                            alpha*=0.5;

                        }
                        while(mPerformLineSearch && alpha > mMinLineSearchStep && trialNormResidual > (1. - alpha) * normResidual);

                        if (alpha > mMinLineSearchStep || !mPerformLineSearch)
                        {
                            //improvement is achieved, go to next Newton step
                            dof_dt0 = trial_dof_dt0;

                            normResidual = trialNormResidual;

                            PrintInfoIteration(normResidual,iteration);
                            iteration++;
                        }
                        else
                        {
                            //and leave
                            iteration = mMaxNumIterations;
                        }

                    } // end of while(normResidual<mToleranceForce && iteration<mMaxNumIterations)



                    MPI_Barrier(MPI_COMM_WORLD);

                    if (normResidual < mToleranceResidual)
                    {
                        //converged solution
                        mStructure->ElementTotalUpdateStaticData();

                        //store converged step
                        lastConverged_dof_dt0 = dof_dt0;

                        prevResidual = residual;

                        mStructure->NodeMergeDofValues(dof_dt0);

                        //update structure time
                        mStructure->SetPrevTime(curTime);

                        mTime+=timeStep;

                        if (structure->mRank == 0)
                        {
                            mStructure->GetLogger() << "Convergence after " << iteration << " iterations at time " << mTime << " (timestep " << timeStep << ").\n";
                            mStructure->GetLogger() << "Residual: \t" << normResidual << "\n\n";
                        }

                        //eventually increase next time step
                        if (mAutomaticTimeStepping && iteration<0.25*mMaxNumIterations)
                        {
                            timeStep*=1.5;
                            if (timeStep>mMaxTimeStep)
                                timeStep = mMaxTimeStep;
                        }

                        //perform Postprocessing
                        PostProcess(prevResidual);

                        if (mCallback && mCallback->Exit(*mStructure))
                            return NuTo::eError::SUCCESSFUL;

                    }
                    else
                    {
                        if (structure->mRank == 0)
                            mStructure->GetLogger() << "No convergence with timestep " << timeStep << "\n\n";

                        //no convergence
                        if (mAutomaticTimeStepping)
                        {
                            //no convergence, reduce the time step and start from scratch
                            curTime -= timeStep;
                            timeStep *= 0.5;
                            if (timeStep < mMinTimeStep) {
                                mStructure->GetLogger() << "The minimal time step achieved, the actual time step is " << timeStep << "\n";
                                throw MechanicsException(__PRETTY_FUNCTION__, "No convergence, the current time step is too short.");
                            }
                        }
                        else
                        {
                            throw MechanicsException(__PRETTY_FUNCTION__, "No convergence with the current maximum number of iterations, either use automatic time stepping, reduce the time step or the minimal line search cut back factor.");
                        }

                    }

                } // end for mStepActiveDofs


            } // end while



        }
        catch (MechanicsException& e)
        {
            e.AddMessage(__PRETTY_FUNCTION__, " ERROR performing FETI.");
            throw e;
        }


        return NuTo::eError::SUCCESSFUL;

    }

private:
    Eigen::SparseQR<Eigen::SparseMatrix<double>,Eigen::COLAMDOrdering<int>> mSolver;
    const double    mCpgTolerance     = 1.0e-6;
    const int       mCpgMaxIterations = 50;
};
}// namespace NuTo
