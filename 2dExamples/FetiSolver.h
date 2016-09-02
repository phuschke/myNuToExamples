#include "nuto/mechanics/timeIntegration/TimeIntegrationBase.h"
#include "nuto/base/ErrorEnum.h"
#include "nuto/mechanics/structures/StructureBase.h"
#include "nuto/mechanics/MechanicsException.h"
#include "nuto/base/Timer.h"

#include "nuto/mechanics/structures/StructureBase.h"
#include "nuto/mechanics/structures/unstructured/Structure.h"
#include "nuto/mechanics/structures/StructureBaseEnum.h"
#include "nuto/mechanics/structures/StructureOutputBlockMatrix.h"
#include "nuto/mechanics/structures/StructureOutputDummy.h"
#include "nuto/math/SparseDirectSolverMUMPS.h"

#include "nuto/base/CallbackInterface.h"
#include "nuto/math/SparseMatrixCSRGeneral.h"

#include "nuto/mechanics/constitutive/ConstitutiveEnum.h"
#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveCalculateStaticData.h"
#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveIOMap.h"
#include "nuto/mechanics/constitutive/inputoutput/ConstitutiveTimeStep.h"
#include "nuto/mechanics/constitutive/staticData/ConstitutiveStaticDataGradientDamage.h"


namespace NuTo
{
class FetiSolver : public TimeIntegrationBase
{
public:
    FetiSolver(StructureBase* rStructure) : TimeIntegrationBase (rStructure)
    {

    }

    //! @brief perform the time integration
    //! @param rTimeDelta ... length of the simulation
    NuTo::eError Solve(double rTimeDelta) override
    {
        if (mAutomaticTimeStepping)
            throw MechanicsException(__PRETTY_FUNCTION__, "Automatic time stepping not supported.");

        try
        {
            NuTo::Timer timerDebug("Init", true, mStructure->GetLogger());

            mStructure->NodeBuildGlobalDofs(__PRETTY_FUNCTION__);

            CalculateStaticAndTimeDependentExternalLoad();

            const DofStatus& dofStatus = mStructure->GetDofStatus();

            if (mStepActiveDofs.empty())
                throw MechanicsException(__PRETTY_FUNCTION__, "Define a set of active dofs for each calculation step. ");

            // deactivate all dof types
            mStructure->DofTypeDeactivateAll();



            /*---------------------------------*\
            |        Allocate Variables         |
            \*---------------------------------*/

            StructureOutputBlockMatrix  hessian0(dofStatus, true);

            StructureOutputBlockVector  delta_dof_dt0(dofStatus, true);

            StructureOutputBlockVector  dof_dt0(dofStatus, true); // e.g. disp
            StructureOutputBlockVector  lastConverged_dof_dt0(dofStatus, true); // e.g. disp

            StructureOutputBlockVector  extForce(dofStatus, true);
            StructureOutputBlockVector  intForce(dofStatus, true);
            StructureOutputBlockVector  residual(dofStatus, true);


            // for constraints
            // ---------------

            BlockFullVector<double> residual_mod(dofStatus);

            /*---------------------------------*\
            |    Declare and fill Output Maps   |
            \*---------------------------------*/

            // Declare output maps
            std::map<NuTo::eStructureOutput, NuTo::StructureOutputBase*> evalUpdateStaticData;
            StructureOutputDummy dummy;
            evalUpdateStaticData                [eStructureOutput::UPDATE_STATIC_DATA] = &dummy;

            std::map<NuTo::eStructureOutput, NuTo::StructureOutputBase*> evalInternalGradient;
            evalInternalGradient                [eStructureOutput::INTERNAL_GRADIENT] = &intForce;

            std::map<NuTo::eStructureOutput, NuTo::StructureOutputBase*> evalInternalGradientAndHessian0;
            evalInternalGradientAndHessian0     [eStructureOutput::INTERNAL_GRADIENT] = &intForce;
            evalInternalGradientAndHessian0     [eStructureOutput::HESSIAN0] = &hessian0;

            /*---------------------------------*\
            |    Declare and fill Input map     |
            \*---------------------------------*/



            ConstitutiveInputMap input;
            input[Constitutive::eInput::TIME_STEP] = std::make_unique<ConstitutiveTimeStep>(2);
            auto& timeStep = *static_cast<ConstitutiveTimeStep*>(input[Constitutive::eInput::TIME_STEP].get());
            input[Constitutive::eInput::CALCULATE_STATIC_DATA] = std::make_unique<ConstitutiveCalculateStaticData>(eCalculateStaticData::USE_PREVIOUS ,1);
            auto& calculateStaticData = *static_cast<ConstitutiveCalculateStaticData*>(input[Constitutive::eInput::CALCULATE_STATIC_DATA].get());

            timeStep.SetCurrentTimeStep(mTimeStep);


            PostProcess(residual);

            int numAcceptedIterations = 0;
            int numRejectedIterations = 0;

            while (mTime < rTimeDelta)
            {
                timerDebug.Reset("\033[1;31m Iteration " + std::to_string(numAcceptedIterations) + " at current time: " + std::to_string(mTime) + "\033[0m");

                mStructure->DofTypeActivateAll();


                // calculate Delta_BRhs and Delta_ExtForce
                auto bRHS = UpdateAndGetConstraintRHS(mTime);
                auto prevExtForce = CalculateCurrentExternalLoad(mTime);


                mTime += mTimeStep;


                auto deltaBRHS = UpdateAndGetConstraintRHS(mTime) - bRHS;
                auto extForce = CalculateCurrentExternalLoad(mTime);

                std::cout << "TimeStep: " << mTimeStep << std::endl;



                for (const auto& activeDofSet : mStepActiveDofs)
                {
                    mStructure->DofTypeSetIsActive(activeDofSet);



                    calculateStaticData.SetCalculateStaticData(eCalculateStaticData::USE_PREVIOUS);
                    calculateStaticData.SetIndexOfPreviousStaticData(0);
                    mStructure->Evaluate(input, evalInternalGradientAndHessian0);

                    delta_dof_dt0.J.SetZero();
                    delta_dof_dt0.K = deltaBRHS;

                    residual = (intForce - extForce) - prevExtForce - extForce;
                    residual += hessian0 * delta_dof_dt0;
                    residual.ApplyCMatrix(mStructure->GetConstraintMatrix());

                    //                    mStructure->GetLogger() << "Initial trial residual:               " << residual.J.CalculateInfNorm() << "\n";

                    hessian0.ApplyCMatrix(mStructure->GetConstraintMatrix());


                    delta_dof_dt0.J =  mStructure->SolveBlockSystem(hessian0.JJ, residual.J);
                    delta_dof_dt0.K = deltaBRHS - mStructure->GetConstraintMatrix()*delta_dof_dt0.J;

                    dof_dt0 += delta_dof_dt0;
                    mStructure->NodeMergeDofValues(dof_dt0);

                } // end for mStepActiveDofs


                if (mCallback && mCallback->Exit(*mStructure))
                    return eError::SUCCESSFUL;

                bool acceptSolution = true;


                if (acceptSolution)
                {
                    // save the new implicit history variables
                    calculateStaticData.SetCalculateStaticData(eCalculateStaticData::EULER_BACKWARD);
                    calculateStaticData.SetIndexOfPreviousStaticData(1);
                    mStructure->Evaluate(input, evalUpdateStaticData);

                    if (mTime + mTimeStep > rTimeDelta)
                        mTimeStep = rTimeDelta - mTime;

                    //                std::cout << "Bfore SaveStaticData" << std::endl;
                    //                mCallback->Exit(*mStructure);
                    mStructure->ElementTotalSaveStaticData();   // shift static data by one to the past
                    timeStep.SetCurrentTimeStep(mTimeStep);     // shift time steps by one to the past

                    //                std::cout << "after SaveStaticData" << std::endl;
                    //                mCallback->Exit(*mStructure);

                    mStructure->DofTypeActivateAll();
                    lastConverged_dof_dt0 = dof_dt0;



                    calculateStaticData.SetCalculateStaticData(eCalculateStaticData::USE_PREVIOUS);
                    mStructure->Evaluate(input, evalInternalGradient);
                    residual = intForce - extForce;
                    PostProcess(residual);

                    numAcceptedIterations++;
                }
                else
                {
                    // do not save static data
                    // do not shift the time but set the new time step
                    timeStep[0] = mTimeStep;


                    // restore the old solution
                    mStructure->DofTypeActivateAll();
                    dof_dt0 = lastConverged_dof_dt0;
                    mStructure->NodeMergeDofValues(dof_dt0);
                    mTime -= mTimeStep;

                    numRejectedIterations++;
                }

            } // end while


            std::cout << "["<<__FUNCTION__<<"] Number of accepted iterations: " << numAcceptedIterations << std::endl;
            std::cout << "["<<__FUNCTION__<<"] Number of rejected iterations: " << numRejectedIterations << std::endl;



        }
        catch (MechanicsException& e)
        {
            e.AddMessage(__PRETTY_FUNCTION__, " ERROR performing FETI.");
            throw e;
        }


        return NuTo::eError::SUCCESSFUL;

    }

private:
};
}// namespace NuTo
