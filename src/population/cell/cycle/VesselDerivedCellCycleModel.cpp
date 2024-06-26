/*

Copyright (c) 2005-2013, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#include "VesselDerivedCellCycleModel.hpp"
#include "Exception.hpp"

#include "MacrophageMutationState.hpp"
#include "VesselCellMutationState.hpp"
#include "StalkCellMutationState.hpp"
#include "TipCellMutationState.hpp"

VesselDerivedCellCycleModel::VesselDerivedCellCycleModel()
    :AbstractCellCycleModel()
{
}

AbstractCellCycleModel* VesselDerivedCellCycleModel::CreateCellCycleModel()
{
    // Create a new cell-cycle model
    VesselDerivedCellCycleModel* p_model = new VesselDerivedCellCycleModel();

    /*
     * Set each member variable of the new cell-cycle model that inherits
     * its value from the parent.
     *
     * Note 1: some of the new cell-cycle model's member variables (namely
     * mBirthTime, mCurrentCellCyclePhase, mReadyToDivide) will already have been
     * correctly initialized in its constructor.
     *
     * Note 2: one or more of the new cell-cycle model's member variables
     * may be set/overwritten as soon as InitialiseDaughterCell() is called on
     * the new cell-cycle model.
     *
     * Note 3: the member variable mDimension remains unset, since this cell-cycle
     * model does not need to know the spatial dimension, so if we were to call
     * SetDimension() on the new cell-cycle model an exception would be triggered;
     * hence we do not set this member variable.
     */
    p_model->SetBirthTime(mBirthTime);
    p_model->SetMinimumGapDuration(mMinimumGapDuration);
    p_model->SetStemCellG1Duration(mStemCellG1Duration);
    p_model->SetTransitCellG1Duration(mTransitCellG1Duration);
    p_model->SetSDuration(mSDuration);
    p_model->SetG2Duration(mG2Duration);
    p_model->SetMDuration(mMDuration);

    return p_model;
}

void VesselDerivedCellCycleModel::InitialiseDaughterCell()
{
	double mProbOfMacrophageFromVessel = 0.000; // Prob of a vessel cell giving rise to a macrophage in 1 hour
 	double mProbOfTipFromVessel = 0.0005; // Prob of a vessel cell giving rise to a tip cell in 1 hour

	if ( mpCell->GetMutationState()->IsType<VesselCellMutationState>() )
	{
		double probability_of_tip = mProbOfTipFromVessel/(mProbOfTipFromVessel + mProbOfMacrophageFromVessel);

		double p_test = RandomNumberGenerator::Instance()->ranf();

		if (p_test < probability_of_tip)
		{
			mpCell->SetMutationState(CellPropertyRegistry::Instance()->Get<TipCellMutationState>());
			mpCell->GetCellData()->SetItem("StalkId",mpCell->GetCellId());
		}
		else
		{
			mpCell->SetMutationState(CellPropertyRegistry::Instance()->Get<MacrophageMutationState>());
		}
	}
	else if ( mpCell->GetMutationState()->IsType<StalkCellMutationState>() )
	{
		mpCell->SetMutationState(CellPropertyRegistry::Instance()->Get<TipCellMutationState>());
	}
	else
	{
		assert(mpCell->GetMutationState()->IsType<TipCellMutationState>() ||
			   mpCell->GetMutationState()->IsType<MacrophageMutationState>() );
		NEVER_REACHED;
	}
}


bool VesselDerivedCellCycleModel::ReadyToDivide()
{
	double mProbOfMacrophageFromVessel = 0.0000; // Prob of a vessel cell giving rise to a macrophage in 1 hour
 	double mProbOfTipFromVessel = 0.0005; // Prob of a vessel cell giving rise to a tip cell in 1 hour

 	double mProbOfTipMoving = 0.05; // Prob of a Tip cell moving in 1 hour


    assert(mpCell != NULL);

    if (!mReadyToDivide)
    {
		// UpdateCellCyclePhase();

		if ( mpCell->GetMutationState()->IsType<VesselCellMutationState>() )
		{
			double probability_of_division = (mProbOfMacrophageFromVessel + mProbOfTipFromVessel) * SimulationTime::Instance()->GetTimeStep();

			double p_test = RandomNumberGenerator::Instance()->ranf();

			if (p_test < probability_of_division)
			{
				mReadyToDivide = true;
			}
		}
		else if ( mpCell->GetMutationState()->IsType<TipCellMutationState>() )
		{
			double probability_of_division = mProbOfTipMoving * SimulationTime::Instance()->GetTimeStep();

			double p_test = RandomNumberGenerator::Instance()->ranf();

			if (p_test < probability_of_division)
			{
				mReadyToDivide = true;
				mpCell->SetMutationState(CellPropertyRegistry::Instance()->Get<StalkCellMutationState>());
			}
		}
		else
		{
			assert(mpCell->GetMutationState()->IsType<StalkCellMutationState>() ||
					mpCell->GetMutationState()->IsType<MacrophageMutationState>() );
		}



//        if ( (mCurrentCellCyclePhase != G_ZERO_PHASE) &&
//             (GetAge() >= GetMDuration() + GetG1Duration() + GetSDuration() + GetG2Duration()) )
//        {
//            mReadyToDivide = true;
//        }
    }
    return mReadyToDivide;
}

void VesselDerivedCellCycleModel::OutputCellCycleModelParameters(out_stream& rParamsFile)
{
    // No new parameters to output, so just call method on direct parent class
    AbstractCellCycleModel::OutputCellCycleModelParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(VesselDerivedCellCycleModel)
