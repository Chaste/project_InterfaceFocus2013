/*

Copyright (c) 2005-2012, University of Oxford.
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

#include "Owen2011OxygenBasedCellCycleModelWithoutOde.hpp"
#include "ApoptoticCellProperty.hpp"
#include "CellPropertyRegistry.hpp"
#include "CellLabel.hpp"
#include "Exception.hpp"
#include "StemCellProliferativeType.hpp"
#include "TransitCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"

#include "MacrophageMutationState.hpp"


Owen2011OxygenBasedCellCycleModelWithoutOde::Owen2011OxygenBasedCellCycleModelWithoutOde()
    :   mPhi(0.0),
        mCurrentQuiescentDuration(0.0),
        mCurrentQuiescenceOnsetTime(0.0),
        mEnterQuiescenceOxygenConcentration(8.9),
        mLeaveQuiescenceOxygenConcentration(9.8),
        mCriticalQuiescentDuration(4000)
{
	mCurrentCellCyclePhase = G_ONE_PHASE;
}

void Owen2011OxygenBasedCellCycleModelWithoutOde::ResetForDivision()
{
	AbstractCellCycleModel::ResetForDivision();
	mPhi=0;
	mBirthTime=mDivideTime;
	SetG1Duration();
}

double Owen2011OxygenBasedCellCycleModelWithoutOde::GetCurrentQuiescentDuration()
{
    return mCurrentQuiescentDuration;
}

double Owen2011OxygenBasedCellCycleModelWithoutOde::GetCurrentQuiescenceOnsetTime()
{
    return mCurrentQuiescenceOnsetTime;
}

void Owen2011OxygenBasedCellCycleModelWithoutOde::UpdateCellCyclePhase()
{
    // mG1Duration is set when value of phi equals 1

	if(mpCell->GetMutationState()->IsType<CancerCellMutationState>())
	{
		CheckAndLabelCell();
	}

	if(mpCell->GetMutationState()->IsType<QuiescentCancerCellMutationState>())
	{
	    UpdateQuiescentDuration();
	}

    if(mpCell->GetMutationState()->IsType<MacrophageMutationState>())
    {
        mCurrentCellCyclePhase = G_ZERO_PHASE;
    }

    double time_duration= SimulationTime::Instance()->GetTime() - mBirthTime;
	double oxygen = GetOxygenConcentration();
	double dt = SimulationTime::Instance()->GetTimeStep();
	mPhi= UpdatePhi(time_duration, oxygen, dt);

    // Update G1 duration when mPhi equals one
    if(mPhi>=1)
    {
      mG1Duration=time_duration;
      mDivideTime= SimulationTime::Instance()->GetTime();
    }
}

double Owen2011OxygenBasedCellCycleModelWithoutOde::GetOxygenConcentration()
{
	// Get cell's oxygen concentration
	double oxygen_concentration;

	oxygen_concentration = mpCell->GetCellData()->GetItem("oxygen");

	return (oxygen_concentration);
}

AbstractCellCycleModel* Owen2011OxygenBasedCellCycleModelWithoutOde::CreateCellCycleModel()
{
    // Create a new cell-cycle model
    Owen2011OxygenBasedCellCycleModelWithoutOde* p_model = new Owen2011OxygenBasedCellCycleModelWithoutOde();

    /*
     * Set each member variable of the new cell-cycle model that inherits
     * its value from the parent.
     *
     * Note 1: some of the new cell-cycle model's member variables (namely
     * mBirthTime, mCurrentCellCyclePhase, mReadyToDivide, mDt, mpOdeSolver)
     * will already have been correctly initialized in its constructor.
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
    p_model->SetEnterQuiescenceOxygenConcentration(mEnterQuiescenceOxygenConcentration);
    p_model->SetLeaveQuiescenceOxygenConcentration(mLeaveQuiescenceOxygenConcentration);
    p_model->SetCriticalQuiescentDuration(mCriticalQuiescentDuration);
    p_model->SetCurrentQuiescenceOnsetTime(mCurrentQuiescenceOnsetTime);

    return p_model;
}

void Owen2011OxygenBasedCellCycleModelWithoutOde::UpdateQuiescentDuration()
{
	assert(mpCell->GetMutationState()->IsType<QuiescentCancerCellMutationState>());
	assert(!(mpCell->HasCellProperty<ApoptoticCellProperty>()));
    assert(!mpCell->HasApoptosisBegun());

    // Get cell's oxygen concentration
    double oxygen_concentration=GetOxygenConcentration();

    if (oxygen_concentration < mLeaveQuiescenceOxygenConcentration)
    {
        // Update the duration of the current period of hypoxia
        mCurrentQuiescentDuration = (SimulationTime::Instance()->GetTime() - mCurrentQuiescenceOnsetTime);

        if (mCurrentQuiescentDuration > mCriticalQuiescentDuration)
        {
            mpCell->AddCellProperty(CellPropertyRegistry::Instance()->Get<CellLabel>());
        }
    }
    else
    {
        // Reset the cell's quiescent duration.
    	 mCurrentQuiescentDuration = 0.0;
    	 mCurrentQuiescenceOnsetTime = 0.0;
    	 mCurrentCellCyclePhase = G_ONE_PHASE;
    	 mpCell->SetMutationState(CellPropertyRegistry::Instance()->Get<CancerCellMutationState>());
    }
}

void Owen2011OxygenBasedCellCycleModelWithoutOde::CheckAndLabelCell()
{
	assert(mpCell->GetMutationState()->IsType<CancerCellMutationState>());
	// Get cell's oxygen concentration
	double oxygen_concentration=GetOxygenConcentration();

	if(oxygen_concentration<mEnterQuiescenceOxygenConcentration)
	{
		mpCell->SetMutationState(CellPropertyRegistry::Instance()->Get<QuiescentCancerCellMutationState>());
		assert(mpCell->GetMutationState()->IsType<QuiescentCancerCellMutationState>());
	    assert(!(mpCell->GetMutationState()->IsType<CancerCellMutationState>()));
		mCurrentQuiescenceOnsetTime=SimulationTime::Instance()->GetTime();
		mCurrentCellCyclePhase = G_ZERO_PHASE;
	}
}

double Owen2011OxygenBasedCellCycleModelWithoutOde::GetEnterQuiescenceOxygenConcentration()
{
    return mEnterQuiescenceOxygenConcentration;
}

void Owen2011OxygenBasedCellCycleModelWithoutOde::SetEnterQuiescenceOxygenConcentration(double enterQuiescenceOxygenConcentration)
{
    assert(enterQuiescenceOxygenConcentration>=0.0);
    mEnterQuiescenceOxygenConcentration = enterQuiescenceOxygenConcentration;
}

double Owen2011OxygenBasedCellCycleModelWithoutOde::GetLeaveQuiescenceOxygenConcentration()
{
    return mLeaveQuiescenceOxygenConcentration;
}

void Owen2011OxygenBasedCellCycleModelWithoutOde::SetLeaveQuiescenceOxygenConcentration(double leaveQuiescenceOxygenConcentration)
{
    assert(leaveQuiescenceOxygenConcentration >= 0.0);
    mLeaveQuiescenceOxygenConcentration = leaveQuiescenceOxygenConcentration;
}

double Owen2011OxygenBasedCellCycleModelWithoutOde::GetCriticalQuiescentDuration()
{
    return mCriticalQuiescentDuration;
}

void Owen2011OxygenBasedCellCycleModelWithoutOde::SetCriticalQuiescentDuration(double criticalQuiescentDuration)
{
    assert(criticalQuiescentDuration >= 0.0);
    mCriticalQuiescentDuration = criticalQuiescentDuration;
}

void Owen2011OxygenBasedCellCycleModelWithoutOde::SetCurrentQuiescenceOnsetTime(double currentQuiescenceOnsetTime)
{
    assert(currentQuiescenceOnsetTime >= 0.0);
    mCurrentQuiescenceOnsetTime = currentQuiescenceOnsetTime;
}

void Owen2011OxygenBasedCellCycleModelWithoutOde::SetPhi(double phi)
{
	mPhi= phi;
}

double Owen2011OxygenBasedCellCycleModelWithoutOde::GetPhi()
{
	return mPhi;
}

double Owen2011OxygenBasedCellCycleModelWithoutOde::UpdatePhi(double current_time, double oxygen, double dt)
{

	// Parameter values are taken from the Owen et al. (2011) paper. Oxygen concentration is measured in mmHg.
    //It is assumed that C=0.045 corresponds to typical dimensional perivascular oxygen tension of 20mmHg.

    	if (mpCell->GetMutationState()->IsType<CancerCellMutationState>())
        {
        	 mTmin = 1600;
        	 mC0 = (1.4*0.045)/20;
        }
        else
        {
        	mTmin=3000;
//        	mTmin=100;
        	mC0=(3*0.045)/20;
        }

    if(current_time!=0)
	{
	 	current_time -= dt;
    	if((current_time-0.0)<1e-6)
    	{
    		current_time=0.0;
    	}
  	 	mPhi= (oxygen/(mTmin * (mC0 + oxygen)) * dt * 60) + UpdatePhi(current_time, oxygen, dt);
	}
	else
	{
		return 0;
	}

	SetPhi(mPhi);
    return mPhi;
}

void Owen2011OxygenBasedCellCycleModelWithoutOde::SetG1Duration()
{
	assert(mpCell != NULL);
	mG1Duration = DBL_MAX;

	/*if (mpCell->GetCellProliferativeType()->IsType<StemCellProliferativeType>())
	{
	   mG1Duration =DBL_MAX;
	}
	else if (mpCell->GetCellProliferativeType()->IsType<TransitCellProliferativeType>())
	{
	   mG1Duration = DBL_MAX;
	}
	else if ( mpCell->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>())
	{
	   mG1Duration = DBL_MAX;
	}
	else
	{
		   NEVER_REACHED;
	}*/
 }

double Owen2011OxygenBasedCellCycleModelWithoutOde::GetSDuration()
{
    /**
     * Tyson & Novak pretends it is running ODEs in just G1,
     * but they really represent the whole cell cycle, so
     * we set the other phases to zero.
     */
    return 0.0;
}

double Owen2011OxygenBasedCellCycleModelWithoutOde::GetG2Duration()
{
    /**
     * Tyson & Novak pretends it is running ODEs in just G1,
     * but they really represent the whole cell cycle so
     * we set the other phases to zero.
     */
    return 0.0;
}

double Owen2011OxygenBasedCellCycleModelWithoutOde::GetMDuration()
{
    /**
     * Tyson & Novak pretends it is running ODEs in just G1,
     * but they really represent the whole cell cycle so
     * we set the other phases to zero.
     */
    return 0.0;
}


void Owen2011OxygenBasedCellCycleModelWithoutOde::SetDivideTime(double dividetime)
{
	mDivideTime = dividetime;
}

double Owen2011OxygenBasedCellCycleModelWithoutOde::GetDivideTime()
{
	return mDivideTime;
}

void Owen2011OxygenBasedCellCycleModelWithoutOde::OutputCellCycleModelParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<EnterQuiescenceOxygenConcentration>" << mEnterQuiescenceOxygenConcentration << "</EnterQuiescenceOxygenConcentration>\n";
    *rParamsFile << "\t\t\t<LeaveQuiescenceOxygenConcentration>" << mLeaveQuiescenceOxygenConcentration << "</LeaveQuiescenceOxygenConcentration>\n";
    *rParamsFile << "\t\t\t<CriticalQuiescentDuration>" << mCriticalQuiescentDuration << "</CriticalQuiescentDuration>\n";

    // Call method on direct parent class
    AbstractSimpleCellCycleModel::OutputCellCycleModelParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(Owen2011OxygenBasedCellCycleModelWithoutOde)
