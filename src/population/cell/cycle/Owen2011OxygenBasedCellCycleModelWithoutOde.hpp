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

#ifndef OWEN2011OXYGENBASEDCELLCYCLEMODELWITHOUTODE_HPP_
#define OWEN2011OXYGENBASEDCELLCYCLEMODELWITHOUTODE_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractSimpleCellCycleModel.hpp"
#include "AbstractCellMutationState.hpp"
#include "CancerCellMutationState.hpp"
#include "QuiescentCancerCellMutationState.hpp"

class Owen2011OxygenBasedCellCycleModelWithoutOde : public AbstractSimpleCellCycleModel
{
private:

	 /**
	   * Constants for the Owen et al. (2011) model
	   */

	 /** Tmin represents the minimum period of the cell cycle (in mins) */
	    double mTmin;

	 /** C0 represents the oxygen concentration at which the speed is half maximal (in mmHg) */
	    double mC0;

	 /** phi (phase representative) */
	     double mPhi;

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the cell-cycle model, never used directly - boost uses this.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractSimpleCellCycleModel>(*this);
        archive & mCurrentQuiescentDuration;
        archive & mCurrentQuiescenceOnsetTime;
        archive & mEnterQuiescenceOxygenConcentration;
        archive & mLeaveQuiescenceOxygenConcentration;
        archive & mCriticalQuiescentDuration;
        archive & mPhi;
        archive &mDivideTime;
     }

protected:

    /**
     * How long the current period of quiescence has lasted.
     * Has units of mins
     */
    double mCurrentQuiescentDuration;

    /**
     * The time when the current period of quiescence began.
     */
    double mCurrentQuiescenceOnsetTime;

    /**
     * Oxygen concentration below which cells enter quiescence.
     * A prolonged period of quiescence causes the cell to become apoptotic.
     */
    double mEnterQuiescenceOxygenConcentration;

    /**
     * Oxygen concentration above which cells leave their state of being quiescent
     */
    double mLeaveQuiescenceOxygenConcentration;

    /**
     * Critical quiescent duration.
     * Has units of mins.
     */
    double mCriticalQuiescentDuration;

    /**
     * The time at which division occurs.
     */
    double mDivideTime;

public:

    /**
      * Constructor.
      */
    Owen2011OxygenBasedCellCycleModelWithoutOde();

    /**
      * Update cell-cycle phase.
      */
    void UpdateCellCyclePhase();

      /**
     * Reset cell-cycle model by calling AbstractSimpleCellCycleModel::ResetForDivision()
     * and setting initial conditions for phi (phase representative).
     */
    void ResetForDivision();

    /**
     * Overridden builder method to create new copies of
     * this cell-cycle model.
     */
    AbstractCellCycleModel* CreateCellCycleModel();

    /**
     * Update the duration for which the cell has been quiescent.
     */
    void UpdateQuiescentDuration();

    /**
     * Check if the oxygen concentration of the cell is below the EnterQuiescenceOxygenConcentration.
     * If it is true the label cells and return true else return false without labelling.
     */
    void CheckAndLabelCell();

    /**
     * @return mCurrentQuiescentDuration
     */
    double GetCurrentQuiescentDuration();

    /**
     * @return mCurrentQuiescenceOnsetTime
     */
    double GetCurrentQuiescenceOnsetTime();

    /**
     * @return mEnterQuiescenceOxygenConcentration
     */
    double GetEnterQuiescenceOxygenConcentration();

    /**
     * Set method for mEnterQuiescenceOxygenConcentration.
     *
     * @param enterQuiescenceOxygenConcentration the new value of mEnterQuiescenceOxygenConcentration
     */
    void SetEnterQuiescenceOxygenConcentration(double enterQuiescenceOxygenConcentration);

    /**
     * @return mLeaveQuiescenceOxygenConcentration
     */
    double GetLeaveQuiescenceOxygenConcentration();

    /**
     * Set method for mLeaveQuiescenceOxygenConcentration.
     *
     * @param leaveQuiescenceOxygenConcentration the new value of mLeaveQuiescenceOxygenConcentration
     */
    void SetLeaveQuiescenceOxygenConcentration(double leaveQuiescenceOxygenConcentration);

    /**
     * @return mCriticalQuiescentDuration
     */
    double GetCriticalQuiescentDuration();

    /**
     * Set method for mCriticalQuiescentDuration.
     *
     * @param criticalQuiescentDuration the new value of mCriticalQuiescentDuration
     */
    void SetCriticalQuiescentDuration(double criticalQuiescentDuration);

    /**
     * Set method for mCurrentQuiescenceOnsetTime.
     *
     * @param currentQuiescenceOnsetTime the new value of mCurrentQuiescenceOnsetTime
     */
    void SetCurrentQuiescenceOnsetTime(double currentQuiescenceOnsetTime);

    /**
      * Set the duration of the cell's G1 phase.
      */
    void SetG1Duration();


    /**
     * Get the duration of the cell's S phase.
     */
    double GetSDuration();

    /**
     * Get the duration of the cell's G2 phase.
     */
    double GetG2Duration();

    /**
     * Get the duration of the cell's M phase.
     */
    double GetMDuration();

    /**
      * Get the value of phi.
      */
    double GetPhi();

    /**
     * Update value of phi
     */
    double UpdatePhi(double current_time, double oxygen, double dt);

    /**
     * Update oxygen concentration
     */
    double GetOxygenConcentration();

    /**
      * Set the value of phi.
      */
    void SetPhi(double phi);

    /**
     * Set the value of mDivideTime.
     */
    void SetDivideTime(double dividetime);

    /**
     * Get the value of mDivideTime.
     */
    double GetDivideTime();

    /**
     * Outputs cell cycle model parameters to file.
     *
public:
     * @param rParamsFile the file stream to which the parameters are output
     */
    virtual void OutputCellCycleModelParameters(out_stream& rParamsFile);
};

// Declare identifier for the serializer
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(Owen2011OxygenBasedCellCycleModelWithoutOde)

#endif /*OWEN2011OXYGENBASEDCELLCYCLEMODELWITHOUTODE_HPP_*/
