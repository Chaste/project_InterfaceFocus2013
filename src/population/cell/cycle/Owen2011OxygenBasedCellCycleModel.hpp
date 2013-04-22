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

#ifndef OWEN2011OXYGENBASEDCELLCYCLEMODEL_HPP_
#define OWEN2011OXYGENBASEDCELLCYCLEMODEL_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include <vector>

#include "AbstractOdeBasedCellCycleModel.hpp"
#include "AbstractCellMutationState.hpp"
#include "CancerCellMutationState.hpp"
#include "Owen2011OxygenBasedCellCycleOdeSystem.hpp"

/**
 * Oxygen-dependent ODE-based cell-cycle model. Published by Owen et al. 2011
 */
class Owen2011OxygenBasedCellCycleModel : public AbstractOdeBasedCellCycleModel
{
private:

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the cell-cycle model and ODE system.
     *
     * @param archive the archive
     * @param version the archive version
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractOdeBasedCellCycleModel>(*this);
        archive & mCurrentQuiescentDuration;
        archive & mCurrentQuiescenceOnsetTime;
        archive & mEnterQuiescenceOxygenConcentration;
        archive & mLeaveQuiescenceOxygenConcentration;
        archive & mCriticalQuiescentDuration;
    }

    /**
     * Adjust any ODE parameters needed before solving until currentTime.
     *
     * @param currentTime  the time up to which the system will be solved.
     */
    void AdjustOdeParameters(double currentTime);

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

public:

    /**
     * Default constructor.
     *
     * @param pOdeSolver An optional pointer to a cell-cycle model ODE solver object (allows the use of different ODE solvers)
     */
    Owen2011OxygenBasedCellCycleModel(boost::shared_ptr<AbstractCellCycleModelOdeSolver> pOdeSolver = boost::shared_ptr<AbstractCellCycleModelOdeSolver>());

    /**
      * Update cell-cycle phase.
      */
    void UpdateCellCyclePhase();

   /**
     * Resets the oxygen-based model to the start of the cell cycle
     * (this model does not cycle naturally). Cells are given a new
     * birth time and cell cycle proteins are reset. Note that the
     * oxygen concentration maintains its current value.
     *
     * Should only be called by the Cell Divide() method.
     */
    virtual void ResetForDivision();

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
     * If it is true the label cells.
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
     * Initialise the cell-cycle model at the start of a simulation.
     *
     * This overridden method sets up a new ODE system.
     */
    void Initialise();

    /**
     * Outputs cell cycle model parameters to files.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    virtual void OutputCellCycleModelParameters(out_stream& rParamsFile);
};

// Declare identifier for the serializer
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(Owen2011OxygenBasedCellCycleModel)
#include "CellCycleModelOdeSolverExportWrapper.hpp"
EXPORT_CELL_CYCLE_MODEL_ODE_SOLVER(Owen2011OxygenBasedCellCycleModel)

#endif /*OWEN2011OXYGENBASEDCELLCYCLEMODEL_HPP_*/
