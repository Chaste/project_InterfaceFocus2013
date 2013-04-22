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

#include "Owen2011OxygenBasedCellCycleOdeSystem.hpp"
#include "CellwiseOdeSystemInformation.hpp"


Owen2011OxygenBasedCellCycleOdeSystem::Owen2011OxygenBasedCellCycleOdeSystem(double oxygenConcentration,
		                                                                     boost::shared_ptr<AbstractCellMutationState> pMutationState,
                                                                             std::vector<double> stateVariables)
    : AbstractOdeSystem(2),
      mOxygenConcentration(oxygenConcentration),
      mpMutationState(pMutationState)
{
    mpSystemInfo.reset(new CellwiseOdeSystemInformation<Owen2011OxygenBasedCellCycleOdeSystem>);

    /**
     * State variables
     *
     * 0. phi = phase of the cell cycle
     * 1. oxygen_concentration
     */

    // Parameter values are taken from the Owen et al. (2011) paper. Oxygen concentration is measured in mmHg.
    //It is assumed that C=0.045 corresponds to typical dimensional perivascular oxygen tension of 20mmHg.
        if (mpMutationState->IsType<CancerCellMutationState>())
        {
        	 mTmin = 1600;
        	 mC0 = (1.4*0.045)/20;
        }
        else
        {
        	mTmin=3000;
        	mC0=(3*0.045)/20;
        }

    // Cell-specific initial conditions
    SetDefaultInitialCondition(1, oxygenConcentration);


    if (stateVariables != std::vector<double>())
    {
        SetStateVariables(stateVariables);
    }
}

Owen2011OxygenBasedCellCycleOdeSystem::~Owen2011OxygenBasedCellCycleOdeSystem()
{
    // Do nothing
}

void Owen2011OxygenBasedCellCycleOdeSystem::EvaluateYDerivatives(double time, const std::vector<double>& rY, std::vector<double>& rDY)
{
    //double phi = rY[0];
    double oxygen_concentration = rY[1];
    /*
     * The variables are
     * 0. phi = phase of the cell cycle
     * 1. oxygen concentration
     */
       double dphi = oxygen_concentration/(mTmin*(mC0+oxygen_concentration));

    // Rescale time to be in hours
    rDY[0] =60*dphi;
    rDY[1]=0.0;
}

bool Owen2011OxygenBasedCellCycleOdeSystem::CalculateStoppingEvent(double time, const std::vector<double>& rY)
{
    return (rY[0] >= 1 );
}

template<>
void CellwiseOdeSystemInformation<Owen2011OxygenBasedCellCycleOdeSystem>::Initialise()
{
    this->mVariableNames.push_back("phase_of_the_cell_cycle");
    this->mVariableUnits.push_back("non_dim");
    this->mInitialConditions.push_back(0.0);

     ///\toDo Oxygen Concentration (C) is calculated through a PDE in the cell population. However, for the moment the PDE is not working.
    this->mVariableNames.push_back("O2");
    this->mVariableUnits.push_back("non_dim");
    this->mInitialConditions.push_back(NAN); // will be filled in later

    this->mInitialised = true;
}

void Owen2011OxygenBasedCellCycleOdeSystem::SetMutationState(boost::shared_ptr<AbstractCellMutationState> pMutationState)
{
    mpMutationState = pMutationState;
}

const boost::shared_ptr<AbstractCellMutationState> Owen2011OxygenBasedCellCycleOdeSystem::GetMutationState() const
{
    return mpMutationState;
}

double Owen2011OxygenBasedCellCycleOdeSystem::GetOxygenConcentration() const
{
    return mOxygenConcentration;
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(Owen2011OxygenBasedCellCycleOdeSystem)
