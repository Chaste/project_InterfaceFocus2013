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

#ifndef OWEN2011MULTIPLECAUPDATERULE_HPP_
#define OWEN2011MULTIPLECAUPDATERULE_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractMultipleCaUpdateRule.hpp"
#include "MultipleCaBasedCellPopulationWithPdeSolutionTracked.hpp"
#include "MultipleCaBasedCellPopulation.hpp"

// Needed here to avoid serialization errors (on Boost<1.37)
#include "CellLabel.hpp"
#include "PottsMesh.hpp"

/**
 * A diffusion update rule for use in cell-based simulations
 * using the cellular MultipleCa model and mock chemotaxis.
 */
template<unsigned DIM>
class Owen2011MultipleCaUpdateRule : public AbstractMultipleCaUpdateRule<DIM>
{
friend class TestMultipleCaUpdateRules;

private:

    /**
     * Diffusion parameter for normal and cancer cells.
     * Set to the default value 0.01 in the constructor.
     * \todo provide units
     */
    double mDiffusionParameter;

    /**
     * Diffusion parameter for macrophages.
     * Set to the default value 1.0 in the constructor.
     * \todo provide units
     */
    double mDiffusionParameterMacrophages;

    double mDt;
    /**
     * Chemotaxis coefficient for normal and cancer cells.
     * Set to the default value 0.0 in the constructor.
     * \todo provide units
     */
    double mChi;

    /**
     * Chemotaxis coefficient for macrophages.
     * Set to the default value 2e-4 in the constructor.
     * \todo provide units
     */
    double mChiMacrophages;

    double mA;
    double mB;

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractMultipleCaUpdateRule<DIM> >(*this);
        archive & mDiffusionParameter;
        archive & mDiffusionParameterMacrophages;
        archive & mDt;
        archive & mChi;
        archive & mChiMacrophages;
    }

public:

    /**
     * Constructor.
     */
    Owen2011MultipleCaUpdateRule();

    /**
     * Destructor.
     */
    ~Owen2011MultipleCaUpdateRule();

    /**
     * Calculate the probability of a given move.
     *
     * Uses random diffusion to each neighbouring node, scaled according to distance.
     *
     * @param currentNodeIndex The index of the current node/lattice site
     * @param targetNodeIndex The index of the target node/lattice site
     * @param rCellPopulation The cell population
     * @param dt is the time interval
     * @param deltaX defines the size of the lattice site
     * @param cell defines the current cell so that the mutation state can be determined.
     * The use of this parameter is necessary when there are multiple cells per lattice sites
     * @return The probability of the cell moving from the current node to the target node
     */
    double EvaluateProbability(unsigned currentNodeIndex,
                               unsigned targetNodeIndex,
                               MultipleCaBasedCellPopulation<DIM>& rCellPopulation,
                               double dt,
                               double deltaX,
                               CellPtr pCell);
    /**
    * Defined to test simulations with chemotaxis. Should be removed once PDE starts working.
    *
    */
    double MockChemotaxis(unsigned currentNodeIndex,
                          unsigned targetNodeIndex,
                          MultipleCaBasedCellPopulation<DIM>& rCellPopulation,
                          CellPtr pCell);

    /**
      * @return mDiffusionParameter
      */
    double GetDiffusionParameter();

    /**
     * Set mDiffusionParameter.
     *
     * @param diffusionParameter the new value of mDiffusionParameter
     */
    void SetDiffusionParameter(double diffusionParameter);


    /**
     * Set mChi.
     *
     * @param chi the new value of mChi
     */
    void SetVEGFParameterChi(double chi);

    /**
     * @return mDiffusionParameterMacrophages
     */
    double GetDiffusionParameterForMacrophages();

    /**
     * Set mDiffusionParameterMacrophages.
     *
     * @param diffusionParameter the new value of mDiffusionParameterMacrophages
     */
    void SetDiffusionParameterForMacrophages(double diffusionParameter);


    /**
     * Set mChiMacrophages.
     *
     * @param chi the new value of mChiMacrophages
     */
    void SetVEGFParameterChiForMacrophages(double chi);

    /**
     * Set mA.
     * This should be removed once the PDE module starts working for the multiple CA population
     * @param A used to calculate the mock VEGF chemotaxis
     */
    void SetVEGFParameterA(double a);

    /**
      * Set mB.
      * This should be removed once the PDE module starts working for the multiple CA population
      * @param B used to calculate the mock VEGF chemotaxis
      */
     void SetVEGFParameterB(double B);

//    /**
//      * @return the time step (mDt) defined in the calculation of the probability of a cell moving to a certain site
//      */
//    double GetTimeStep();
//
//
//    /**
//      * @define the time step (dt) defined in the calculation of the probability of a cell moving to a certain site
//      */
//    void SetTimeStep(double dt);

    /**
     * Overridden OutputUpdateRuleParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputUpdateRuleParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(Owen2011MultipleCaUpdateRule)

#endif /*OWEN2011MULTIPLECAUPDATERULE_HPP_*/
