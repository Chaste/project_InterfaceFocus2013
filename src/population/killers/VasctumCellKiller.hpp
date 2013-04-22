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

#ifndef VASCTUMCELLKILLER_HPP_
#define VASCTUMCELLKILLER_HPP_

#include "AbstractCellKiller.hpp"

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

/**
 * Cancer cells die if they are quiescent for a long time.
 * Normal cells die if the oxygen concentration falls below a certain threshold.
 * The threshold takes a higher when a normal cell is surrounded by more cancer cells
 * and a lower value otherwise.
 */
template<unsigned DIM>
class VasctumCellKiller : public AbstractCellKiller<DIM>
{
private:

	/**
	 * Threshold value of the ratio of normal cells to normal and cancer cells.
	 */
	double mrho;
	/**
     * Low threshold value of oxygen concentration.
	 */
	double mlowvalue;

	/**
	 * High threshold value of oxygen concentration.
	 */
	double mhighvalue;

    /**
     * Mean life span of the macrophage.
     * Has units of mins.
     */
    double mMacrophageMeanLifeSpan;

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the object.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellKiller<DIM> >(*this);
    }

public:

    /**
     * Default constructor.
     *
     * @param pCellPopulation pointer to the cell population
     */
    VasctumCellKiller(AbstractCellPopulation<DIM>* pCellPopulation);

    /**
     * Value returned by CalculateRatio() is then compared with 'mrho'
     * to decide the threshold oxygen concentration.
     */
    double ThresholdOxygenConcentration(double ratio);

    /*
     *
     */
    void CheckAndLabelSingleCellForApoptosis(CellPtr pCell);

    /**
     * Calculate the local ratio of normal cells to normal cells and cancer cells.
     * For this purpose we consider the neighbourhood of the cell's lattice site.
     * If the cell's lattice site contains more than one cell, then neighbourhood
     * is simply that lattice site and otherwise the Moore neighbourhood is considered.
     */
    double CalculateRatio(CellPtr pCell);

    /**
     * Loop over cells and kill according to the rule.
     */
    virtual void CheckAndLabelCellsForApoptosisOrDeath();

    /**
     * @return mMacrophageMeanLifeSpan
     */
    double GetMacrophageMeanLifeSpan();

    /**
     * Set method for mMacrophageMeanLifeSpan.
     *
     * @param macrophageMeanLifeSpan the new value of mMacrophageMeanLifeSpan
     */
    void SetMacrophageMeanLifeSpan(double macrophageMeanLifeSpan);

    /**
     * Overridden OutputCellKillerParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputCellKillerParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(VasctumCellKiller)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a RandomCellKiller.
 */
template<class Archive, unsigned DIM>
inline void save_construct_data(
    Archive & ar, const VasctumCellKiller<DIM> * t, const BOOST_PFTO unsigned int file_version)
{
    // Save data required to construct instance
    const AbstractCellPopulation<DIM>* const p_cell_population = t->GetCellPopulation();
    ar << p_cell_population;
}

/**
 * De-serialize constructor parameters and initialise a RandomCellKiller.
 */
template<class Archive, unsigned DIM>
inline void load_construct_data(
    Archive & ar, VasctumCellKiller<DIM> * t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    AbstractCellPopulation<DIM>* p_cell_population;
    ar >> p_cell_population;

    // Invoke inplace constructor to initialise instance
    ::new(t)VasctumCellKiller<DIM>(p_cell_population);
}
}
} // namespace ...

#endif /*VASCTUMCELLKILLER_HPP_*/
