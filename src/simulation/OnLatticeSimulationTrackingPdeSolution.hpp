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

#ifndef ONLATTICESIMULATIONTRACKINGPDESOLUTION_HPP_
#define ONLATTICESIMULATIONTRACKINGPDESOLUTION_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "OnLatticeSimulation.hpp"
#include "MultipleCaBasedCellPopulationWithPdeSolutionTracked.hpp"

template<unsigned DIM>
class OnLatticeSimulationTrackingPdeSolution : public OnLatticeSimulation<DIM>
{
private:

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Serialize the object and any member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellBasedSimulation<DIM> >(*this);
    }


public:

    /**
     * Constructor.
     *
     * @param rCellPopulation A cell population object
     * @param deleteCellPopulationInDestructor Whether to delete the cell population on destruction to
     *        free up memory (defaults to false)
     * @param initialiseCells Whether to initialise cells (defaults to true, set to false when loading
     * from an archive)
     */
    OnLatticeSimulationTrackingPdeSolution(AbstractCellPopulation<DIM>& rCellPopulation,
                                           bool deleteCellPopulationInDestructor=false,
                                           bool initialiseCells=true);


    /**
     *  This method loops over all lattice sites from a PottsMesh, gets the solution to PDE,
     *  and stores that in the MultipleCaBasedCellPopulationWithPDESolutionTracked class when the
     *  simulation starts. This is used to initialize the PDE vectors in the population
     */
    void SetupSolve();

    /**
     *  This method loops over all lattice sites from a PottsMesh, gets the solution to PDE,
     *  and stores that in the MultipleCaBasedCellPopulationWithPDESolutionTracked class.
     */
    void UpdateAtEndOfTimeStep();

};

// Serialization for Boost >= 1.36
#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(OnLatticeSimulationTrackingPdeSolution)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct an OnLatticeSimulation.
 */
template<class Archive, unsigned DIM>
inline void save_construct_data(
    Archive & ar, const OnLatticeSimulationTrackingPdeSolution<DIM> * t, const BOOST_PFTO unsigned int file_version)
{
    // Save data required to construct instance
    const AbstractCellPopulation<DIM>* p_cell_population = &(t->rGetCellPopulation());
    ar & p_cell_population;
}

/**
 * De-serialize constructor parameters and initialise an OnLatticeSimulationTrackingPdeSolution.
 */
template<class Archive, unsigned DIM>
inline void load_construct_data(
    Archive & ar, OnLatticeSimulationTrackingPdeSolution<DIM> * t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    AbstractCellPopulation<DIM>* p_cell_population;
    ar >> p_cell_population;

    // Invoke inplace constructor to initialise instance, last two variables set extra
    // member variables to be deleted as they are loaded from archive and to not initialise sells.
    ::new(t)OnLatticeSimulationTrackingPdeSolution<DIM>(*p_cell_population, true, false);
}
}
} // namespace

#endif /*ONLATTICESIMULATIONTRACKINGPDESOLUTION_HPP_*/
