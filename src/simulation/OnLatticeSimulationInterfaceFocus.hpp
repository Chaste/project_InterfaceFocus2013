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

#ifndef ONLATTICESIMULATIONINTERFACEFOCUS_HPP_
#define ONLATTICESIMULATIONINTERFACEFOCUS_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractCellBasedSimulationInterfaceFocus.hpp"
#include "AbstractPottsUpdateRule.hpp"
#include "AbstractMultipleCaUpdateRule.hpp"

/**
 * Run an on-lattice 2D or 3D cell-based simulation.
 *
 * The OnLatticeSimulation is constructed with a CellPopulation, which
 * updates the correspondence between each Cell and its spatial representation
 * and handles cell division (governed by the CellCycleModel associated
 * with each cell). Once constructed, one or more Update rules may be passed
 * to the OnLatticeSimulation object, to define the processes which update
 * cells in the CellPopulation. Similarly, one or more CellKillers may be passed
 * to the OnLatticeSimulation object to specify conditions in which Cells
 * may die.
 */
template<unsigned DIM>
class OnLatticeSimulationInterfaceFocus : public AbstractCellBasedSimulationInterfaceFocus<DIM>
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
        archive & mOutputCellVelocities;
    }

protected:

    /**
     * Whether to write the cell velocities to a file.
     * Initialised to false in constuctor.
     */
    bool mOutputCellVelocities;

    /** Results file cell velocities. */
    out_stream mpCellVelocitiesFile;

    /**
     * Overridden UpdateCellPopulation() method.
     *
     * If using a CaBasedCellPopulation, this method does nothing if at the start of a simulation that
     * has just been loaded, to ensure consistency in random number generation.
     */
    void UpdateCellPopulation();

    /**
     * Overridden UpdateCellLocationsAndTopology() method.
     *
     * If using a PottsBasedCellPopulation, this method performs Monte Carlo sampling.
     */
    void UpdateCellLocationsAndTopology();

    /**
     * Overridden CalculateCellDivisionVector() method for determining how cell division occurs.
     * This method returns a vector which is then passed into the CellPopulation method AddCell().
     * This method may be overridden by subclasses.
     *
     * @param pParentCell the parent cell
     *
     * @return a vector containing information on cell division.
     */
    virtual c_vector<double, DIM> CalculateCellDivisionVector(CellPtr pParentCell);

    /**
     * Overridden WriteVisualizerSetupFile() method.
     */
    virtual void WriteVisualizerSetupFile();


public:

    /**
     * Constructor.
     *
     * @param rCellPopulation A cell population object
     * @param deleteCellPopulationInDestructor Whether to delete the cell population on destruction to
     *     free up memory (defaults to false)
     * @param initialiseCells Whether to initialise cells (defaults to true, set to false when loading
     * from an archive)
     */
    OnLatticeSimulationInterfaceFocus(AbstractCellPopulation<DIM>& rCellPopulation,
                                      bool deleteCellPopulationInDestructor=false,
                                      bool initialiseCells=true);

    /**
     * Add an update rule to be used in this simulation.
     *
     * @param pUpdateRule shared pointer to a multiple CA update rule law
     */
    void AddMultipleCaUpdateRule(boost::shared_ptr<AbstractMultipleCaUpdateRule<DIM> > pUpdateRule);

    /**
     * Method to remove all the MultipleCaUpdateRules
     */
    void RemoveAllMultipleCaUpdateRules();

    /**
     * Add an update rule to be used in this simulation (use this to set the Hamiltonian).
     *
     * @param pUpdateRule shared pointer to a Potts update rule law
     */
    void AddPottsUpdateRule(boost::shared_ptr<AbstractPottsUpdateRule<DIM> > pUpdateRule);

    /**
     * Method to remove all the PottsUpdateRules
     */
    void RemoveAllPottsUpdateRules();

    /**
     * Overridden OutputAdditionalSimulationSetup() method.
     * Outputs the update rule information.
     */
    bool GetOutputCellVelocities();

    /**
     * Set mOutputCellVelocities.
     *
     * @param outputCellVelocities the new value of mOutputCellVelocities
     */
    void SetOutputCellVelocities(bool outputCellVelocities);

    /**
     * Overridden SetupSolve() method to setup the cell velocities file.
     */
    virtual void SetupSolve();

    /**
     * Overridden UpdateAtEndOfSolve() method to close the cell velocities file.
     */
    virtual void UpdateAtEndOfSolve();

    /**
     * Overridden OutputAdditionalSimulationSetup() method to output the force and cell
     * population boundary condition information.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputAdditionalSimulationSetup(out_stream& rParamsFile);

    /**
     * Overridden OutputSimulationParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputSimulationParameters(out_stream& rParamsFile);
};

// Serialization for Boost >= 1.36
#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(OnLatticeSimulationInterfaceFocus)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct an OnLatticeSimulation.
 */
template<class Archive, unsigned DIM>
inline void save_construct_data(
    Archive & ar, const OnLatticeSimulationInterfaceFocus<DIM> * t, const BOOST_PFTO unsigned int file_version)
{
    // Save data required to construct instance
    const AbstractCellPopulation<DIM>* p_cell_population = &(t->rGetCellPopulation());
    ar & p_cell_population;
}

/**
 * De-serialize constructor parameters and initialise an OnLatticeSimulation.
 */
template<class Archive, unsigned DIM>
inline void load_construct_data(
    Archive & ar, OnLatticeSimulationInterfaceFocus<DIM> * t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    AbstractCellPopulation<DIM>* p_cell_population;
    ar >> p_cell_population;

    // Invoke inplace constructor to initialise instance, last two variables set extra
    // member variables to be deleted as they are loaded from archive and to not initialise sells.
    ::new(t)OnLatticeSimulationInterfaceFocus<DIM>(*p_cell_population, true, false);
}
}
} // namespace

#endif /*ONLATTICESIMULATIONINTERFACEFOCUS_HPP_*/
