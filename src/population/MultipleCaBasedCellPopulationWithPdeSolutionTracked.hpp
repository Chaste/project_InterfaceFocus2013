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

#ifndef MULTIPLECABASEDCELLPOPULATIONWITHPDESOLUTIONTRACKED_HPP_
#define MULTIPLECABASEDCELLPOPULATIONWITHPDESOLUTIONTRACKED_HPP_

#include "AbstractOnLatticeCellPopulation.hpp"
#include "PottsMesh.hpp"
#include "VertexMesh.hpp"
#include "AbstractMultipleCaUpdateRule.hpp"

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/vector.hpp>

// Needed here to avoid serialization errors (on Boost<1.37)
#include "WildTypeCellMutationState.hpp"
#include "MultipleCaBasedCellPopulation.hpp"


template<unsigned DIM>
class MultipleCaBasedCellPopulationWithPdeSolutionTracked : public MultipleCaBasedCellPopulation<DIM>
{
    friend class TestMultipleCaBasedCellPopulationWithPdeSolutionTracked;

private:
    /** Results file for cell locations. */
    //out_stream mpVizLocationsFile;

    /** Records for each node the Pde solution for oxygen. */
    std::vector<double> mOxygenPdeSolutionsAtNodes;

    /** Records for each node the Pde solution for VEGF. */
    std::vector<double> mVEGFPdeSolutionsAtNodes;

    friend class boost::serialization::access;
    /**
     * Serialize the object and its member variables.
     *
     * Note that serialization of the mesh and cells is handled by load/save_construct_data.
     *
     * Note also that member data related to writers is not saved - output must
     * be set up again by the caller after a restart.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
#define COVERAGE_IGNORE
        archive & boost::serialization::base_object<MultipleCaBasedCellPopulation<DIM> >(*this);
        archive & mOxygenPdeSolutionsAtNodes;
        archive & mVEGFPdeSolutionsAtNodes;
 #undef COVERAGE_IGNORE
    }

    /**
     * Overridden WriteVtkResultsToFile() method.
     */
    //void WriteVtkResultsToFile();

public:

    /**
     * Create a new cell population facade from a mesh and collection of cells.
     *
     * There must be precisely one CellPtr for each PottsElement in
     * the mesh.
     *
     * @param rMesh reference to a PottsMesh
     * @param rCells reference to a vector of CellPtrs
     * @param locationIndices an optional vector of location indices that correspond to real cells
     * @param latticeCarryingCapacity an optional parameter to allow more than one cell per site
     * @param deleteMesh set to true if you want the cell population to free the mesh memory on destruction
     *                   (defaults to false)
     * @param validate whether to validate the cell population when it is created (defaults to false as not used in CA simulations)
     */
    MultipleCaBasedCellPopulationWithPdeSolutionTracked(PottsMesh<DIM>& rMesh,
                                  std::vector<CellPtr>& rCells,
                                  const std::vector<unsigned> locationIndices=std::vector<unsigned>(),
                                  unsigned latticeCarryingCapacity=1u,
                                  bool deleteMesh=false,
                                  bool validate=false);

    /**
     * Constructor for use by the de-serializer.
     *
     * @param rMesh a vertex mesh.
     */
    MultipleCaBasedCellPopulationWithPdeSolutionTracked(PottsMesh<DIM>& rMesh);

    /**
     * Destructor, which frees any memory allocated by the constructor.
     */
    virtual ~MultipleCaBasedCellPopulationWithPdeSolutionTracked();

    /**
     * Sets the vector mOxygenPdeSolutionsAtNodes with the PDE results for the parameter oxygen.
     * This method sets all the values at once.
     *
     * @param rOxygenPdeSolutions is a vector passing the solution of the PDEs on the nodes
     */
    void SetOxygenPdeSolutions(const std::vector<double> rOxygenPdeSolutions);

    /**
     * Sets the vector mVEGFPdeSolutionsAtNodes with the PDE results for the parameter VEGF.
     * This method sets all the values at once.
     *
     * @param rVEGFPdeSolutions is a vector passing the solution of the PDEs on the nodes
     */
    void SetVEGFPdeSolutions(const std::vector<double> rVEGFPdeSolutions);

    /**
     * Gets the vector mOxygenPdeSolutionsAtNodes with the PDE results for the parameter Oxygen.
     */
    std::vector<double>& rGetOxygenPdeSolutions();

    /**
     * Gets the vector mVEGFPdeSolutionsAtNodes with the PDE results for the parameter VEGF.
     */
    std::vector<double>& rGetVEGFPdeSolutions();

    /**
     * Sets the vector mOxygenPdeSolutionsAtNodes with the PDE results for oxygen.
     * This method sets the value at a specific node.
     *
     * @param index defines the node index which value will be set.
     * @param solution defined the value to be set for the pde solution at the node
     */
    void SetOxygenPdeSolutionsAtNode(unsigned index, double solution);

    /**
      * Sets the vector mOxygenPdeSolutionsAtNodes with the PDE results for VEGF.
      * This method sets the value at a specific node.
      *
      * @param index defines the node index which value will be set.
      * @param solution defined the value to be set for the pde solution at the node
      */
    void SetVEGFPdeSolutionsAtNode(unsigned index, double solution);

    /**
     * Gets the PDE result for Oxygen at a certain node defined by index.
     *
     * @@param index defined the node index which value will be returned.
     */
    double GetOxygenPdeSolutionsAtNode(unsigned index);

    /**
      * Gets the PDE result for VEGF at a certain node defined by index.
      *
      * @@param index defined the node index which value will be returned.
      */
    double GetVEGFPdeSolutionsAtNode(unsigned index);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(MultipleCaBasedCellPopulationWithPdeSolutionTracked)

// No archiving yet so untested
#define COVERAGE_IGNORE
namespace boost
{
namespace serialization
{

// Serialize information required to construct a MultipleCaBasedCellPopulationWithPdeSolutionTracked.

template<class Archive, unsigned DIM>
inline void save_construct_data(
    Archive & ar, const MultipleCaBasedCellPopulationWithPdeSolutionTracked<DIM> * t, const BOOST_PFTO unsigned int file_version)
{
    // Save data required to construct instance
    const PottsMesh<DIM>* p_mesh = &(t->rGetMesh());
    ar & p_mesh;
}

/**
 * De-serialize constructor parameters and initialise a MultipleCaBasedCellPopulationWithPdeSolutionTracked.
 * Loads the mesh from separate files.
 */
template<class Archive, unsigned DIM>
inline void load_construct_data(
    Archive & ar, MultipleCaBasedCellPopulationWithPdeSolutionTracked<DIM> * t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    PottsMesh<DIM>* p_mesh;
    ar >> p_mesh;

    // Invoke inplace constructor to initialise instance
    ::new(t)MultipleCaBasedCellPopulationWithPdeSolutionTracked<DIM>(*p_mesh);
}
}
} // namespace ...
#undef COVERAGE_IGNORE

#endif /*MULTIPLECABASEDCELLPOPULATIONWITHPDESOLUTIONTRACKED_HPP_*/
