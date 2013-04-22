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

#ifndef MULTIPLECABASEDCELLPOPULATIONWITHVASCULATURE_HPP_
#define MULTIPLECABASEDCELLPOPULATIONWITHVASCULATURE_HPP_

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
class MultipleCaBasedCellPopulationWithVasculature : public MultipleCaBasedCellPopulation<DIM>
{
    friend class TestMultipleCaBasedCellPopulationWithVasculature;

private:

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
 #undef COVERAGE_IGNORE
    }

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
    MultipleCaBasedCellPopulationWithVasculature(PottsMesh<DIM>& rMesh,
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
    MultipleCaBasedCellPopulationWithVasculature(PottsMesh<DIM>& rMesh);

    /**
     * Destructor, which frees any memory allocated by the constructor.
     */
    virtual ~MultipleCaBasedCellPopulationWithVasculature();

    /**
      * Overridden EvaluateDivisionPropensity method to implement chemotaxis of tip cells
      *
      * Calculate the propensity of a dividing into a given site.
      *
      * @param currentNodeIndex The index of the current node/lattice site
      * @param targetNodeIndex The index of the target node/lattice site
      * @param cell a pointer to the cell (needed if more than one cell per lattice site
      * @return The probability of the cell dividing from the current node to the target node
      */
     double virtual EvaluateDivisionPropensity(unsigned currentNodeIndex,
                                        unsigned targetNodeIndex,
                                        CellPtr pCell);

    /**
     * Overridden IsSiteAvailable method to implement propper tip movement.
     *
     * Find if a given node has space available.
     *
     * @param index the global index of a specified node
     * @param the cell wanting to divide into the lattice site (defaults to NULL)
     *
     * @return whether the node is an empty site
     */
    virtual bool IsSiteAvailable(unsigned index, CellPtr pCell);


};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(MultipleCaBasedCellPopulationWithVasculature)

// No archiving yet so untested
#define COVERAGE_IGNORE
namespace boost
{
namespace serialization
{

// Serialize information required to construct a MultipleCaBasedCellPopulationWithVasculature.

template<class Archive, unsigned DIM>
inline void save_construct_data(
    Archive & ar, const MultipleCaBasedCellPopulationWithVasculature<DIM> * t, const BOOST_PFTO unsigned int file_version)
{
    // Save data required to construct instance
    const PottsMesh<DIM>* p_mesh = &(t->rGetMesh());
    ar & p_mesh;
}

/**
 * De-serialize constructor parameters and initialise a MultipleCaBasedCellPopulationWithVasculature.
 * Loads the mesh from separate files.
 */
template<class Archive, unsigned DIM>
inline void load_construct_data(
    Archive & ar, MultipleCaBasedCellPopulationWithVasculature<DIM> * t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    PottsMesh<DIM>* p_mesh;
    ar >> p_mesh;

    // Invoke inplace constructor to initialise instance
    ::new(t)MultipleCaBasedCellPopulationWithVasculature<DIM>(*p_mesh);
}
}
} // namespace ...
#undef COVERAGE_IGNORE

#endif /*MULTIPLECABASEDCELLPOPULATIONWITHVASCULATURE_HPP_*/
