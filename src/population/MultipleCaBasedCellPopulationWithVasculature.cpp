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

#include "MultipleCaBasedCellPopulationWithVasculature.hpp"
#include "RandomNumberGenerator.hpp"
#include "StalkCellMutationState.hpp"
#include "VesselCellMutationState.hpp"
#include "TipCellMutationState.hpp"

#include "Debug.hpp"

template<unsigned DIM>
MultipleCaBasedCellPopulationWithVasculature<DIM>::MultipleCaBasedCellPopulationWithVasculature(PottsMesh<DIM>& rMesh,
                                                        std::vector<CellPtr>& rCells,
                                                        const std::vector<unsigned> locationIndices,
                                                        unsigned latticeCarryingCapacity,
                                                        bool deleteMesh,
                                                        bool validate)
    : MultipleCaBasedCellPopulation<DIM>(rMesh, rCells, locationIndices, latticeCarryingCapacity, deleteMesh, validate)
{
}

template<unsigned DIM>
MultipleCaBasedCellPopulationWithVasculature<DIM>::MultipleCaBasedCellPopulationWithVasculature(PottsMesh<DIM>& rMesh)
    : MultipleCaBasedCellPopulation<DIM>(rMesh)
{
}

template<unsigned DIM>
MultipleCaBasedCellPopulationWithVasculature<DIM>::~MultipleCaBasedCellPopulationWithVasculature()
{
    if (this->mDeleteMesh)
    {
        delete &this->mrMesh;
    }
}

template<unsigned DIM>
double MultipleCaBasedCellPopulationWithVasculature<DIM>:: EvaluateDivisionPropensity(unsigned currentNodeIndex,
    																				   unsigned targetNodeIndex,
																					   CellPtr pCell)
{
	double propensity = 1.0;

	if (pCell->GetMutationState()->IsType<StalkCellMutationState>())
	{
		c_vector<double, DIM> node_index_location = this->GetNode(currentNodeIndex)->rGetLocation();
		c_vector<double, DIM> node_neighbour_location = this->GetNode(targetNodeIndex)->rGetLocation();

		if (node_index_location(1) < node_neighbour_location(1))
		{
			propensity = 1.0;
		}
		else
		{
			propensity = 0.1;
		}
	}

	return propensity;
}


template<unsigned DIM>
bool MultipleCaBasedCellPopulationWithVasculature<DIM>::IsSiteAvailable(unsigned index, CellPtr pCell)
{
	bool space_availiable = (this->rGetAvailableSpaces()[index] != 0);

	if (space_availiable)
	{
		if (pCell->GetMutationState()->IsType<StalkCellMutationState>())
		{
			std::set<CellPtr> occupying_cells =  this->GetCellsUsingLocationIndex(index);

			if (occupying_cells.size() > 0)
			{
				for (std::set<CellPtr>::iterator iter = occupying_cells.begin();
					iter != occupying_cells.end();
					iter++)
				{
					if (pCell != (*iter))
					{
						// check if room for new tip cell
						//if ((*iter)->GetCellData()->GetItem("StalkId")>0)
						if ((*iter)->GetMutationState()->IsType<StalkCellMutationState>() )
						{
							return false;
						}
						else if ((*iter)->GetMutationState()->IsType<TipCellMutationState>() )
						{
							return false;
						}
						else if ( (*iter)->GetMutationState()->IsType<VesselCellMutationState>() )
						{
							return true;
						}

					}
				}
			}
		}
	}
	return space_availiable;
}


/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class MultipleCaBasedCellPopulationWithVasculature<1>;
template class MultipleCaBasedCellPopulationWithVasculature<2>;
template class MultipleCaBasedCellPopulationWithVasculature<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(MultipleCaBasedCellPopulationWithVasculature)
