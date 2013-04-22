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

#include "OnLatticeSimulationWithVasculature.hpp"
#include "CellBasedPdeHandler.hpp"
#include "StalkCellMutationState.hpp"
#include "TipCellMutationState.hpp"
#include "VesselCellMutationState.hpp"
#include "Debug.hpp"

template<unsigned DIM>
OnLatticeSimulationWithVasculature<DIM>::OnLatticeSimulationWithVasculature(AbstractCellPopulation<DIM>& rCellPopulation,
                                           bool deleteCellPopulationInDestructor,
                                           bool initialiseCells)
	: OnLatticeSimulation<DIM>(rCellPopulation, deleteCellPopulationInDestructor, initialiseCells)
{

}


template <unsigned DIM>
void OnLatticeSimulationWithVasculature<DIM>::UpdateAtEndOfTimeStep()
{
	// Find if any tip cells and vessel cells are in the same Lattice site.
	bool anastomosis_present = true;

	while(anastomosis_present)
	{
		bool is_iterator_valid = true;

		// Loop over cells to check Anastomosis
		for (typename AbstractCellPopulation<DIM>::Iterator tip_cell_iter = this->mrCellPopulation.Begin();
						tip_cell_iter != this->mrCellPopulation.End();
						 ++tip_cell_iter)
		{
			assert(is_iterator_valid);

			CellPtr p_tip_cell = (*tip_cell_iter);

			if (p_tip_cell->GetMutationState()->IsType<TipCellMutationState>() )
			{
				unsigned location_index = this->mrCellPopulation.GetLocationIndexUsingCell(p_tip_cell);

				std::set<CellPtr> all_cells_at_site =  this->mrCellPopulation.GetCellsUsingLocationIndex(location_index);

				assert(all_cells_at_site.size() > 0);

				if (all_cells_at_site.size() > 1)
				{
					for (std::set<CellPtr>::iterator cell_iter = all_cells_at_site.begin();
						 cell_iter != all_cells_at_site.end();
						 cell_iter++)
					{
						if ((*cell_iter) != p_tip_cell)
						{
							CellPtr p_cell = *cell_iter;

							if (p_cell->GetMutationState()->IsType<VesselCellMutationState>())
							{
								double protovessel_id = (p_tip_cell->GetCellData()->GetItem("StalkId"));

								std::vector<CellPtr> protovessel;

								assert(is_iterator_valid);
								for (typename AbstractCellPopulation<DIM>::Iterator stalk_cell_iter = this->mrCellPopulation.Begin();
									 stalk_cell_iter != this->mrCellPopulation.End();
									 ++stalk_cell_iter)
								{
									assert(is_iterator_valid);

									CellPtr p_stalk_cell = (*stalk_cell_iter);

									 if (p_stalk_cell->GetMutationState()->IsType<StalkCellMutationState>())
									 {
										 if (p_stalk_cell->GetCellData()->GetItem("StalkId") == protovessel_id)
										 {
											 protovessel.push_back(p_stalk_cell);
										 }
									 }
								}

								// Want to remove the tip cell. Note dont use Kill() as this would count the cell as dying rather than forming a vessel
	//							unsigned node_index = this->mrCellPopulation.GetLocationIndexUsingCell(p_tip_cell);
	//							this->mrCellPopulation.RemoveCellUsingLocationIndex(node_index, p_tip_cell);
								p_tip_cell->Kill();
								//p_tip_cell->SetMutationState(CellPropertyRegistry::Instance()->Get<VesselCellMutationState>());

								unsigned protovessel_length = protovessel.size();
								unsigned min_protovessel_length = 5u;

								for (unsigned i = 0; i < protovessel_length; i++)
								{
									if (protovessel_length < min_protovessel_length)
									{
										// Kill all stalk cells in protovessel
										protovessel[i]->Kill();
									}
									else
									{
										// Convert all stalk cells in protovessel to vessel cell
										protovessel[i]->SetMutationState(CellPropertyRegistry::Instance()->Get<VesselCellMutationState>());
										protovessel[i]->GetCellData()->SetItem("StalkId", 0);
									}
								}
								is_iterator_valid = false;
								break;
							}
						}
					}
				}
			}
			/*
			 *  Killing cells will invalidate the iterator so restart the cell iteration.
			 */
			if (!is_iterator_valid)
			{
				break;
			}
		}

		// No cell killing has happened. No Anastomosis has been found.
		if (is_iterator_valid)
		{
			anastomosis_present = false;
		}
	}


	// Need to call RemoveDeadCells now as otherwise it wouuldnt be called till next timestep and there would be dead cells lying around causing memory leaks.
	unsigned num_cells_killed = this->mrCellPopulation.RemoveDeadCells();
	if (num_cells_killed>0)
	{
		PRINT_VARIABLE(num_cells_killed);
	}

}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class OnLatticeSimulationWithVasculature<1>;
template class OnLatticeSimulationWithVasculature<2>;
template class OnLatticeSimulationWithVasculature<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(OnLatticeSimulationWithVasculature)
