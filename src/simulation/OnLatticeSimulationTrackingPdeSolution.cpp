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

#include "OnLatticeSimulationTrackingPdeSolution.hpp"
#include "CellBasedPdeHandler.hpp"

template<unsigned DIM>
OnLatticeSimulationTrackingPdeSolution<DIM>::OnLatticeSimulationTrackingPdeSolution(AbstractCellPopulation<DIM>& rCellPopulation,
                                           bool deleteCellPopulationInDestructor,
                                           bool initialiseCells)
	: OnLatticeSimulation<DIM>(rCellPopulation, deleteCellPopulationInDestructor, initialiseCells)
{

}

template <unsigned DIM>
void OnLatticeSimulationTrackingPdeSolution<DIM>::SetupSolve()
{
	UpdateAtEndOfTimeStep();
}

template <unsigned DIM>
void OnLatticeSimulationTrackingPdeSolution<DIM>::UpdateAtEndOfTimeStep()
{
	CellBasedPdeHandler<DIM> *pde_handler = this->GetCellBasedPdeHandler();

	if(static_cast<MultipleCaBasedCellPopulationWithPdeSolutionTracked<DIM>*>(&(this->mrCellPopulation)))
	{

		std::vector<double>VEGFSolutionAtNodes;
		std::vector<double>OxygenSolutionAtNodes;

		if (pde_handler != NULL)
		{
			// This casting is necessary because mrCellPopulation is of the Abstract type
			MultipleCaBasedCellPopulationWithPdeSolutionTracked<DIM>* p_static_cast_cell_population = static_cast<MultipleCaBasedCellPopulationWithPdeSolutionTracked<DIM>*>(&(this->mrCellPopulation));

			for(typename PottsMesh<DIM>::NodeIterator node_iter = p_static_cast_cell_population->rGetMesh().GetNodeIteratorBegin();
					 node_iter != p_static_cast_cell_population->rGetMesh().GetNodeIteratorEnd();
					 ++node_iter)
			{

				c_vector<double,DIM> r_location = node_iter->rGetLocation();

				VEGFSolutionAtNodes.push_back(pde_handler->GetPdeSolutionAtPoint(r_location, "vegf"));
				OxygenSolutionAtNodes.push_back(pde_handler->GetPdeSolutionAtPoint(r_location, "oxygen"));
			}

			// Setting the solutions in the population
			p_static_cast_cell_population->SetVEGFPdeSolutions(VEGFSolutionAtNodes);
			p_static_cast_cell_population->SetOxygenPdeSolutions(OxygenSolutionAtNodes);
		}
		else
		{
			EXCEPTION("A pde handler must be defined for the simulation with PDE solution tracked.");
		}
	}
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class OnLatticeSimulationTrackingPdeSolution<1>;
template class OnLatticeSimulationTrackingPdeSolution<2>;
template class OnLatticeSimulationTrackingPdeSolution<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(OnLatticeSimulationTrackingPdeSolution)
