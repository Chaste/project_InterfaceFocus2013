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

#include "VasctumCellKiller.hpp"
#include "CancerCellMutationState.hpp"
#include "QuiescentCancerCellMutationState.hpp"
#include "MacrophageMutationState.hpp"
#include "Owen2011OxygenBasedCellCycleModel.hpp"
#include "Owen2011OxygenBasedCellCycleModelWithoutOde.hpp"
#include "ApoptoticCellProperty.hpp"
#include "PottsMesh.hpp"

template<unsigned DIM>
VasctumCellKiller<DIM>::VasctumCellKiller(AbstractCellPopulation<DIM>* pCellPopulation)
: AbstractCellKiller<DIM>(pCellPopulation)
{
  mrho=0.75;
  mlowvalue=1.5;
  mhighvalue=DBL_MAX;
  mMacrophageMeanLifeSpan=300; // 90 days=129600
}

template<unsigned DIM>
void VasctumCellKiller<DIM>::CheckAndLabelSingleCellForApoptosis(CellPtr pCell)
{

	/* A normal cell is marked for apoptosis if the oxygen concentration within its area is below a certain threshold.
	    The cell ratio is calculated by the formula normal_count/(normal_count+cancer_count)
	   The threshold is defined as:
	 	if(ratio>mrho)
		  threshold = mlowvalue;
	    else threshold =  mhighvalue;
	    where:
	      mrho=0.75;
          mlowvalue=1.5;
          mhighvalue=DBL_MAX;
	*/
	if (!pCell->GetMutationState()->IsType<CancerCellMutationState>() && !pCell->GetMutationState()->IsType<MacrophageMutationState>())
	{
	   double ratio = CalculateRatio(pCell);
	   double threshold_oxygen = ThresholdOxygenConcentration(ratio);

	   if (pCell->GetCellData()->GetItem("oxygen") < threshold_oxygen)
	   {
		   pCell->Kill();
	   }
	}

	else if (pCell->GetMutationState()->IsType<QuiescentCancerCellMutationState>())
	{
		if(pCell->HasCellProperty<CellLabel>())
		{
			pCell->Kill();
	    }
	}

	else if (pCell->GetMutationState()->IsType<MacrophageMutationState>())
    {
        if (pCell->GetAge() > mMacrophageMeanLifeSpan)
	    {
	       pCell->Kill();
        }
    }
}

template<unsigned DIM>
double VasctumCellKiller<DIM>::CalculateRatio(CellPtr pCell)
{
	double normal_count=1;
	double cancer_count=0;

	//Get location index of the cell
    unsigned location_index = this->mpCellPopulation->GetLocationIndexUsingCell(pCell);

    // Counting the number of normal cells and cancer cells in the lattice site if there is
    //  more than one cell

   	if (this->mpCellPopulation->GetCellsUsingLocationIndex(location_index).size() > 1u)
    {
    	//Get the set of the Cell Ptr for cell's lattice site.
   		std::set<CellPtr> cells_on_lattice = this->mpCellPopulation->GetCellsUsingLocationIndex(location_index);
   		for (std::set<CellPtr>::iterator iter = cells_on_lattice.begin(); iter != cells_on_lattice.end(); iter++)
    	{
    	   	if ((*iter)->GetMutationState()->IsType<CancerCellMutationState>() || (*iter)->GetMutationState()->IsType<QuiescentCancerCellMutationState>() )
    	   		cancer_count++;
    		else
    			normal_count++;
    	}
    }
    else
    {
    	//Get the set of neighbouring cell indices
    	std::set<unsigned> neighbouring_node_indices = static_cast<PottsMesh<DIM>& >((this->mpCellPopulation)->rGetMesh()).GetMooreNeighbouringNodeIndices(location_index);

   		//Iterate over the set of neighbouring sites to estimate the total number of normal cells and cancer cells surrounding the target cell.
   		for (std::set<unsigned>::iterator neighbour_iter = neighbouring_node_indices.begin();
		     neighbour_iter != neighbouring_node_indices.end();
   		     ++neighbour_iter)
   		{
			//Get the set of the Cell Ptr for the current site.
   			std::set<CellPtr> cells_on_lattice = this->mpCellPopulation->GetCellsUsingLocationIndex(*neighbour_iter);

   			for (std::set<CellPtr>::iterator iter = cells_on_lattice.begin();
   				 iter != cells_on_lattice.end();
   				 iter++)
			{
    			if ((*iter)->GetMutationState()->IsType<CancerCellMutationState>() || (*iter)->GetMutationState()->IsType<QuiescentCancerCellMutationState>() )
    				cancer_count++;
    			else
    				normal_count++;
    		}

   		}
    }
    return (normal_count/(normal_count+cancer_count));
}

template<unsigned DIM>
double VasctumCellKiller<DIM>::ThresholdOxygenConcentration(double ratio)
{
	if(ratio>mrho)
		return mlowvalue;
	return mhighvalue;
}

template<unsigned SPACE_DIM>
void VasctumCellKiller<SPACE_DIM>::CheckAndLabelCellsForApoptosisOrDeath()
{
    for (typename AbstractCellPopulation<SPACE_DIM>::Iterator cell_iter = this->mpCellPopulation->Begin();
        cell_iter != this->mpCellPopulation->End();
        ++cell_iter)
    {
        CheckAndLabelSingleCellForApoptosis(*cell_iter);
    }
}

template<unsigned DIM>
void VasctumCellKiller<DIM>::OutputCellKillerParameters(out_stream& rParamsFile)
{

    // Call method on direct parent class
    AbstractCellKiller<DIM>::OutputCellKillerParameters(rParamsFile);
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class VasctumCellKiller<1>;
template class VasctumCellKiller<2>;
template class VasctumCellKiller<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(VasctumCellKiller)
