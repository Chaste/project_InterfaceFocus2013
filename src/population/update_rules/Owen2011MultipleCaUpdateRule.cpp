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

#include "Owen2011MultipleCaUpdateRule.hpp"

#include "WildTypeCellMutationState.hpp"
#include "CancerCellMutationState.hpp"
#include "MacrophageMutationState.hpp"
#include "VesselCellMutationState.hpp"
#include "StalkCellMutationState.hpp"
#include "TipCellMutationState.hpp"
#include "QuiescentCancerCellMutationState.hpp"


template<unsigned DIM>
Owen2011MultipleCaUpdateRule<DIM>::Owen2011MultipleCaUpdateRule()
    : AbstractMultipleCaUpdateRule<DIM>(),
      mDiffusionParameter(0.01),
      mDiffusionParameterMacrophages(1.0),
      mChi(0),
      mChiMacrophages(2e-4),
      mA(0.5),
      mB(0.3)
{
}

template<unsigned DIM>
Owen2011MultipleCaUpdateRule<DIM>::~Owen2011MultipleCaUpdateRule()
{
}

template<unsigned DIM>
double Owen2011MultipleCaUpdateRule<DIM>::EvaluateProbability(unsigned currentNodeIndex,
                                                              unsigned targetNodeIndex,
                                                              MultipleCaBasedCellPopulation<DIM>& rCellPopulation,
                                                              double dt,
                                                              double deltaX,
                                                              CellPtr pCell)
{
	if ( pCell->GetMutationState()->IsType<VesselCellMutationState>() ||
	     pCell->GetMutationState()->IsType<TipCellMutationState>() ||
	     pCell->GetMutationState()->IsType<StalkCellMutationState>() )
	{
		return 0.0;
	}
	else if (pCell->GetMutationState()->IsType<MacrophageMutationState>())
	{
		c_vector<double, DIM> node_index_location = rCellPopulation.GetNode(currentNodeIndex)->rGetLocation();
		c_vector<double, DIM> node_neighbour_location = rCellPopulation.GetNode(targetNodeIndex)->rGetLocation();

		double diffusion = mDiffusionParameterMacrophages * dt/(2 * pow(norm_2(rCellPopulation.rGetMesh().GetVectorFromAtoB(node_index_location, node_neighbour_location)), 2));

		if(dynamic_cast<MultipleCaBasedCellPopulationWithPdeSolutionTracked<DIM>*>(&(rCellPopulation)))
		{
			MultipleCaBasedCellPopulationWithPdeSolutionTracked<DIM>* p_static_cast_cell_population = static_cast<MultipleCaBasedCellPopulationWithPdeSolutionTracked<DIM>*>(&(rCellPopulation));

			double vegfLevelAtCurrentNode = p_static_cast_cell_population->GetVEGFPdeSolutionsAtNode(currentNodeIndex);
			double vegfLevelAtNeighbourNode = p_static_cast_cell_population->GetVEGFPdeSolutionsAtNode(targetNodeIndex);

			return(diffusion * deltaX * (1 + mChi * (vegfLevelAtNeighbourNode - vegfLevelAtCurrentNode)/(2 * mDiffusionParameter)));
		}
		else // For now: if it is not tracking the PDE solution in the population, the mock chemotaxis will be used for macrophages.
		{
			double mock = MockChemotaxis(currentNodeIndex, targetNodeIndex, rCellPopulation, pCell);
			return(diffusion * deltaX * (1 + mock/(2 * mDiffusionParameterMacrophages)));
		}
	}
	else
	{
		assert(pCell->GetMutationState()->IsType<WildTypeCellMutationState>() ||
			   pCell->GetMutationState()->IsType<CancerCellMutationState>() ||
			   pCell->GetMutationState()->IsType<QuiescentCancerCellMutationState>());

		c_vector<double, DIM> node_index_location = rCellPopulation.GetNode(currentNodeIndex)->rGetLocation();
		c_vector<double, DIM> node_neighbour_location = rCellPopulation.GetNode(targetNodeIndex)->rGetLocation();

		// Same diffusion parameter for normal and cancer cells
		double diffusion = mDiffusionParameter*dt/(2* pow(norm_2(rCellPopulation.rGetMesh().GetVectorFromAtoB(node_index_location, node_neighbour_location)), 2));

		//if (rCellPopulation.GetCellUsingLocationIndex(currentNodeIndex)->GetMutationState()!=p_cancer_mutation)
		//{

			if(dynamic_cast<MultipleCaBasedCellPopulationWithPdeSolutionTracked<DIM>*>(&rCellPopulation))
			{
				MultipleCaBasedCellPopulationWithPdeSolutionTracked<DIM>* p_static_cast_cell_population = static_cast<MultipleCaBasedCellPopulationWithPdeSolutionTracked<DIM>*>(&(rCellPopulation));
				double vegfLevelAtCurrentNode = p_static_cast_cell_population->GetVEGFPdeSolutionsAtNode(currentNodeIndex);
				double vegfLevelAtNeighbourNode = p_static_cast_cell_population->GetVEGFPdeSolutionsAtNode(targetNodeIndex);

				//PRINT_VARIABLE(vegfLevelAtCurrentNode);
				//PRINT_VARIABLE(vegfLevelAtNeighbourNode);
				//PRINT_VARIABLE(diffusion);
				//PRINT_VARIABLE(diffusion * deltaX * (1 + mChi * (vegfLevelAtNeighbourNode - vegfLevelAtCurrentNode)/(2 * mDiffusionParameter)));

				double advection_diffusion  = diffusion * deltaX * (1 + mChi * (vegfLevelAtNeighbourNode - vegfLevelAtCurrentNode)/(2 * mDiffusionParameter));
				return advection_diffusion;
			}
			// For now: if it is not tracking the PDE solution in the population, the mock chemotaxis will be used for normal cells
			else
			{
				double mock = MockChemotaxis(currentNodeIndex, targetNodeIndex, rCellPopulation, pCell);
				double advection_diffusion  = diffusion * deltaX * (1 + mock/(2 * mDiffusionParameter));
				return advection_diffusion;

			}
		//}
		//else // cancer cells do not have chemotaxis
		//{
		//    return(diffusion);
		//}

	}
}

template<unsigned DIM>
double Owen2011MultipleCaUpdateRule<DIM>::MockChemotaxis(unsigned currentNodeIndex,
                                                          unsigned targetNodeIndex,
		                                                  MultipleCaBasedCellPopulation<DIM>& rCellPopulation,
		                                                  CellPtr pCell)
{
    c_vector<double, DIM> node_index_location = rCellPopulation.GetNode(currentNodeIndex)->rGetLocation();
	c_vector<double, DIM> node_neighbour_location = rCellPopulation.GetNode(targetNodeIndex)->rGetLocation();

	double VCurrentNode = mA * node_index_location[0] + mB;
    double VTargetNode = mA * node_neighbour_location[0]+ mB;

    if(pCell->GetMutationState()->IsType<MacrophageMutationState>())
    	return(mChiMacrophages *(VTargetNode - VCurrentNode));
    else
    	return(mChi *(VTargetNode - VCurrentNode));
}


template<unsigned DIM>
double Owen2011MultipleCaUpdateRule<DIM>::GetDiffusionParameter()
{
    return mDiffusionParameter;
}

template<unsigned DIM>
void Owen2011MultipleCaUpdateRule<DIM>::SetDiffusionParameter(double diffusionParameter)
{
	mDiffusionParameter = diffusionParameter;
}

template<unsigned DIM>
double Owen2011MultipleCaUpdateRule<DIM>::GetDiffusionParameterForMacrophages()
{
    return mDiffusionParameterMacrophages;
}

template<unsigned DIM>
void Owen2011MultipleCaUpdateRule<DIM>::SetDiffusionParameterForMacrophages(double diffusionParameter)
{
    mDiffusionParameterMacrophages = diffusionParameter;
}

template<unsigned DIM>
void Owen2011MultipleCaUpdateRule<DIM>::SetVEGFParameterChi(double chi)
{
	mChi = chi;
}

template<unsigned DIM>
void Owen2011MultipleCaUpdateRule<DIM>::SetVEGFParameterChiForMacrophages(double chi)
{
    mChiMacrophages = chi;
}

template<unsigned DIM>
void Owen2011MultipleCaUpdateRule<DIM>::SetVEGFParameterA(double a)
{
	mA = a;
}

template<unsigned DIM>
void Owen2011MultipleCaUpdateRule<DIM>::SetVEGFParameterB(double b)
{
	mB = b;
}

//template<unsigned DIM>
//double DiffusionMultipleCaUpdateRule<DIM>::GetTimeStep()
//{
//    return mDt;
//}
//
//template<unsigned DIM>
//void DiffusionMultipleCaUpdateRule<DIM>::SetTimeStep(double dt)
//{
//	mDt = dt;
//}

template<unsigned DIM>
void Owen2011MultipleCaUpdateRule<DIM>::OutputUpdateRuleParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<DiffusionParameter>" << mDiffusionParameter << "</DiffusionParameter>\n";
    *rParamsFile << "\t\t\t<DiffusionParameterForMacrophages>" << mDiffusionParameterMacrophages << "</DiffusionParameterForMacrophages>\n";
    *rParamsFile << "\t\t\t<ChemotaxisCoefficient>" << mChi << "</ChemotaxisCoefficient>\n";
    *rParamsFile << "\t\t\t<ChemotaxisCoefficientForMacrophagesr>" << mChiMacrophages << "</ChemotaxisCoefficientForMacrophages>\n";

    // Call method on direct parent class
    AbstractMultipleCaUpdateRule<DIM>::OutputUpdateRuleParameters(rParamsFile);
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class Owen2011MultipleCaUpdateRule<1u>;
template class Owen2011MultipleCaUpdateRule<2u>;
template class Owen2011MultipleCaUpdateRule<3u>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(Owen2011MultipleCaUpdateRule)
