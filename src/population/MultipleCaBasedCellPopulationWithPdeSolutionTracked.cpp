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

#include "MultipleCaBasedCellPopulationWithPdeSolutionTracked.hpp"
#include "RandomNumberGenerator.hpp"

// Needed to convert mesh in order to write nodes to VTK (visualize as glyphs)
#include "VtkMeshWriter.hpp"
#include "NodesOnlyMesh.hpp"
#include "Exception.hpp"
#include "Debug.hpp"


template<unsigned DIM>
MultipleCaBasedCellPopulationWithPdeSolutionTracked<DIM>::MultipleCaBasedCellPopulationWithPdeSolutionTracked(PottsMesh<DIM>& rMesh,
                                                        std::vector<CellPtr>& rCells,
                                                        const std::vector<unsigned> locationIndices,
                                                        unsigned latticeCarryingCapacity,
                                                        bool deleteMesh,
                                                        bool validate)
    : MultipleCaBasedCellPopulation<DIM>(rMesh, rCells, locationIndices, latticeCarryingCapacity, deleteMesh, validate)
{
    if (validate)
    {
        EXCEPTION("There is no validation for MultipleCaBasedCellPopulationWithPdes.");
    }
}

template<unsigned DIM>
MultipleCaBasedCellPopulationWithPdeSolutionTracked<DIM>::MultipleCaBasedCellPopulationWithPdeSolutionTracked(PottsMesh<DIM>& rMesh)
    : MultipleCaBasedCellPopulation<DIM>(rMesh)
{
}

template<unsigned DIM>
MultipleCaBasedCellPopulationWithPdeSolutionTracked<DIM>::~MultipleCaBasedCellPopulationWithPdeSolutionTracked()
{
    if (this->mDeleteMesh)
    {
        delete &this->mrMesh;
    }
}

template<unsigned DIM>
void MultipleCaBasedCellPopulationWithPdeSolutionTracked<DIM>::SetOxygenPdeSolutions(const std::vector<double> rOxygenPdeSolutions)
{
	mOxygenPdeSolutionsAtNodes = rOxygenPdeSolutions;
}

template<unsigned DIM>
void MultipleCaBasedCellPopulationWithPdeSolutionTracked<DIM>::SetVEGFPdeSolutions(const std::vector<double> rVEGFPdeSolutions)
{
	mVEGFPdeSolutionsAtNodes = rVEGFPdeSolutions;
}

template<unsigned DIM>
std::vector<double>& MultipleCaBasedCellPopulationWithPdeSolutionTracked<DIM>::rGetOxygenPdeSolutions()
{
   return mOxygenPdeSolutionsAtNodes;
}

template<unsigned DIM>
std::vector<double>& MultipleCaBasedCellPopulationWithPdeSolutionTracked<DIM>::rGetVEGFPdeSolutions()
{
   return mVEGFPdeSolutionsAtNodes;
}

template<unsigned DIM>
void MultipleCaBasedCellPopulationWithPdeSolutionTracked<DIM>::SetOxygenPdeSolutionsAtNode(unsigned index, double solution)
{
	mOxygenPdeSolutionsAtNodes[index] = solution;
}

template<unsigned DIM>
void MultipleCaBasedCellPopulationWithPdeSolutionTracked<DIM>::SetVEGFPdeSolutionsAtNode(unsigned index, double solution)
{
	mVEGFPdeSolutionsAtNodes[index] = solution;
}

template<unsigned DIM>
double MultipleCaBasedCellPopulationWithPdeSolutionTracked<DIM>::GetOxygenPdeSolutionsAtNode(unsigned index)
{
	return mOxygenPdeSolutionsAtNodes[index];
}

template<unsigned DIM>
double MultipleCaBasedCellPopulationWithPdeSolutionTracked<DIM>::GetVEGFPdeSolutionsAtNode(unsigned index)
{
	return mVEGFPdeSolutionsAtNodes[index];
}

//template<unsigned DIM>
//void MultipleCaBasedCellPopulationWithPdeSolutionTracked<DIM>::WriteResultsToFiles()
//{
//    AbstractCellPopulation<DIM>::WriteResultsToFiles();
//
//    SimulationTime* p_time = SimulationTime::Instance();
//
//    // Write location data to file
//    *mpVizLocationsFile << p_time->GetTime() << "\t";
//
//    // Loop over cells and find associated nodes so in the same order as the cells in output files
//    for (std::list<CellPtr>::iterator cell_iter = this->mCells.begin();
//         cell_iter != this->mCells.end();
//         ++cell_iter)
//    {
//        unsigned node_index = this->GetLocationIndexUsingCell(*cell_iter);
//
//        // Write the index of the of Node the cell is associated with.
//        *mpVizLocationsFile << node_index << " ";
//    }
//    *mpVizLocationsFile << "\n";
//}

//template<unsigned DIM>
//void MultipleCaBasedCellPopulationWithPdeSolutionTracked<DIM>::WriteVtkResultsToFile()
//{
//#ifdef CHASTE_VTK
//    ///\todo #2032 Compare with MeshBasedCellPopulation::WriteVtkResultsToFile etc.
//    std::stringstream time;
//    time << SimulationTime::Instance()->GetTimeStepsElapsed();
//    VtkMeshWriter<DIM, DIM> mesh_writer(this->mDirPath, "results_"+time.str(), false);
//
//    unsigned num_cells = this->GetNumRealCells();
//    std::vector<double> cell_types(num_cells, -1.0);
//    std::vector<double> cell_mutation_states(num_cells, -1.0);
//    std::vector<double> cell_labels(num_cells, -1.0);
//    std::vector<double> cell_ids(num_cells, -1.0);
//    std::vector<double> cell_ancestors(num_cells, -1.0);
//    std::vector<double> cell_ages(num_cells, -1.0);
//    std::vector<double> cell_cycle_phases(num_cells, -1.0);
//    std::vector<Node<DIM>*> nodes;
//
//    unsigned cell = 0;
//
//    // Counter to keep track of how many cells are at a lattice site
//    unsigned num_sites = this->mrMesh.GetNumNodes();
//    unsigned number_of_cells_at_site[num_sites];
//    for (unsigned i=0; i<num_sites; i++)
//    {
//    	number_of_cells_at_site[i] = 0;
//    }
//
//    for (std::list<CellPtr>::iterator iter = this->mCells.begin();
//         iter != this->mCells.end();
//         ++iter)
//    {
//        CellPtr cell_ptr = *iter;
//        cell_ids[cell] = cell_ptr->GetCellId();
//
//        unsigned location_index = this->GetLocationIndexUsingCell(*iter);
//
//        number_of_cells_at_site[location_index]++;
//        assert(number_of_cells_at_site[location_index]<=mLatticeCarryingCapacity);
//
//		c_vector<double, DIM> coords = this->mrMesh.GetNode(location_index)->rGetLocation();
//
//        // Move the coordinate slightly so that we can visualise all cells in a lattice site if there is more than one per site
//        if (mLatticeCarryingCapacity > 1)
//        {
//        	c_vector<double, DIM> offset;
//
//        	if (DIM == 2)
//        	{
//        		double angle = (double)number_of_cells_at_site[location_index]*2.0*M_PI/(double)mLatticeCarryingCapacity;
//        		offset[0] = 0.2*sin(angle);
//        		offset[1] = 0.2*cos(angle);
//        	}
//        	else
//        	{
//                RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
//
//        		for (unsigned i=0; i<DIM; i++)
//				{
//					offset[i] = p_gen->ranf(); // This assumes that all sites are 1 apart
//				}
//        	}
//
//    		for (unsigned i=0; i<DIM; i++)
//			{
//				coords[i] += offset[i];
//			}
//
//        }
//
//        nodes.push_back(new Node<DIM>(cell, coords, false));
//
//        if (this->mOutputCellAncestors)
//        {
//            double ancestor_index = (cell_ptr->GetAncestor() == UNSIGNED_UNSET) ? (-1.0) : (double)cell_ptr->GetAncestor();
//            cell_ancestors[cell] = ancestor_index;
//        }
//        if (this->mOutputCellProliferativeTypes)
//        {
//            cell_types[cell] = cell_ptr->GetCellProliferativeType()->GetColour();
//        }
//        if (this->mOutputCellMutationStates)
//        {
//            cell_mutation_states[cell] = cell_ptr->GetMutationState()->GetColour();
//        }
//        if (this->mOutputCellAges)
//        {
//            cell_ages[cell] = cell_ptr->GetAge();
//        }
//        if (this->mOutputCellCyclePhases)
//        {
//            cell_cycle_phases[cell] = cell_ptr->GetCellCycleModel()->GetCurrentCellCyclePhase();
//        }
//
//        cell ++;
//
//        ///\todo #2032 Add CellData
//    }
//
//    // Cell IDs can be used to threshold out the empty lattice sites (which have ID=-1)
//    mesh_writer.AddPointData("Cell ids", cell_ids);
//
//    if (this->mOutputCellProliferativeTypes)
//    {
//        mesh_writer.AddPointData("Cell types", cell_types);
//    }
//    if (this->mOutputCellAncestors)
//    {
//        mesh_writer.AddPointData("Ancestors", cell_ancestors);
//    }
//    if (this->mOutputCellMutationStates)
//    {
//        mesh_writer.AddPointData("Mutation states", cell_mutation_states);
//    }
//    if (this->mOutputCellAges)
//    {
//        mesh_writer.AddPointData("Ages", cell_ages);
//    }
//    if (this->mOutputCellCyclePhases)
//    {
//        mesh_writer.AddPointData("Cycle phases", cell_cycle_phases);
//    }
//    if (this->mOutputCellMutationStates)
//    {
//        mesh_writer.AddPointData("Mutation states", cell_mutation_states);
//        mesh_writer.AddPointData("Cell labels", cell_labels);
//    }
//
//    /*
//     * The current VTK writer can only write things which inherit from AbstractTetrahedralMeshWriter.
//     * For now, we do an explicit conversion to NodesOnlyMesh. This can be written to VTK then visualized as glyphs.
//     */
//    NodesOnlyMesh<DIM> temp_mesh;
//    temp_mesh.ConstructNodesWithoutMesh(nodes);
//    mesh_writer.WriteFilesUsingMesh(temp_mesh);
//
//    *(this->mpVtkMetaFile) << "        <DataSet timestep=\"";
//    *(this->mpVtkMetaFile) << time.str();
//    *(this->mpVtkMetaFile) << "\" group=\"\" part=\"0\" file=\"results_";
//    *(this->mpVtkMetaFile) << time.str();
//    *(this->mpVtkMetaFile) << ".vtu\"/>\n";
//
//    // Tidy up
//    for (unsigned i=0; i<nodes.size(); i++)
//    {
//        delete nodes[i];
//    }
//#endif //CHASTE_VTK
//}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class MultipleCaBasedCellPopulationWithPdeSolutionTracked<1>;
template class MultipleCaBasedCellPopulationWithPdeSolutionTracked<2>;
template class MultipleCaBasedCellPopulationWithPdeSolutionTracked<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(MultipleCaBasedCellPopulationWithPdeSolutionTracked)
