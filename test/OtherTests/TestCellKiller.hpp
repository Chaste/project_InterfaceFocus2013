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

#ifndef TESTCELLKILLER_HPP_
#define TESTCELLKILLER_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "Owen2011OxygenBasedCellCycleModel.hpp"
#include "CellsGenerator.hpp"
#include "VasctumCellKiller.hpp"
#include "WildTypeCellMutationState.hpp"
#include "CancerCellMutationState.hpp"
#include "QuiescentCancerCellMutationState.hpp"
#include "ApoptoticCellProperty.hpp"
#include "CellPropertyRegistry.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "PottsBasedCellPopulation.hpp"
#include "AbstractCellPopulation.hpp"
#include "PottsMeshGenerator.hpp"
#include "CellsGenerator.hpp"
#include "MultipleCaBasedCellPopulation.hpp"
#include "OnLatticeSimulation.hpp"
#include "Warnings.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "SmartPointers.hpp"

/**
 * This class contains tests for methods on classes
 * inheriting from AbstractCellKiller.
 */
class TestCellKiller : public AbstractCellBasedTestSuite
{
public:

    void TestVasctumCellKiller() throw(Exception)
    {
    	// Set up singleton classes
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(5000, 5000);

    	// Create a simple 2D PottsMesh
        PottsMeshGenerator<2> generator(5, 0, 0, 5, 0, 0);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        // Create cells
		std::vector<CellPtr> cells;
		MAKE_PTR(StemCellProliferativeType, p_stem_type);
		CellsGenerator<Owen2011OxygenBasedCellCycleModel, 2> cells_generator;
		cells_generator.GenerateBasicRandom(cells, 15, p_stem_type);


       	// Create a cell mutation state
        boost::shared_ptr<AbstractCellProperty> p_cancer_mutation(CellPropertyRegistry::Instance()->Get<CancerCellMutationState>());
        for (unsigned i = 0; i<cells.size(); i++)
        {
            if (i==0 || i==1 || i==2 || i==4 || i==5 || i==6 || i==7)
            {
            	cells[i]->SetMutationState(p_cancer_mutation);
            }
        }

        // Set up oxygen_concentration
        double oxygen_concentration = 100.0;
        double hypoxic_concentration = 0.1;

        //Pass CellData to all cells
        for (unsigned i = 0; i<cells.size(); i++)
        {
        	if(i==6 || i==7 || i==10 || i==11 || i==12)
        	cells[i]->GetCellData()->SetItem("oxygen", hypoxic_concentration);
        	else
        	cells[i]->GetCellData()->SetItem("oxygen", oxygen_concentration);
        }

        // Specify where cells lie
        std::vector<unsigned> location_indices;
        location_indices.push_back(1u);
        location_indices.push_back(1u);
        location_indices.push_back(1u);
        location_indices.push_back(3u);
        location_indices.push_back(3u);
        location_indices.push_back(3u);
        location_indices.push_back(6u);
        location_indices.push_back(6u);
        location_indices.push_back(7u);
        location_indices.push_back(11u);
        location_indices.push_back(13u);
        location_indices.push_back(13u);
        location_indices.push_back(13u);
        location_indices.push_back(16u);
        location_indices.push_back(24u);

        // Create cell population
        MultipleCaBasedCellPopulation<2> cell_population(*p_mesh, cells, location_indices, 10);
        cell_population.SetOutputCellIdData(true);
        cell_population.SetOutputCellMutationStates(true);
        cell_population.SetOutputCellProliferativeTypes(true);
        cell_population.SetOutputCellCyclePhases(true);
        cell_population.SetOutputCellAncestors(true);
        cell_population.SetOutputCellAges(true);

        // Get a reference to the cells held in cell population
        std::list<CellPtr>& r_cells = cell_population.rGetCells();

        //Add Cells using the location index
        for(unsigned i=0; i<r_cells.size(); i++)
        cell_population.AddCellUsingLocationIndex(location_indices[i],cells[i]);

        // Create cell killer
        VasctumCellKiller<2> cell_killer(&cell_population);

        TS_ASSERT_EQUALS(cell_killer.GetIdentifier(), "VasctumCellKiller-2");

        TS_ASSERT_THROWS_NOTHING(cell_killer.CheckAndLabelSingleCellForApoptosis(*r_cells.begin()));

        //Test to check death of normal cells based on it's environment and oxygen concentration.

        //Surrounded by tumour microenvironment, hence must be dead
        cell_killer.CheckAndLabelSingleCellForApoptosis(*cell_population.GetCellsUsingLocationIndex(7u).begin());
        TS_ASSERT_EQUALS((*(cell_population.GetCellsUsingLocationIndex(7u)).begin())->IsDead(), true);

        //No cancer cells around but oxygen concentration is below hypoxic threshold, hence must be dead
        cell_killer.CheckAndLabelSingleCellForApoptosis(*cell_population.GetCellsUsingLocationIndex(13u).begin());
        TS_ASSERT_EQUALS((*(cell_population.GetCellsUsingLocationIndex(13u)).begin())->IsDead(), true);

        //Surrounded by normal cells and oxygen concentration above hypoxic threshold, hence must not be dead
        cell_killer.CheckAndLabelSingleCellForApoptosis(*cell_population.GetCellsUsingLocationIndex(16u).begin());
        TS_ASSERT_EQUALS((*(cell_population.GetCellsUsingLocationIndex(16u)).begin())->IsDead(), false);

        // Surrounded by no cells, oxygen concentration above hypoxic threshold, hence must not be dead
        cell_killer.CheckAndLabelSingleCellForApoptosis(*cell_population.GetCellsUsingLocationIndex(24u).begin());
        TS_ASSERT_EQUALS((*(cell_population.GetCellsUsingLocationIndex(24u)).begin())->IsDead(), false);


        //Test to check death of cancer cell after being in quiescent state for a prolonged period.
        (*(cell_population.GetCellsUsingLocationIndex(6u)).begin())->InitialiseCellCycleModel();
        TS_ASSERT_EQUALS((*(cell_population.GetCellsUsingLocationIndex(6u)).begin())->GetCellCycleModel()->GetCurrentCellCyclePhase(), M_PHASE);

        for(unsigned i=1; i<=4000; i++)
        {
          p_simulation_time->IncrementTimeOneStep();
          (*(cell_population.GetCellsUsingLocationIndex(6u)).begin())->GetCellCycleModel()->UpdateCellCyclePhase();
          TS_ASSERT_EQUALS((*(cell_population.GetCellsUsingLocationIndex(6u)).begin())->GetCellCycleModel()->GetCurrentCellCyclePhase(), G_ZERO_PHASE);
          TS_ASSERT_EQUALS((*(cell_population.GetCellsUsingLocationIndex(6u)).begin())->GetMutationState()->IsType<QuiescentCancerCellMutationState>(), true);
          TS_ASSERT_EQUALS((*(cell_population.GetCellsUsingLocationIndex(6u)).begin())->HasCellProperty<CellLabel>(), false);
        }

        p_simulation_time->IncrementTimeOneStep();
        (*(cell_population.GetCellsUsingLocationIndex(6u)).begin())->GetCellCycleModel()->UpdateCellCyclePhase();
        TS_ASSERT_EQUALS((*(cell_population.GetCellsUsingLocationIndex(6u)).begin())->GetCellCycleModel()->GetCurrentCellCyclePhase(), G_ZERO_PHASE);
        TS_ASSERT_EQUALS((*(cell_population.GetCellsUsingLocationIndex(6u)).begin())->HasCellProperty<CellLabel>(), true);

        cell_killer.CheckAndLabelSingleCellForApoptosis(*cell_population.GetCellsUsingLocationIndex(6u).begin());

        // Cancer cell with apoptotic cell property, hence must be dead
        TS_ASSERT_EQUALS((*(cell_population.GetCellsUsingLocationIndex(6u)).begin())->IsDead(), true);

        std::set<double> old_locations;
        // Store 'locations' of cells which are not dead
        for (std::list<CellPtr>::iterator cell_iter = r_cells.begin();
             cell_iter != r_cells.end();
             ++cell_iter)
        {
            if (!(*cell_iter)->IsDead())
            {
                Node<2>* p_node = cell_population.GetNodeCorrespondingToCell(*cell_iter);
                c_vector<double, 2> location = p_node->rGetLocation();
                old_locations.insert(location[0] + location[1]*1000);
            }
        }

        // Remove dead cells
        cell_population.RemoveDeadCells();

        // Check that dead cells are removed from the mesh
        std::set< double > new_locations;
        for (std::list<CellPtr>::iterator cell_iter = r_cells.begin();
             cell_iter != r_cells.end();
             ++cell_iter)
        {
            TS_ASSERT_EQUALS((*cell_iter)->IsDead(), false);
            Node<2>* p_node = cell_population.GetNodeCorrespondingToCell(*cell_iter);
            c_vector<double, 2> location = p_node->rGetLocation();
            new_locations.insert(location[0] + location[1]*1000);
        }

        TS_ASSERT(new_locations == old_locations);
   }
};
      ///\todo TestArchiveVasctumCellKiller()

//    void TestArchivingOfTargetedCellKiller() throw (Exception)
//    {
//        // Set up singleton classes
//        OutputFileHandler handler("archive", false);    // don't erase contents of folder
//        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "single_cell_killer.arch";
//
//        {
//            // Create an output archive
//             TargetedCellKiller<2> cell_killer(NULL, 1u);
//
//            std::ofstream ofs(archive_filename.c_str());
//            boost::archive::text_oarchive output_arch(ofs);
//
//            // Serialize via pointer
//            TargetedCellKiller<2>* const p_cell_killer = &cell_killer;
//            output_arch << p_cell_killer;
//
//            TS_ASSERT_EQUALS(p_cell_killer->GetTargetIndex(), 1u);
//        }
//
//        {
//            // Create an input archive
//            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
//            boost::archive::text_iarchive input_arch(ifs);
//
//            TargetedCellKiller<2>* p_cell_killer;
//
//            // Restore from the archive
//            input_arch >> p_cell_killer;
//
//            // Test we have restored the Target Cell correctly
//            TS_ASSERT_EQUALS(p_cell_killer->GetTargetIndex(), 1u);
//
//            delete p_cell_killer;
//       }
//    }
//
//    void TestCellKillersOutputParameters()
//    {
//        std::string output_directory = "TestCellKillersOutputParameters";
//        OutputFileHandler output_file_handler(output_directory, false);
//
//        // Test with TargetedCellKiller
//        TargetedCellKiller<2> targeted_cell_killer(NULL, 1u);
//        TS_ASSERT_EQUALS(targeted_cell_killer.GetIdentifier(), "TargetedCellKiller-2");
//
//        out_stream targeted_cell_killer_parameter_file = output_file_handler.OpenOutputFile("targeted_results.parameters");
//        targeted_cell_killer.OutputCellKillerParameters(targeted_cell_killer_parameter_file);
//        targeted_cell_killer_parameter_file->close();
//
//        std::string targeted_cell_killer_results_dir = output_file_handler.GetOutputDirectoryFullPath();
//        TS_ASSERT_EQUALS(system(("diff " + targeted_cell_killer_results_dir + "targeted_results.parameters cell_based/test/data/TestCellKillers/targeted_results.parameters").c_str()), 0);
//    }
// };

#endif /*TESTCELLKILLER_HPP_*/
