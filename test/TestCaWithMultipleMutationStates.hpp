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

#ifndef TESTCAWITHMULTIPLEMUTATIONSTATES_HPP_
#define TESTCAWITHMULTIPLEMUTATIONSTATES_HPP_

#include <cxxtest/TestSuite.h>

// Must be included before other cell_based headers
#include "CellBasedSimulationArchiver.hpp"

#include "CellwiseSourcePde.hpp"
#include "ConstBoundaryCondition.hpp"
#include "PetscSetupAndFinalize.hpp"
//#include "ReplicatableVector.hpp"
#include "NumericFileComparison.hpp"
//#include "FunctionalBoundaryCondition.hpp"
#include "AveragedSourcePde.hpp"
#include "OxygenPde.hpp"


#include "PottsBasedCellPopulation.hpp"

#include "AbstractCellPopulation.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "PottsMeshGenerator.hpp"
#include "CellsGenerator.hpp"
#include "StemCellProliferativeType.hpp"
#include "TransitCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "Owen2011OxygenBasedCellCycleModel.hpp"
#include "Owen2011OxygenBasedCellCycleModelWithoutOde.hpp"
#include "VesselDerivedCellCycleModel.hpp"
#include "StochasticDurationCellCycleModel.hpp"

#include "WildTypeCellMutationState.hpp"
#include "CancerCellMutationState.hpp"
#include "MacrophageMutationState.hpp"
#include "VesselCellMutationState.hpp"
#include "StalkCellMutationState.hpp"
#include "TipCellMutationState.hpp"

#include "MultipleCaBasedCellPopulationWithVasculature.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "PlaneBasedCellKiller.hpp"
#include "OnLatticeSimulationWithVasculature.hpp"
#include "Warnings.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "SmartPointers.hpp"

#include "Owen2011MultipleCaUpdateRule.hpp"
#include "DiffusionMultipleCaUpdateRule.hpp"

#include "VasctumCellKiller.hpp"
#include "AbstractMultipleCaUpdateRule.hpp"

//  class SimplePdeForTesting : public AbstractLinearEllipticPde<2,2>
//    {
//    public:
//        double ComputeConstantInUSourceTerm(const ChastePoint<2>&, Element<2,2>* pElement)
//        {
//            return 0;
//        }
//
//        double ComputeLinearInUCoeffInSourceTerm(const ChastePoint<2>&, Element<2,2>*)
//        {
//            return 1;
//        }
//
//        c_matrix<double,2,2> ComputeDiffusionTerm(const ChastePoint<2>& )
//        {
//            return identity_matrix<double>(2u);
//        }
//    };

class TestCaWithMultipleMutationStates : public AbstractCellBasedTestSuite
{
private:

    double mLastStartTime;
    void setUp()
    {
        mLastStartTime = std::clock();
        AbstractCellBasedTestSuite::setUp();
    }
    void tearDown()
    {
        double time = std::clock();
        double elapsed_time = (time - mLastStartTime)/(CLOCKS_PER_SEC);
        std::cout << "Elapsed time: " << elapsed_time << std::endl;
        AbstractCellBasedTestSuite::tearDown();
    }

public:

    void TestCaMonolayerWithMultipleCellTypes() throw (Exception)
    {
        EXIT_IF_PARALLEL;

        // Create a simple 2D PottsMesh
        PottsMeshGenerator<2> generator(10, 0, 0, 10, 0, 0);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        MAKE_PTR(StemCellProliferativeType, p_stem_type);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);

        // Create some cell mutation states
        boost::shared_ptr<AbstractCellProperty> p_cancer_mutation(CellPropertyRegistry::Instance()->Get<CancerCellMutationState>());
        boost::shared_ptr<AbstractCellProperty> p_normal_mutation(CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());
        boost::shared_ptr<AbstractCellProperty> p_macrophage_mutation(CellPropertyRegistry::Instance()->Get<MacrophageMutationState>());
        boost::shared_ptr<AbstractCellProperty> p_vessel_mutation(CellPropertyRegistry::Instance()->Get<VesselCellMutationState>());
        boost::shared_ptr<AbstractCellProperty> p_tip_mutation(CellPropertyRegistry::Instance()->Get<TipCellMutationState>());
        boost::shared_ptr<AbstractCellProperty> p_stalk_mutation(CellPropertyRegistry::Instance()->Get<StalkCellMutationState>());

        boost::shared_ptr<AbstractCellProperty> p_state;

        // Create cells
		std::vector<CellPtr> cells;
        for (unsigned i=0; i<20u; i++)
        {
            // Way to store which stalk and tip cells are in each vessel stalk
            unsigned stalk_id = 0;

        	// Specify cell type and mutation
        //	if (i<7 || i>41)
        	{
        		p_state=p_vessel_mutation;
        	}
//            else if (i==8)
//            {
//            	p_state=p_tip_mutation;
//            	stalk_id = 1;
//            }
//            else if (i==11 || i==18)
//			{
//            	p_state=p_stalk_mutation;
//            	stalk_id = 2;
//			}
//            else if (i==25)
//            {
//            	p_state=p_tip_mutation;
//            	stalk_id=2;
//            }
//            else if (i==23 || i==24)
//			{
//            	p_state=p_cancer_mutation;
//			}
//            else
//            {
//            	p_state=p_normal_mutation;
//            }


        	AbstractCellCycleModel* p_cycle_model;
        	boost::shared_ptr<AbstractCellProliferativeType> p_cell_proliferative_type;

        	if (p_state==p_normal_mutation || p_state==p_cancer_mutation)
			{
        		//p_cycle_model = new Owen2011OxygenBasedCellCycleModel();
        		//p_cell_proliferative_type = p_stem_type;
        		p_cycle_model = new StochasticDurationCellCycleModel();
        		p_cell_proliferative_type = p_differentiated_type;

			}
        	else
        	{
        		p_cycle_model = new VesselDerivedCellCycleModel();
        		p_cell_proliferative_type = p_differentiated_type;
        	}

            p_cycle_model->SetDimension(2);
            CellPtr p_cell(new Cell(p_state, p_cycle_model));

            p_cell->SetCellProliferativeType(p_cell_proliferative_type);
        	p_cell->SetMutationState(p_state);
            //p_cell->InitialiseCellCycleModel();

        	// Use cell data to give a stalk ID to all cells
        	p_cell->GetCellData()->SetItem("StalkId", stalk_id);

        	cells.push_back(p_cell);


        }

        // Specify where cells lie
//        std::vector<unsigned> location_indices;
//        unsigned i=1074;
//        while(i<=1374)
//        {
//         for(unsigned j=0; j<7; j++)
//         {
//
//            location_indices.push_back(i+j);
//
//         }
//            i=i+50;
//        }
        std::vector<unsigned> location_indices;
        for(unsigned i=0; i<10u; i++)
        {
            location_indices.push_back(i);
            location_indices.push_back(90+i);
        }


        // Create cell population
        MultipleCaBasedCellPopulationWithVasculature<2> cell_population(*p_mesh, cells, location_indices,4);
        cell_population.SetOutputCellIdData(true);
        cell_population.SetOutputCellMutationStates(true);
        cell_population.SetOutputCellProliferativeTypes(true);
        cell_population.SetOutputCellCyclePhases(true);
        cell_population.SetOutputCellAncestors(true);
        cell_population.SetOutputCellAges(true);

        //Set CellData to all cells
	    cell_population.SetDataOnAllCells("oxygen", 100.0);
//		cell_population.SetDataOnAllCells("vegf", 0.0);

        // Set up cell-based simulation
	    OnLatticeSimulationWithVasculature<2> simulator(cell_population);
        std::string output_directory = "TestMultipleCaMonolayerWithMultipleMutationStates";
        simulator.SetOutputDirectory(output_directory);
        simulator.SetDt(0.1);
        simulator.SetSamplingTimestepMultiple(10);
        simulator.SetEndTime(500.0);

        // Set up PDE and pass to simulation via handler
//		OxygenPde<2> pde_1(cell_population, 0.1);
//		ConstBoundaryCondition<2> bc_1(100.0);
//		PdeAndBoundaryConditions<2> pde_and_bc_1(&pde_1, &bc_1, false);
//		pde_and_bc_1.SetDependentVariableName("oxygen");

//		VegfPde<2> pde_2(cell_population, 0.01, 0.01);
//		ConstBoundaryCondition<2> bc_2(0.0);
//		PdeAndBoundaryConditions<2> pde_and_bc_2(&pde_2, &bc_2, true);
//		pde_and_bc_2.SetDependentVariableName("vegf");

//		CellBasedPdeHandler<2> pde_handler(&cell_population);
//		pde_handler.AddPdeAndBc(&pde_and_bc_1);
//		pde_handler.AddPdeAndBc(&pde_and_bc_2);
//		ChastePoint<2> lower(0, 0);
//		ChastePoint<2> upper(49, 49);
//		ChasteCuboid<2> cuboid(lower, upper);
//		pde_handler.UseCoarsePdeMesh(1.0, cuboid);
//		pde_handler.SetImposeBcsOnCoarseBoundary(true);

//		simulator.SetCellBasedPdeHandler(&pde_handler);

		// Adding update rule(s)
        MAKE_PTR(Owen2011MultipleCaUpdateRule<2u>, p_diffusion_update_rule);
        p_diffusion_update_rule->SetDiffusionParameter(0.001);
        p_diffusion_update_rule->SetVEGFParameterChi(0.0001);
        p_diffusion_update_rule->SetVEGFParameterChiForMacrophages(0.0001);
        p_diffusion_update_rule->SetVEGFParameterA(0.0001);
        p_diffusion_update_rule->SetVEGFParameterB(0.0001);

        simulator.AddMultipleCaUpdateRule(p_diffusion_update_rule);

//		MAKE_PTR(DiffusionMultipleCaUpdateRule<2u>, p_diffusion_update_rule);
//		p_diffusion_update_rule->SetDiffusionParameter(0.005);
//		simulator.AddMultipleCaUpdateRule(p_diffusion_update_rule);

        // Create cell killer
//        MAKE_PTR_ARGS(VasctumCellKiller<2>, cell_killer,(&cell_population));
//
//        simulator.AddCellKiller(cell_killer);

        // Run simulation
        simulator.Solve();

        // Check the number of cells
//        TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetNumRealCells(), 17u);///\todo #2066 Check this!

        // Test no deaths and some births
//        TS_ASSERT_EQUALS(simulator.GetNumBirths(), 15u);///\todo #2066 Check this!
//        TS_ASSERT_EQUALS(simulator.GetNumDeaths(), 49u);

#ifdef CHASTE_VTK
//Test that VTK writer has produced a file
OutputFileHandler output_file_handler(output_directory, false);
std::string results_dir = output_file_handler.GetOutputDirectoryFullPath();

// Initial condition file
FileFinder vtk_file(results_dir + "results_from_time_0/results_0.vtu", RelativeTo::Absolute);
TS_ASSERT(vtk_file.Exists());

// Final file
FileFinder vtk_file2(results_dir + "results_from_time_0/results_200.vtu", RelativeTo::Absolute);
TS_ASSERT(vtk_file2.Exists());
#endif //CHASTE_VTK

    }

    void noTestCaMonolayerWithMultipleCellTypesIn3D() throw (Exception)
	{
		EXIT_IF_PARALLEL;

		// Create a simple 3D PottsMesh
		PottsMeshGenerator<3> generator(20, 0, 0, 20, 0, 0, 20, 0, 0);
		PottsMesh<3>* p_mesh = generator.GetMesh();

		// Create cells
		std::vector<CellPtr> cells;
		MAKE_PTR(StemCellProliferativeType, p_stem_type);
		CellsGenerator<Owen2011OxygenBasedCellCycleModel, 3> cells_generator;
		cells_generator.GenerateBasicRandom(cells, 1u, p_stem_type);

		// Create a cell mutation state
		boost::shared_ptr<AbstractCellProperty> p_cancer_mutation(CellPropertyRegistry::Instance()->Get<CancerCellMutationState>());
		boost::shared_ptr<AbstractCellProperty> p_normal_mutation(CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());
		boost::shared_ptr<AbstractCellProperty> p_macrophage(CellPropertyRegistry::Instance()->Get<MacrophageMutationState>());

		for (unsigned i = 0; i<cells.size(); i++)
		{
//            if (i==23 || i==24)
//              cells[i]->SetMutationState(p_cancer_mutation);
//            else if (i==48 || i==47)
//                cells[i]->SetMutationState(p_macrophage);
//            else
			  cells[i]->SetMutationState(p_normal_mutation);
		}

		// Specify where cells lie
		std::vector<unsigned> location_indices;
//		unsigned k=400;
//		while(k<=500)
//		{
//			unsigned i=50;
//			while(i<=60)
//			{
//				for(unsigned j=5; j<7; j++)
//				{
//					location_indices.push_back(k+i+j);
//				}
//					i=i+10;
//			}
//				k=k+100;
//		}
		location_indices.push_back(3368);

		// Create cell population
		MultipleCaBasedCellPopulation<3> cell_population(*p_mesh, cells, location_indices);
		cell_population.SetOutputCellIdData(true);
		cell_population.SetOutputCellMutationStates(true);
		cell_population.SetOutputCellProliferativeTypes(true);
		cell_population.SetOutputCellCyclePhases(true);
		cell_population.SetOutputCellAncestors(true);
		cell_population.SetOutputCellAges(true);

		//Set CellData to all cells
		cell_population.SetDataOnAllCells("oxygen",30.0);
//		cell_population.SetDataOnAllCells("vegf", 0.0);

		// Set up cell-based simulation
		OnLatticeSimulation<3> simulator(cell_population);
		std::string output_directory = "TestMultipleCaMonolayerWithMultipleMutationStatesIn3D";
		simulator.SetOutputDirectory(output_directory);
		simulator.SetDt(0.5);
		simulator.SetSamplingTimestepMultiple(10);
		simulator.SetEndTime(10.0);

		// Set up PDE and pass to simulation via handler
		OxygenPde<3> pde_1(cell_population, 1.0);
		ConstBoundaryCondition<3> bc_1(30.0);
		PdeAndBoundaryConditions<3> pde_and_bc_1(&pde_1, &bc_1, false);
		pde_and_bc_1.SetDependentVariableName("oxygen");

//		VegfPde<2> pde_2(cell_population, 0.01, 0.01);
//		ConstBoundaryCondition<2> bc_2(0.0);
//		PdeAndBoundaryConditions<2> pde_and_bc_2(&pde_2, &bc_2, true);
//		pde_and_bc_2.SetDependentVariableName("vegf");

		CellBasedPdeHandler<3> pde_handler(&cell_population);
		pde_handler.AddPdeAndBc(&pde_and_bc_1);
//		pde_handler.AddPdeAndBc(&pde_and_bc_2);
		ChastePoint<3> lower(-0.5, -0.5, -0.5);
		ChastePoint<3> upper(19.5, 19.5, 19.5);
		ChasteCuboid<3> cuboid(lower, upper);
		pde_handler.UseCoarsePdeMesh(1.0, cuboid);
		pde_handler.SetImposeBcsOnCoarseBoundary(true);

		simulator.SetCellBasedPdeHandler(&pde_handler);

//		// Adding update rule(s)
//        MAKE_PTR(Owen2011MultipleCaUpdateRule<2u>, p_diffusion_update_rule);
//        p_diffusion_update_rule->SetDiffusionParameter(0.05);
//        p_diffusion_update_rule->SetVEGFParameterChi(0.05);
//        p_diffusion_update_rule->SetVEGFParameterChiForMacrophages(0.6);
//        p_diffusion_update_rule->SetVEGFParameterA(0.5);
//        p_diffusion_update_rule->SetVEGFParameterB(0.03);
//
//        simulator.AddMultipleCaUpdateRule(p_diffusion_update_rule);

		// Create cell killer
		MAKE_PTR_ARGS(VasctumCellKiller<3>, cell_killer,(&cell_population));

		simulator.AddCellKiller(cell_killer);

		// Run simulation
		simulator.Solve();

		// Check the number of cells
//        TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetNumRealCells(), 17u);///\todo #2066 Check this!

		// Test no deaths and some births
//        TS_ASSERT_EQUALS(simulator.GetNumBirths(), 15u);///\todo #2066 Check this!
//        TS_ASSERT_EQUALS(simulator.GetNumDeaths(), 49u);

	#ifdef CHASTE_VTK
//Test that VTK writer has produced a file
OutputFileHandler output_file_handler(output_directory, false);
std::string results_dir = output_file_handler.GetOutputDirectoryFullPath();

// Initial condition file
FileFinder vtk_file(results_dir + "results_from_time_0/results_0.vtu", RelativeTo::Absolute);
TS_ASSERT(vtk_file.Exists());

// Final file
FileFinder vtk_file2(results_dir + "results_from_time_0/results_200.vtu", RelativeTo::Absolute);
TS_ASSERT(vtk_file2.Exists());
 #endif //CHASTE_VTK

	}
};

#endif /*TESTCAWITHMULTIPLEMUTATIONSTATES_HPP_*/
