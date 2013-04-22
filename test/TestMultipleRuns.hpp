#ifndef TESTMULTIPLERUNS_HPP_
#define TESTMULTIPLERUNS_HPP_

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


#include "PottsBasedCellPopulation.hpp"

#include "AbstractCellPopulation.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "PottsMeshGenerator.hpp"
#include "CellsGenerator.hpp"
#include "FixedDurationGenerationBasedCellCycleModel.hpp"
#include "Owen2011OxygenBasedCellCycleModel.hpp"
//#include "StochasticDurationGenerationBasedCellCycleModel.hpp"
#include "WildTypeCellMutationState.hpp"
#include "MultipleCaBasedCellPopulation.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "PlaneBasedCellKiller.hpp"
#include "OnLatticeSimulation.hpp"
#include "Warnings.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "SmartPointers.hpp"
#include "Owen2011MultipleCaUpdateRule.hpp"
#include <stdlib.h>
#include "Debug.hpp"


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

    /**
     * For use in TestOnLatticeSimulationWithPdes::TestWithBoundaryConditionVaryingInTime.
     */
    double bc_func(const ChastePoint<2>& p)
    {
        double value = SimulationTime::Instance()->GetTime();
        return value;
    }

class TestMultipleRuns : public AbstractCellBasedTestSuite
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

    void RandomlyLabelCells(std::vector<CellPtr>& rCells, boost::shared_ptr<AbstractCellProperty> pLabel, double labelledRatio)
    {
        for (unsigned i = 0; i<rCells.size(); i++)
        {
            if (RandomNumberGenerator::Instance()->ranf() < labelledRatio)
            {
                rCells[i]->AddCellProperty(pLabel);
            }
        }
    }

public:

    void CaMonolayerWithBirthForMultipleRuns(char filename[100], int i) throw (Exception)
    {
    	//SimulationTime::Instance();
    	if (i > 0)
    	   SimulationTime::Instance()->SetStartTime(0.0);

     	// Create a simple 2D PottsMesh
		PottsMeshGenerator<2> generator(100, 0, 0, 100, 0, 0);
		PottsMesh<2>* p_mesh = generator.GetMesh();

		// Create cells
		std::vector<CellPtr> cells;
        MAKE_PTR(StemCellProliferativeType, p_stem_type);
		CellsGenerator<Owen2011OxygenBasedCellCycleModel, 2> cells_generator;
		cells_generator.GenerateBasicRandom(cells, 1u, p_stem_type);

		 // Set up oxygen_concentration
		 double oxygen_concentration = 1.0;

		// Specify where cells lie
		std::vector<unsigned> location_indices;
		location_indices.push_back(5050u);

		// Create cell population
		MultipleCaBasedCellPopulation<2> cell_population(*p_mesh, cells, location_indices);
        cell_population.SetOutputCellIdData(true);
        cell_population.SetOutputCellMutationStates(true);
        cell_population.SetOutputCellProliferativeTypes(true);
        cell_population.SetOutputCellCyclePhases(true);
        cell_population.SetOutputCellAncestors(true);
        cell_population.SetOutputCellAges(true);

        //Set CellData to all cells
        cell_population.SetDataOnAllCells("oxygen", oxygen_concentration);

        // Set up cell-based simulation
		OnLatticeSimulation<2> simulator(cell_population);
        std::string output_directory = filename;
		simulator.SetOutputDirectory(output_directory);
        simulator.SetDt(0.5);
        simulator.SetEndTime(480);

        simulator.SetSamplingTimestepMultiple(80);

        // Adding update rule(s).
        MAKE_PTR(Owen2011MultipleCaUpdateRule<2u>, p_diffusion_update_rule);
        p_diffusion_update_rule->SetDiffusionParameter((0.3e-9)*60); // converting time units
        p_diffusion_update_rule->SetVEGFParameterChi(0); // No chemotaxis now
        p_diffusion_update_rule->SetVEGFParameterA(0);
        p_diffusion_update_rule->SetVEGFParameterB(0);

        simulator.AddMultipleCaUpdateRule(p_diffusion_update_rule);

        // Run simulation
        simulator.Solve();

        /*
		 * SimulationTime::Destroy() must be called at the end of the test. If not, when
		 * SimulationTime::Instance()->SetStartTime(0.0); is called at the beginning of the next test in
		 * this file, an assertion will be triggered.
		 */
		SimulationTime::Destroy();
  }

    //////////////////////////////////////////////////////////////////////////////////////////

  	void TestCaMonolayerWithBirthMultipleRuns() throw (Exception)
	{
		for(unsigned int i = 0; i <10; i++)
		{
			/* Creating the names for the multiple result folders
			 * the first folder is TestMultipleCaMonolayerWithBirthMultipleRuns_0
			 * second TestMultipleCaMonolayerWithBirthMultipleRuns_1,... etc
			 * There will be a matlab script that will read the results.vizcellphases for each folder
			 * and print the mean of all runs in order to validate chaste results with the vasctum implementation
			 */
			char buffer[10];
			sprintf (buffer, "%d", i);
			char c[100] = "TestMultipleCaMonolayerWithBirthMultipleRuns_";
			strcat(c,buffer);

			// Calling the simulation
			CaMonolayerWithBirthForMultipleRuns(c,i);
		}
   }

  	//////////////////////////////////////////////////////////////////////////////////////////////

    void MultipleCaSingleCellRandomMovement(char filename[100], int i) throw (Exception)
    {
    	//SimulationTime::Instance();
    	if (i > 0)
    	   SimulationTime::Instance()->SetStartTime(0.0);

    	// Create a simple 2D PottsMesh
        PottsMeshGenerator<2> generator(50, 0, 0, 50, 0, 0);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        CellsGenerator<Owen2011OxygenBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, 1u, p_diff_type);

        // Set up oxygen_concentration
        double oxygen_concentration = 0.0;

        // Specify where cells lie
        std::vector<unsigned> location_indices;
        location_indices.push_back(1250u);

        // Create cell population
        MultipleCaBasedCellPopulation<2> cell_population(*p_mesh, cells, location_indices);

        //Set CellData to all cells
        cell_population.SetDataOnAllCells("oxygen", oxygen_concentration);

        // Set up cell-based simulation
 		OnLatticeSimulation<2> simulator(cell_population);
        std::string output_directory = filename;
 		simulator.SetOutputDirectory(output_directory);
        simulator.SetDt(0.35);
        simulator.SetEndTime(400);

        /*
         * Adding update rule(s).
         */
        MAKE_PTR(Owen2011MultipleCaUpdateRule<2u>, p_diffusion_update_rule);
        p_diffusion_update_rule->SetDiffusionParameter(0.03);
        simulator.AddMultipleCaUpdateRule(p_diffusion_update_rule);

        // Run simulation
        simulator.Solve();

        // Check that the same number of cells
        TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetNumRealCells(), 1u);

        // Test no births or deaths
        TS_ASSERT_EQUALS(simulator.GetNumBirths(), 0u);
        TS_ASSERT_EQUALS(simulator.GetNumDeaths(), 0u);

		/*
		 * SimulationTime::Destroy() must be called at the end of the test. If not, when
		 * SimulationTime::Instance()->SetStartTime(0.0); is called at the beginning of the next test in
		 * this file, an assertion will be triggered.
		 */
		SimulationTime::Destroy();
    }


 	void TestMultipleCaSingleCellRandomMovementMultipleTimes() throw (Exception)
	{
		for(unsigned int i = 0; i < 10; i++)
		{
			char buffer[10];
			sprintf (buffer, "%d", i);
			char c[100] = "TestMultipleCaSingleCellRandomMovementMultipleTimes_";
			strcat(c,buffer);

			// Calling the simulation
			MultipleCaSingleCellRandomMovement(c,i);
		}
   }

};

#endif /* TESTMULTIPLERUNS_HPP_ */
