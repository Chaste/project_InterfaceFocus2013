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

#ifndef TESTSIMULATIONWITHPDES_HPP_
#define TESTSIMULATIONWITHPDES_HPP_

#define PI 3.141592653589793238462643383279502884197169399375105820974944592308

#include <cxxtest/TestSuite.h>

// Must be included before other cell_based headers
#include "CellBasedSimulationArchiver.hpp"

#include "AbstractCellPopulation.hpp"
#include "OnLatticeSimulation.hpp"
#include "CellsGenerator.hpp"
#include "FixedDurationGenerationBasedCellCycleModel.hpp"
#include "ConstBoundaryCondition.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "ReplicatableVector.hpp"
#include "NumericFileComparison.hpp"
#include "WildTypeCellMutationState.hpp"
#include "OxygenPde.hpp"
#include "AveragedSourcePde.hpp"
#include "SmartPointers.hpp"
#include "PottsMeshGenerator.hpp"
#include "MultipleCaBasedCellPopulation.hpp"
#include "DiffusionMultipleCaUpdateRule.hpp"
#include "Warnings.hpp"

#include "Debug.hpp"

class SimplePdeForTesting : public AbstractLinearEllipticPde<2,2>
{
public:
    double ComputeConstantInUSourceTerm(const ChastePoint<2>&, Element<2,2>* pElement)
    {
        return -1.0;
    }

    double ComputeLinearInUCoeffInSourceTerm(const ChastePoint<2>&, Element<2,2>*)
    {
        return 0.0;
    }

    c_matrix<double,2,2> ComputeDiffusionTerm(const ChastePoint<2>& )
    {
        return identity_matrix<double>(2);
    }
};

class TestSimulationWithPdes : public AbstractCellBasedTestSuite
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

    void TestMultipleCaBasedWithCoarseMesh() throw(Exception)
    {

    	EXIT_IF_PARALLEL;

        // Create a simple 2D PottsMesh
        PottsMeshGenerator<2> generator(5, 0, 0, 5, 0, 0);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, 1);

        std::vector<unsigned> location_indices;
        location_indices.push_back(12);
//        location_indices.push_back(11);

        // Create cell population
        MultipleCaBasedCellPopulation<2> cell_population(*p_mesh, cells, location_indices);

        // Initialize CellData
        /*
         * Since values are first passed in to CellData before it is updated in UpdateAtEndOfTimeStep(),
         * we need to pass it some initial conditions.
         */
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            cell_iter->GetCellData()->SetItem("nutrient", 1);
        }

        // Set up cell-based simulation
        OnLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestMultipleCaBasedCellPopulationWithPdes");
//        simulator.SetSamplingTimestepMultiple(10);
        simulator.SetEndTime(0.1);


        // Set up PDE and pass to simulation via handler (zero uptake to check analytic solution)
        OxygenPde<2> pde(cell_population, 0.0);
        ConstBoundaryCondition<2> bc(0.0);
        PdeAndBoundaryConditions<2> pde_and_bc(&pde, &bc, true);
        pde_and_bc.SetDependentVariableName("nutrient");

        CellBasedPdeHandler<2> pde_handler(&cell_population);
        pde_handler.AddPdeAndBc(&pde_and_bc);

        // Create coarse mesh and centre it on the centre of the Potts Mesh
        c_vector<double,2> centre_of_potts_mesh = zero_vector<double>(2);
        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            centre_of_potts_mesh += p_mesh->GetNode(i)->rGetLocation();
        }
        centre_of_potts_mesh /= p_mesh->GetNumNodes();
        ChastePoint<2> lower(0.0, 0.0);
        ChastePoint<2> upper(4.0, 4.0);

        c_vector<double, 2> translation = 0.5*(lower.rGetLocation() + upper.rGetLocation()) - centre_of_potts_mesh;
        lower.rGetLocation() -= translation;
        upper.rGetLocation() -= translation;
        ChasteCuboid<2> cuboid(lower, upper);
        pde_handler.UseCoarsePdeMesh(1.0, cuboid);
        pde_handler.SetImposeBcsOnCoarseBoundary(true);
        simulator.SetCellBasedPdeHandler(&pde_handler);

        MAKE_PTR(DiffusionMultipleCaUpdateRule<2>, p_diffusion_update_rule);
        p_diffusion_update_rule->SetDiffusionParameter(0.5);
        simulator.AddMultipleCaUpdateRule(p_diffusion_update_rule);

        // Solve the system
        simulator.Solve();

//        // Test solution is constant
//        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
//             cell_iter != cell_population.End();
//             ++cell_iter)
//        {
//            double analytic_solution = 1.0;
//            // Test that PDE solver is working correctly
//            TS_ASSERT_DELTA(cell_iter->GetCellData()->GetItem("nutrient"), analytic_solution, 1e-2);
//        }

        // Find centre of coarse PDE mesh
        c_vector<double,2> centre_of_coarse_pde_mesh = zero_vector<double>(2);
        TetrahedralMesh<2,2>* p_coarse_mesh = simulator.GetCellBasedPdeHandler()->GetCoarsePdeMesh();
        for (unsigned i=0; i<p_coarse_mesh->GetNumNodes(); i++)
        {
            centre_of_coarse_pde_mesh += p_coarse_mesh->GetNode(i)->rGetLocation();
        }
        centre_of_coarse_pde_mesh /= p_coarse_mesh->GetNumNodes();
        c_vector<double,2> centre_diff = centre_of_coarse_pde_mesh - centre_of_potts_mesh;

        // Test that the two centres match
        TS_ASSERT_DELTA(norm_2(centre_diff), 0.0, 1e-4);

        // Test FindCoarseElementContainingCell() and initialisation of mCellPdeElementMap
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
            cell_iter != cell_population.End();
            ++cell_iter)
        {
            unsigned containing_element_index = simulator.GetCellBasedPdeHandler()->mCellPdeElementMap[*cell_iter];
            TS_ASSERT_LESS_THAN(containing_element_index, p_coarse_mesh->GetNumElements());
            TS_ASSERT_EQUALS(containing_element_index, simulator.GetCellBasedPdeHandler()->FindCoarseElementContainingCell(*cell_iter));
        }

//        // Test solution is constant
//        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
//             cell_iter != cell_population.End();
//             ++cell_iter)
//        {
//            double analytic_solution = 1.0;
//
//            // Test that PDE solver is working correctly on both pdes
//            TS_ASSERT_DELTA(cell_iter->GetCellData()->GetItem("nutrient"), analytic_solution, 1e-2);
//        }
    }

    void NoTestMultipleCaBasedWithCoarseMeshTwoEquations() throw(Exception)
    {
        EXIT_IF_PARALLEL;

        // Create a simple 2D PottsMesh
        PottsMeshGenerator<2> generator(5, 0, 0, 5, 0, 0);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells,1);

		std::vector<unsigned> location_indices;
		location_indices.push_back(12);
		//        location_indices.push_back(11);

		// Create cell population
		MultipleCaBasedCellPopulation<2> cell_population(*p_mesh, cells, location_indices);

        cell_population.SetDataOnAllCells("quantity 1", 1.0);
        cell_population.SetDataOnAllCells("quantity 2", 1.0);

        // Set up cell-based simulation
        OnLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestMultipleCaBasedCellPopulationWithTwoPdes");
        simulator.SetEndTime(0.1);

        // Set up PDE and pass to simulation via handler (zero uptake to check analytic solution)
        OxygenPde<2> pde_1(cell_population, 0.0);
        ConstBoundaryCondition<2> bc_1(0.0);
        PdeAndBoundaryConditions<2> pde_and_bc_1(&pde_1, &bc_1, true);
        pde_and_bc_1.SetDependentVariableName("quantity 1");

        VegfPde<2> pde_2(cell_population, 0.0, 0.0);
        ConstBoundaryCondition<2> bc_2(0.0);
        PdeAndBoundaryConditions<2> pde_and_bc_2(&pde_2, &bc_2, true);
        pde_and_bc_2.SetDependentVariableName("quantity 2");

        CellBasedPdeHandler<2> pde_handler(&cell_population);
        pde_handler.AddPdeAndBc(&pde_and_bc_1);
        pde_handler.AddPdeAndBc(&pde_and_bc_2);
        ChastePoint<2> lower(0.0, 0.0);
        ChastePoint<2> upper(4.0, 4.0);
        ChasteCuboid<2> cuboid(lower, upper);
        pde_handler.UseCoarsePdeMesh(0.1, cuboid);
        pde_handler.SetImposeBcsOnCoarseBoundary(true);

        simulator.SetCellBasedPdeHandler(&pde_handler);

        MAKE_PTR(DiffusionMultipleCaUpdateRule<2>, p_diffusion_update_rule);
        p_diffusion_update_rule->SetDiffusionParameter(0.5);
        simulator.AddMultipleCaUpdateRule(p_diffusion_update_rule);

        // Solve the system
        simulator.Solve();

        // Test solution is constant
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            double analytic_solution = 1.0;

            // Test that PDE solver is working correctly on both pdes
            TS_ASSERT_DELTA(cell_iter->GetCellData()->GetItem("quantity 1"), analytic_solution, 1e-2);
            TS_ASSERT_DELTA(cell_iter->GetCellData()->GetItem("quantity 2"), analytic_solution, 1e-2);
        }
    }

    void NoTestPdeOutput() throw(Exception)
    {
    	EXIT_IF_PARALLEL;

        OutputFileHandler handler("TestPottsBasedCellPopulationWithPdes", false);
        std::string results_dir = handler.GetOutputDirectoryFullPath() + "results_from_time_0";

        NumericFileComparison comp_nut(results_dir + "/results.vizcoarsepdesolution", "cell_based/test/data/TestPottsBasedCellPopulationWithPdes/results.vizcoarsepdesolution");
        TS_ASSERT(comp_nut.CompareFiles());
        TS_ASSERT_EQUALS(system(("diff " + results_dir + "/results.vizcoarsepdesolution cell_based/test/data/TestPottsBasedCellPopulationWithPdes/results.vizcoarsepdesolution").c_str()), 0);
    }


};

#endif /*TESTSIMULATIONWITHPDES_HPP_*/