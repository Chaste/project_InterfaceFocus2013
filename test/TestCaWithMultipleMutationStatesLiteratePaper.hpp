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

#ifndef TESTCAWITHMULTIPLEMUTATIONSTATESINTERFACEFOCUS_HPP_
#define TESTCAWITHMULTIPLEMUTATIONSTATESINTERFACEFOCUS_HPP_

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

#include "WildTypeCellMutationState.hpp"
#include "CancerCellMutationState.hpp"
#include "MacrophageMutationState.hpp"

#include "MultipleCaBasedCellPopulation.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "PlaneBasedCellKiller.hpp"
#include "CheckpointArchiveTypes.hpp"
#include "OnLatticeSimulationInterfaceFocus.hpp"
#include "Warnings.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "SmartPointers.hpp"
#include "Owen2011MultipleCaUpdateRule.hpp"
#include "VasctumCellKiller.hpp"
#include "CellBasedPdeHandlerInterfaceFocus.hpp"


class TestCaWithMultipleMutationStatesInterfaceFocus : public AbstractCellBasedTestSuite
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


	void TestWithMultipleCellTypesAndSingleOccupancy() throw (Exception)
	{
	   EXIT_IF_PARALLEL;

	   // Create a simple 2D PottsMesh
	   PottsMeshGenerator<2> generator(50, 0, 0, 50, 0, 0);
	   PottsMesh<2>* p_mesh = generator.GetMesh();

	   // Create cells
	   std::vector<CellPtr> cells;
	   MAKE_PTR(StemCellProliferativeType, p_stem_type);
	   CellsGenerator<Owen2011OxygenBasedCellCycleModel, 2> cells_generator;
	   cells_generator.GenerateBasicRandom(cells, 49u, p_stem_type);

	   // Create a cell mutation state
	   boost::shared_ptr<AbstractCellProperty> p_cancer_mutation(CellPropertyRegistry::Instance()->Get<CancerCellMutationState>());
	   boost::shared_ptr<AbstractCellProperty> p_normal_mutation(CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());
	   boost::shared_ptr<AbstractCellProperty> p_macrophage_mutation(CellPropertyRegistry::Instance()->Get<MacrophageMutationState>());

	   for (unsigned i = 0; i<cells.size(); i++)
	   {
		  if (i==23 || i==24 || i==25)
			 cells[i]->SetMutationState(p_cancer_mutation);
		  else if (i==48 || i==47)
			 cells[i]->SetMutationState(p_macrophage_mutation);
		  else
			 cells[i]->SetMutationState(p_normal_mutation);
	   }

	   // Specify where cells lie
	   std::vector<unsigned> location_indices;

	   for (unsigned i=1074; i<=1374; i=i+50)
	   {
		  for(unsigned j=0; j<7; j++)
		  {
			 location_indices.push_back(i+j);
		  }
	   }

	   // Create cell population
	   MultipleCaBasedCellPopulation<2> cell_population(*p_mesh, cells, location_indices);
	   cell_population.SetOutputCellIdData(true);
	   cell_population.SetOutputCellMutationStates(true);
	   cell_population.SetOutputCellProliferativeTypes(true);
	   cell_population.SetOutputCellCyclePhases(true);
	   cell_population.SetOutputCellAncestors(true);
	   cell_population.SetOutputCellAges(true);

	   cell_population.SetDataOnAllCells("oxygen", 1.0);

	   // Set up cell-based simulation
	   OnLatticeSimulationInterfaceFocus<2> simulator(cell_population);
	   std::string output_directory = "TestWithMultipleMutationStatesAndSingleOccupancy";
	   simulator.SetOutputDirectory(output_directory);
	   simulator.SetDt(0.25);
	   simulator.SetSamplingTimestepMultiple(10);
	   simulator.SetEndTime(500.0);

	   // Set up PDE and pass to simulation via handler
	   OxygenPde<2> pde_1(cell_population, 0.1);
	   ConstBoundaryCondition<2> bc_1(100.0);
	   PdeAndBoundaryConditions<2> pde_and_bc_1(&pde_1, &bc_1, false);
	   pde_and_bc_1.SetDependentVariableName("oxygen");

	   CellBasedPdeHandlerInterfaceFocus<2> pde_handler(&cell_population);
	   pde_handler.AddPdeAndBc(&pde_and_bc_1);
	   ChastePoint<2> lower(0, 0);
	   ChastePoint<2> upper(49, 49);
	   ChasteCuboid<2> cuboid(lower, upper);
	   pde_handler.UseCoarsePdeMesh(1.0, cuboid);
	   pde_handler.SetImposeBcsOnCoarseBoundary(true);

	   simulator.SetCellBasedPdeHandler(&pde_handler);

	   // Adding update rule(s)
	   MAKE_PTR(Owen2011MultipleCaUpdateRule<2u>, p_diffusion_update_rule);
	   p_diffusion_update_rule->SetDiffusionParameter(0.02);

	   simulator.AddMultipleCaUpdateRule(p_diffusion_update_rule);

	   // Create cell killer
	   MAKE_PTR_ARGS(VasctumCellKiller<2>, cell_killer,(&cell_population));

	   simulator.AddCellKiller(cell_killer);

	   // Run simulation
	   simulator.Solve();

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

	void TestWithMultipleCellTypesAndMultipleOccupancy() throw (Exception)
	{
	   EXIT_IF_PARALLEL;

	   // Create a simple 2D PottsMesh
	   PottsMeshGenerator<2> generator(50, 0, 0, 50, 0, 0);
	   PottsMesh<2>* p_mesh = generator.GetMesh();

	   // Create cells
	   std::vector<CellPtr> cells;
	   MAKE_PTR(StemCellProliferativeType, p_stem_type);
	   CellsGenerator<Owen2011OxygenBasedCellCycleModel, 2> cells_generator;
	   cells_generator.GenerateBasicRandom(cells, 49u, p_stem_type);

	   // Create a cell mutation state
	   boost::shared_ptr<AbstractCellProperty> p_cancer_mutation(CellPropertyRegistry::Instance()->Get<CancerCellMutationState>());
	   boost::shared_ptr<AbstractCellProperty> p_normal_mutation(CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());
	   boost::shared_ptr<AbstractCellProperty> p_macrophage_mutation(CellPropertyRegistry::Instance()->Get<MacrophageMutationState>());

	   // Defining the types of cells in the lattice
	   for (unsigned i = 0; i<cells.size(); i++)
	   {
		  if (i==23 || i==24 || i==25)
			 cells[i]->SetMutationState(p_cancer_mutation);
		  else
			 if (i==48 || i==47)
			   cells[i]->SetMutationState(p_macrophage_mutation);
			 else
				cells[i]->SetMutationState(p_normal_mutation);
	   }

	   // Specify where cells lie
	   std::vector<unsigned> location_indices;

	   for(unsigned i = 1074; i<=1374; i+= 50)
	   {
		  for(unsigned j=0; j<7; j++)
		  {
			 location_indices.push_back(i+j);
		  }
		}

		// Create cell population
		MultipleCaBasedCellPopulation<2> cell_population(*p_mesh, cells, location_indices, 4);
		cell_population.SetOutputCellIdData(true);
		cell_population.SetOutputCellMutationStates(true);
		cell_population.SetOutputCellProliferativeTypes(true);
		cell_population.SetOutputCellCyclePhases(true);
		cell_population.SetOutputCellAncestors(true);
		cell_population.SetOutputCellAges(true);

		cell_population.SetDataOnAllCells("oxygen", 1.0);

		// Set up cell-based simulation
		OnLatticeSimulationInterfaceFocus<2> simulator(cell_population);
		std::string output_directory = "TestWithMultipleMutationStatesAndMultipleOccupancy";
		simulator.SetOutputDirectory(output_directory);
		simulator.SetDt(0.5);
		simulator.SetSamplingTimestepMultiple(10);
		simulator.SetEndTime(500.0);

		// Set up PDE and pass to simulation via handler
		OxygenPde<2> pde_1(cell_population, 0.1);
		ConstBoundaryCondition<2> bc_1(100.0);
		PdeAndBoundaryConditions<2> pde_and_bc_1(&pde_1, &bc_1, false);
		pde_and_bc_1.SetDependentVariableName("oxygen");

		CellBasedPdeHandlerInterfaceFocus<2> pde_handler(&cell_population);
		pde_handler.AddPdeAndBc(&pde_and_bc_1);
		//pde_handler.AddPdeAndBc(&pde_and_bc_2);
		ChastePoint<2> lower(0, 0);
		ChastePoint<2> upper(49, 49);
		ChasteCuboid<2> cuboid(lower, upper);
		pde_handler.UseCoarsePdeMesh(1.0, cuboid);
		pde_handler.SetImposeBcsOnCoarseBoundary(true);

		simulator.SetCellBasedPdeHandler(&pde_handler);

		// Create cell killer
		MAKE_PTR_ARGS(VasctumCellKiller<2>, cell_killer,(&cell_population));

		simulator.AddCellKiller(cell_killer);

		// Run simulation
		simulator.Solve();

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

   void TestCaMonolayerWithMultipleCellTypesIn3D() throw (Exception)
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
      boost::shared_ptr<AbstractCellProperty> p_macrophage_mutation(CellPropertyRegistry::Instance()->Get<MacrophageMutationState>());

      for (unsigned i = 0; i<cells.size(); i++)
      {
         cells[i]->SetMutationState(p_normal_mutation);
      }

      // Specify where cells lie
      std::vector<unsigned> location_indices;
      location_indices.push_back(3368);

      // Create cell population
      MultipleCaBasedCellPopulation<3> cell_population(*p_mesh, cells, location_indices);
      cell_population.SetOutputCellIdData(true);
 	  cell_population.SetOutputCellMutationStates(true);
      cell_population.SetOutputCellProliferativeTypes(true);
      cell_population.SetOutputCellCyclePhases(true);
      cell_population.SetOutputCellAncestors(true);
      cell_population.SetOutputCellAges(true);

      cell_population.SetDataOnAllCells("oxygen", 1.0);

      // Set up cell-based simulation
      OnLatticeSimulationInterfaceFocus<3> simulator(cell_population);
      std::string output_directory = "TestWithMultipleMutationStatesIn3D";
      simulator.SetOutputDirectory(output_directory);
      simulator.SetDt(0.5);
      simulator.SetSamplingTimestepMultiple(10);
      simulator.SetEndTime(500.0);

      // Set up PDE and pass to simulation via handler
      OxygenPde<3> pde_1(cell_population, 1.0);
      ConstBoundaryCondition<3> bc_1(30.0);
      PdeAndBoundaryConditions<3> pde_and_bc_1(&pde_1, &bc_1, false);
      pde_and_bc_1.SetDependentVariableName("oxygen");

      CellBasedPdeHandlerInterfaceFocus<3> pde_handler(&cell_population);
      pde_handler.AddPdeAndBc(&pde_and_bc_1);
      ChastePoint<3> lower(-0.5, -0.5, -0.5);
      ChastePoint<3> upper(19.5, 19.5, 19.5);
      ChasteCuboid<3> cuboid(lower, upper);
      pde_handler.UseCoarsePdeMesh(1.0, cuboid);
      pde_handler.SetImposeBcsOnCoarseBoundary(true);

      simulator.SetCellBasedPdeHandler(&pde_handler);

      // Create cell killer
      MAKE_PTR_ARGS(VasctumCellKiller<3>, cell_killer,(&cell_population));

      simulator.AddCellKiller(cell_killer);

      // Run simulation
      simulator.Solve();

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

#endif /*TESTCAWITHMULTIPLEMUTATIONSTATESINTERFACEFOCUS_HPP_*/
