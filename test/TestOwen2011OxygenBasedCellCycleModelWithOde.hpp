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

#ifndef TESTODEBASEDCELLCYCLEMODELS_HPP_
#define TESTODEBASEDCELLCYCLEMODELS_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <fstream>
#include <boost/shared_ptr.hpp>

#include "Owen2011OxygenBasedCellCycleModel.hpp"

#include "AbstractCellMutationState.hpp"
#include "WildTypeCellMutationState.hpp"
#include "CancerCellMutationState.hpp"
#include "StemCellProliferativeType.hpp"
#include "ApcTwoHitCellMutationState.hpp"
#include "BetaCateninOneHitCellMutationState.hpp"

#include "OutputFileHandler.hpp"
#include "CheckReadyToDivideAndPhaseIsUpdated.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "SmartPointers.hpp"

/**
 * This class contains tests for methods on classes
 * inheriting from AbstractOdeBasedCellCycleModel.
 */
class TestOwen2011OxygenBasedCellCycleModelWithOde : public AbstractCellBasedTestSuite
{
public:

    void TestOwen2011OxygenBasedCellCycleModelWithOdeForNormalCells() throw(Exception)
    {
        // Set up SimulationTime
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(75.0, 150);

        // Set up oxygen_concentration
        double oxygen_concentration = 1.0;

        // Create cell-cycle models
        Owen2011OxygenBasedCellCycleModel* p_model_1d = new Owen2011OxygenBasedCellCycleModel();
        p_model_1d->SetDimension(1);

        Owen2011OxygenBasedCellCycleModel* p_model_2d = new Owen2011OxygenBasedCellCycleModel();
        p_model_2d->SetDimension(2);

        Owen2011OxygenBasedCellCycleModel* p_model_3d = new Owen2011OxygenBasedCellCycleModel();
        p_model_3d->SetDimension(3);

        // Create cells
        MAKE_PTR(WildTypeCellMutationState, p_state);
		MAKE_PTR(StemCellProliferativeType, p_stem_type);

        CellPtr p_cell_1d(new Cell(p_state, p_model_1d));
        p_cell_1d->GetCellData()->SetItem("oxygen", oxygen_concentration);
        p_cell_1d->SetCellProliferativeType(p_stem_type);
        p_cell_1d->InitialiseCellCycleModel();

        CellPtr p_cell_2d(new Cell(p_state, p_model_2d));
        p_cell_2d->GetCellData()->SetItem("oxygen", oxygen_concentration);
        p_cell_2d->SetCellProliferativeType(p_stem_type);
        p_cell_2d->InitialiseCellCycleModel();

        CellPtr p_cell_3d(new Cell(p_state, p_model_3d));
        p_cell_3d->GetCellData()->SetItem("oxygen", oxygen_concentration);
        p_cell_3d->SetCellProliferativeType(p_stem_type);
        p_cell_3d->InitialiseCellCycleModel();

        // For coverage, we create another cell-cycle model that is identical to p_model_2d except for the ODE solver
        boost::shared_ptr<CellCycleModelOdeSolver<Owen2011OxygenBasedCellCycleModel, RungeKutta4IvpOdeSolver> >
        p_solver(CellCycleModelOdeSolver<Owen2011OxygenBasedCellCycleModel, RungeKutta4IvpOdeSolver>::Instance());
        p_solver->Initialise();

        Owen2011OxygenBasedCellCycleModel* p_other_model_2d = new Owen2011OxygenBasedCellCycleModel(p_solver);
        p_other_model_2d->SetDimension(2);

        CellPtr p_other_cell_2d(new Cell(p_state, p_other_model_2d));
        p_other_cell_2d->GetCellData()->SetItem("oxygen", oxygen_concentration);
        p_other_cell_2d->SetCellProliferativeType(p_stem_type);
        p_other_cell_2d->InitialiseCellCycleModel();

        // Check oxygen concentration is correct in cell-cycle model
        TS_ASSERT_DELTA(p_model_2d->GetProteinConcentrations()[1], 1.0, 1e-5);
        TS_ASSERT_EQUALS(p_model_2d->ReadyToDivide(), false);

        TS_ASSERT_DELTA(p_other_model_2d->GetProteinConcentrations()[1], 1.0, 1e-5);
        TS_ASSERT_EQUALS(p_other_model_2d->ReadyToDivide(), false);

        // Divide the cells
        Owen2011OxygenBasedCellCycleModel* p_model_1d_2 = static_cast<Owen2011OxygenBasedCellCycleModel*> (p_model_1d->CreateCellCycleModel());
        CellPtr p_cell_1d_2(new Cell(p_state, p_model_1d_2));
        p_cell_1d_2->GetCellData()->SetItem("oxygen", oxygen_concentration);
        p_cell_1d_2->SetCellProliferativeType(p_stem_type);

        Owen2011OxygenBasedCellCycleModel* p_model_2d_2 = static_cast<Owen2011OxygenBasedCellCycleModel*> (p_model_2d->CreateCellCycleModel());
        CellPtr p_cell_2d_2(new Cell(p_state, p_model_2d_2));
        p_cell_2d_2->GetCellData()->SetItem("oxygen", oxygen_concentration);;
        p_cell_2d_2->SetCellProliferativeType(p_stem_type);

        Owen2011OxygenBasedCellCycleModel* p_model_3d_2 = static_cast<Owen2011OxygenBasedCellCycleModel*> (p_model_3d->CreateCellCycleModel());
        CellPtr p_cell_3d_2(new Cell(p_state, p_model_3d_2));
        p_cell_3d_2->GetCellData()->SetItem("oxygen", oxygen_concentration);
        p_cell_3d_2->SetCellProliferativeType(p_stem_type);

        for (unsigned i=0; i<100; i++)
        {
            p_simulation_time->IncrementTimeOneStep();
            TS_ASSERT_EQUALS(p_model_1d->ReadyToDivide(), false);
            TS_ASSERT_EQUALS(p_model_2d->ReadyToDivide(), false);
            TS_ASSERT_EQUALS(p_model_3d->ReadyToDivide(), false);
        }

        p_simulation_time->IncrementTimeOneStep();

        TS_ASSERT_EQUALS(p_model_1d->ReadyToDivide(), true);
        TS_ASSERT_EQUALS(p_model_2d->ReadyToDivide(), true);
        TS_ASSERT_EQUALS(p_model_3d->ReadyToDivide(), true);


        TS_ASSERT_THROWS_NOTHING(p_model_2d->ResetForDivision());

        // For coverage, create a 1D model
        Owen2011OxygenBasedCellCycleModel* p_cell_model3 = new Owen2011OxygenBasedCellCycleModel();
        p_cell_model3->SetDimension(1);

        CellPtr p_cell3(new Cell(p_state, p_cell_model3));
        p_cell3->GetCellData()->SetItem("oxygen", oxygen_concentration);
        p_cell3->SetCellProliferativeType(p_stem_type);
        p_cell3->InitialiseCellCycleModel();

        TS_ASSERT_DELTA(p_cell_model3->GetProteinConcentrations()[1], 1.0, 1e-5);
        TS_ASSERT_EQUALS(p_cell_model3->ReadyToDivide(), false);
        p_simulation_time->IncrementTimeOneStep();
        TS_ASSERT_EQUALS(p_cell_model3->ReadyToDivide(), false);
    }

    void TestOwen2011OxygenBasedCellCycleModelWithOdeForCancerCells() throw(Exception)
    {
           // Set up SimulationTime
           SimulationTime* p_simulation_time = SimulationTime::Instance();
           p_simulation_time->SetEndTimeAndNumberOfTimeSteps(75.0, 150);

           // Set up oxygen_concentration
           double oxygen_concentration = 10.0;

           // Create cell-cycle models
           Owen2011OxygenBasedCellCycleModel* p_model_1d = new Owen2011OxygenBasedCellCycleModel();
           p_model_1d->SetDimension(1);

           Owen2011OxygenBasedCellCycleModel* p_model_2d = new Owen2011OxygenBasedCellCycleModel();
           p_model_2d->SetDimension(2);

           Owen2011OxygenBasedCellCycleModel* p_model_3d = new Owen2011OxygenBasedCellCycleModel();
           p_model_3d->SetDimension(3);

           // Create cells
           MAKE_PTR(CancerCellMutationState, p_state);
       	   MAKE_PTR(StemCellProliferativeType, p_stem_type);

           CellPtr p_cell_1d(new Cell(p_state, p_model_1d));
           p_cell_1d->GetCellData()->SetItem("oxygen", oxygen_concentration);
           p_cell_1d->SetCellProliferativeType(p_stem_type);
           p_cell_1d->InitialiseCellCycleModel();
           TS_ASSERT_EQUALS(p_state->IsType<CancerCellMutationState>(), true);

           CellPtr p_cell_2d(new Cell(p_state, p_model_2d));
           p_cell_2d->GetCellData()->SetItem("oxygen", oxygen_concentration);
           p_cell_2d->SetCellProliferativeType(p_stem_type);
           p_cell_2d->InitialiseCellCycleModel();
           TS_ASSERT_EQUALS(p_state->IsType<CancerCellMutationState>(), true);

           CellPtr p_cell_3d(new Cell(p_state, p_model_3d));
           p_cell_3d->GetCellData()->SetItem("oxygen", oxygen_concentration);
           p_cell_3d->SetCellProliferativeType(p_stem_type);
           p_cell_3d->InitialiseCellCycleModel();
           TS_ASSERT_EQUALS(p_state->IsType<CancerCellMutationState>(), true);

           // For coverage, we create another cell-cycle model that is identical to p_model_2d except for the ODE solver
           boost::shared_ptr<CellCycleModelOdeSolver<Owen2011OxygenBasedCellCycleModel, RungeKutta4IvpOdeSolver> >
           p_solver(CellCycleModelOdeSolver<Owen2011OxygenBasedCellCycleModel, RungeKutta4IvpOdeSolver>::Instance());
           p_solver->Initialise();

           Owen2011OxygenBasedCellCycleModel* p_other_model_2d = new Owen2011OxygenBasedCellCycleModel(p_solver);
           p_other_model_2d->SetDimension(2);

           CellPtr p_other_cell_2d(new Cell(p_state, p_other_model_2d));
           p_other_cell_2d->GetCellData()->SetItem("oxygen", oxygen_concentration);
           p_other_cell_2d->SetCellProliferativeType(p_stem_type);
           p_other_cell_2d->InitialiseCellCycleModel();

           // Check oxygen concentration is correct in cell-cycle model
           TS_ASSERT_DELTA(p_model_2d->GetProteinConcentrations()[1], 10.0, 1e-5);
           TS_ASSERT_EQUALS(p_model_2d->ReadyToDivide(), false);

           TS_ASSERT_DELTA(p_other_model_2d->GetProteinConcentrations()[1], 10.0, 1e-5);
           TS_ASSERT_EQUALS(p_other_model_2d->ReadyToDivide(), false);

           // Divide the cells
           Owen2011OxygenBasedCellCycleModel* p_model_1d_2 = static_cast<Owen2011OxygenBasedCellCycleModel*> (p_model_1d->CreateCellCycleModel());
           CellPtr p_cell_1d_2(new Cell(p_state, p_model_1d_2));
           p_cell_1d_2->GetCellData()->SetItem("oxygen", oxygen_concentration);
           p_cell_1d_2->SetCellProliferativeType(p_stem_type);

           Owen2011OxygenBasedCellCycleModel* p_model_2d_2 = static_cast<Owen2011OxygenBasedCellCycleModel*> (p_model_2d->CreateCellCycleModel());
           CellPtr p_cell_2d_2(new Cell(p_state, p_model_2d_2));
           p_cell_2d_2->GetCellData()->SetItem("oxygen", oxygen_concentration);
           p_cell_2d_2->SetCellProliferativeType(p_stem_type);

           Owen2011OxygenBasedCellCycleModel* p_model_3d_2 = static_cast<Owen2011OxygenBasedCellCycleModel*> (p_model_3d->CreateCellCycleModel());
           CellPtr p_cell_3d_2(new Cell(p_state, p_model_3d_2));
           p_cell_3d_2->GetCellData()->SetItem("oxygen", oxygen_concentration);
           p_cell_3d_2->SetCellProliferativeType(p_stem_type);

           for (unsigned i=0; i<53; i++)
           {
        	   p_simulation_time->IncrementTimeOneStep();
        	   TS_ASSERT_EQUALS(p_model_1d->ReadyToDivide(), false);
        	   TS_ASSERT_EQUALS(p_model_2d->ReadyToDivide(), false);
        	   TS_ASSERT_EQUALS(p_model_3d->ReadyToDivide(), false);
           }

           p_simulation_time->IncrementTimeOneStep();

           TS_ASSERT_EQUALS(p_model_1d->ReadyToDivide(), true);
           TS_ASSERT_EQUALS(p_model_2d->ReadyToDivide(), true);
           TS_ASSERT_EQUALS(p_model_3d->ReadyToDivide(), true);


           TS_ASSERT_THROWS_NOTHING(p_model_2d->ResetForDivision());
           TS_ASSERT_EQUALS(p_model_2d->ReadyToDivide(), false);

           for (unsigned i=0; i<53; i++)
           {
        	   p_simulation_time->IncrementTimeOneStep();
        	   TS_ASSERT_EQUALS(p_model_1d->ReadyToDivide(), true);
        	   TS_ASSERT_EQUALS(p_model_2d->ReadyToDivide(), false);
           }

           p_simulation_time->IncrementTimeOneStep();

           TS_ASSERT_EQUALS(p_model_1d->ReadyToDivide(), true);
           TS_ASSERT_EQUALS(p_model_2d->ReadyToDivide(), true);

           // For coverage, create a 1D model
           Owen2011OxygenBasedCellCycleModel* p_cell_model3 = new Owen2011OxygenBasedCellCycleModel();
           p_cell_model3->SetDimension(1);

           CellPtr p_cell3(new Cell(p_state, p_cell_model3));
           p_cell3->GetCellData()->SetItem("oxygen", oxygen_concentration);
           p_cell3->InitialiseCellCycleModel();
           p_cell3->SetCellProliferativeType(p_stem_type);

           TS_ASSERT_DELTA(p_cell_model3->GetProteinConcentrations()[1], 10.0, 1e-5);
           TS_ASSERT_EQUALS(p_cell_model3->ReadyToDivide(), false);
           p_simulation_time->IncrementTimeOneStep();
           TS_ASSERT_EQUALS(p_cell_model3->ReadyToDivide(), false);
	}

    void TestOwen2011OxygenBasedCellCycleModelWithOdeForQuiescentCancerCells() throw(Exception)
    {
         //Check that oxygen concentration is set up correctly
    	 SimulationTime* p_simulation_time = SimulationTime::Instance();
         p_simulation_time->SetEndTimeAndNumberOfTimeSteps(3.0, 3);

         Owen2011OxygenBasedCellCycleModel* p_model1 = new Owen2011OxygenBasedCellCycleModel();
         p_model1->SetDimension(2);

         MAKE_PTR(CancerCellMutationState, p_state);
     	 MAKE_PTR(StemCellProliferativeType, p_stem_type);
         CellPtr p_cell1(new Cell(p_state, p_model1));

    	 double lo_oxygen = 1.0;
         double hi_oxygen = 10.0;

         p_cell1->GetCellData()->SetItem("oxygen", lo_oxygen);
    	 p_cell1->SetCellProliferativeType(p_stem_type);
         p_cell1->InitialiseCellCycleModel();

    	 p_model1->ReadyToDivide();
    	 TS_ASSERT_DELTA(p_model1->GetCurrentQuiescentDuration(), 0.0, 1e-12);
    	 TS_ASSERT_DELTA(p_model1->GetCurrentQuiescenceOnsetTime(), 0.0, 1e-12);

    	 p_simulation_time->IncrementTimeOneStep(); // t=1.0
    	 p_model1->ReadyToDivide();
    	 TS_ASSERT_DELTA(p_model1->GetCurrentQuiescentDuration(), 1.0, 1e-12);
    	 TS_ASSERT_DELTA(p_model1->GetCurrentQuiescenceOnsetTime(), 0.0, 1e-12);

    	 p_cell1->GetCellData()->SetItem("oxygen", hi_oxygen);
    	 p_simulation_time->IncrementTimeOneStep(); // t=2.0
    	 p_model1->ReadyToDivide();
    	 TS_ASSERT_DELTA(p_model1->GetCurrentQuiescentDuration(), 0.0, 1e-12);
    	 TS_ASSERT_DELTA(p_model1->GetCurrentQuiescenceOnsetTime(), 0.0, 1e-12);

    	 p_cell1->GetCellData()->SetItem("oxygen", lo_oxygen);
    	 p_simulation_time->IncrementTimeOneStep(); // t=3.0
    	 p_model1->ReadyToDivide();
    	 TS_ASSERT_DELTA(p_model1->GetCurrentQuiescentDuration(), 0.0, 1e-12);
    	 TS_ASSERT_DELTA(p_model1->GetCurrentQuiescenceOnsetTime(), 3.0, 1e-12);

    	 // Set up SimulationTime
    	 SimulationTime::Destroy();
    	 p_simulation_time = SimulationTime::Instance();
    	 p_simulation_time->SetStartTime(0.0);
    	 p_simulation_time->SetEndTimeAndNumberOfTimeSteps(75.0, 150);

        // Create cell-cycle models and cells
        Owen2011OxygenBasedCellCycleModel* p_model = new Owen2011OxygenBasedCellCycleModel();
    	p_model->SetDimension(2);
    	p_model->SetBirthTime(0.0);

        CellPtr p_cell(new Cell(p_state, p_model));
        p_cell->GetCellData()->SetItem("oxygen", hi_oxygen);
   	    p_cell->SetCellProliferativeType(p_stem_type);
        p_cell->InitialiseCellCycleModel();

    	// Check that the cell cycle phase and ready to divide are updated correctly
        TS_ASSERT_EQUALS(p_model->ReadyToDivide(), false);
        TS_ASSERT_EQUALS(p_model->GetCurrentCellCyclePhase(),G_ONE_PHASE);

    	for (unsigned i=0; i<53; i++)
    	{
    	  p_simulation_time->IncrementTimeOneStep();

          // Note that we need to pass in the updated G1 duration
          CheckReadyToDivideAndPhaseIsUpdated(p_model, 27);
        }

    	TS_ASSERT_DELTA(p_model->GetAge(), p_simulation_time->GetTime(), 1e-9);

    	// Check that cell division correctly resets the cell cycle phase
        p_simulation_time->IncrementTimeOneStep();
        TS_ASSERT_EQUALS(p_cell->ReadyToDivide(), true);
        CellPtr p_cell2 = p_cell->Divide();
        Owen2011OxygenBasedCellCycleModel* p_model2 = static_cast <Owen2011OxygenBasedCellCycleModel*>(p_cell2->GetCellCycleModel());

        TS_ASSERT_EQUALS(p_model2->GetCurrentCellCyclePhase(), M_PHASE);
        TS_ASSERT_EQUALS(p_model2->ReadyToDivide(), false);
        TS_ASSERT_EQUALS(p_model2->GetCurrentCellCyclePhase(), G_ONE_PHASE);

        p_cell->GetCellData()->SetItem("oxygen", lo_oxygen);

        TS_ASSERT_EQUALS(p_model->ReadyToDivide(), false);
        TS_ASSERT_EQUALS(p_model->GetCurrentCellCyclePhase(), G_ZERO_PHASE);

        p_simulation_time->IncrementTimeOneStep();
        TS_ASSERT_EQUALS(p_model->ReadyToDivide(), false);
    	TS_ASSERT_EQUALS(p_model->GetCurrentCellCyclePhase(),G_ZERO_PHASE);

    	p_simulation_time->IncrementTimeOneStep();
    	TS_ASSERT_EQUALS(p_model->ReadyToDivide(), false);
    	TS_ASSERT_EQUALS(p_model->GetCurrentCellCyclePhase(),G_ZERO_PHASE);

    	p_cell->GetCellData()->SetItem("oxygen", hi_oxygen);

    	TS_ASSERT_EQUALS(p_model->ReadyToDivide(), false);
    	TS_ASSERT_EQUALS(p_model->GetCurrentCellCyclePhase(),G_ONE_PHASE);

    	for (unsigned i=0; i<54; i++)
    	{
    	  p_simulation_time->IncrementTimeOneStep();

          // Note that we need to pass in the updated G1 duration
          CheckReadyToDivideAndPhaseIsUpdated(p_model, 27);
        }

       	// Check that cell division correctly resets the cell cycle phase
       	TS_ASSERT_EQUALS(p_cell->ReadyToDivide(), true);
        CellPtr p_cell3 = p_cell->Divide();
        Owen2011OxygenBasedCellCycleModel* p_model3 = static_cast <Owen2011OxygenBasedCellCycleModel*>(p_cell3->GetCellCycleModel());

        TS_ASSERT_EQUALS(p_model3->ReadyToDivide(), false);
        TS_ASSERT_EQUALS(p_model3->GetCurrentCellCyclePhase(), G_ONE_PHASE);

       	// For coverage, create a 1D model

        Owen2011OxygenBasedCellCycleModel* p_cell_model1d = new Owen2011OxygenBasedCellCycleModel();
    	p_cell_model1d->SetDimension(1);

    	CellPtr p_cell1d(new Cell(p_state, p_cell_model1d));
    	p_cell1d->GetCellData()->SetItem("oxygen", hi_oxygen);
   	    p_cell1d->SetCellProliferativeType(p_stem_type);
        p_cell1d->InitialiseCellCycleModel();

    	TS_ASSERT_EQUALS(p_cell_model1d->ReadyToDivide(), false);

        // For coverage, create a 3D model

        Owen2011OxygenBasedCellCycleModel* p_cell_model3d = new Owen2011OxygenBasedCellCycleModel();
    	p_cell_model3d->SetDimension(3);

    	CellPtr p_cell3d(new Cell(p_state, p_cell_model3d));
    	p_cell3d->GetCellData()->SetItem("oxygen", hi_oxygen);
   	    p_cell3d->SetCellProliferativeType(p_stem_type);
        p_cell3d->InitialiseCellCycleModel();

        TS_ASSERT_EQUALS(p_cell_model3d->ReadyToDivide(), false);
    }

    void TestArchiveOwen2011OxygenBasedCellCycleModelWithOdeForNormalCells()
    {
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "Owen2011OxygenBasedCellCycleModel.arch";

        // Set up oxygen_concentration
        double oxygen_concentration = 1.0;

        {
            // We must set up SimulationTime to avoid memory leaks
            SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);

            // As usual, we archive via a pointer to the most abstract class possible
            AbstractCellCycleModel* const p_model = new Owen2011OxygenBasedCellCycleModel;

            p_model->SetDimension(1);
            p_model->SetBirthTime(-1.5);
            static_cast<Owen2011OxygenBasedCellCycleModel*>(p_model)->SetDt(0.5);

            // We must create a cell to be able to initialise the cell cycle model's ODE system
            MAKE_PTR(WildTypeCellMutationState, p_healthy_state);
        	MAKE_PTR(StemCellProliferativeType, p_stem_type);
            CellPtr p_cell(new Cell(p_healthy_state, p_model));
            p_cell->GetCellData()->SetItem("oxygen", oxygen_concentration);
       	    p_cell->SetCellProliferativeType(p_stem_type);
            p_cell->InitialiseCellCycleModel();

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            output_arch << p_model;

            // Note that here, deletion of the cell-cycle model is handled by the cell destructor
            SimulationTime::Destroy();
        }

        {
            // We must set SimulationTime::mStartTime here to avoid tripping an assertion
            SimulationTime::Instance()->SetStartTime(0.0);

            AbstractCellCycleModel* p_model2;

            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            input_arch >> p_model2;

            TS_ASSERT_EQUALS(p_model2->GetDimension(), 1u);
            TS_ASSERT_DELTA(p_model2->GetBirthTime(), -1.5, 1e-12);
            TS_ASSERT_DELTA(static_cast<Owen2011OxygenBasedCellCycleModel*>(p_model2)->GetDt(), 0.5, 1e-3);

            Owen2011OxygenBasedCellCycleModel* p_static_cast_model =
                static_cast<Owen2011OxygenBasedCellCycleModel*>(p_model2);

            Owen2011OxygenBasedCellCycleOdeSystem* p_ode_system =
                static_cast<Owen2011OxygenBasedCellCycleOdeSystem*>(p_static_cast_model->GetOdeSystem());

            TS_ASSERT(p_ode_system != NULL);

            // Avoid memory leaks
            delete p_model2;
        }
    }

    void TestArchiveOwen2011OxygenBasedCellCycleModelWithOdeForCancerCells()
        {
            OutputFileHandler handler("archive", false);
            std::string archive_filename = handler.GetOutputDirectoryFullPath() + "Owen2011OxygenBasedCellCycleModel.arch";

            // Set up oxygen_concentration
            double oxygen_concentration = 1.0;

            {
                // We must set up SimulationTime to avoid memory leaks
                SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);

                // As usual, we archive via a pointer to the most abstract class possible
                AbstractCellCycleModel* const p_model = new Owen2011OxygenBasedCellCycleModel;

                p_model->SetDimension(1);
                p_model->SetBirthTime(-1.5);
                static_cast<Owen2011OxygenBasedCellCycleModel*>(p_model)->SetDt(0.5);

                // We must create a cell to be able to initialise the cell cycle model's ODE system
                MAKE_PTR(CancerCellMutationState, p_cancer_state);
            	MAKE_PTR(StemCellProliferativeType, p_stem_type);
                CellPtr p_cell(new Cell(p_cancer_state, p_model));
                p_cell->GetCellData()->SetItem("oxygen", oxygen_concentration);
           	    p_cell->SetCellProliferativeType(p_stem_type);
                p_cell->InitialiseCellCycleModel();

                std::ofstream ofs(archive_filename.c_str());
                boost::archive::text_oarchive output_arch(ofs);

                output_arch << p_model;

                // Note that here, deletion of the cell-cycle model is handled by the cell destructor
                SimulationTime::Destroy();
            }

            {
                // We must set SimulationTime::mStartTime here to avoid tripping an assertion
                SimulationTime::Instance()->SetStartTime(0.0);

                AbstractCellCycleModel* p_model2;

                std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
                boost::archive::text_iarchive input_arch(ifs);

                input_arch >> p_model2;

                TS_ASSERT_EQUALS(p_model2->GetDimension(), 1u);
                TS_ASSERT_DELTA(p_model2->GetBirthTime(), -1.5, 1e-12);
                TS_ASSERT_DELTA(static_cast<Owen2011OxygenBasedCellCycleModel*>(p_model2)->GetDt(), 0.5, 1e-3);

                Owen2011OxygenBasedCellCycleModel* p_static_cast_model =
                    static_cast<Owen2011OxygenBasedCellCycleModel*>(p_model2);

                Owen2011OxygenBasedCellCycleOdeSystem* p_ode_system =
                    static_cast<Owen2011OxygenBasedCellCycleOdeSystem*>(p_static_cast_model->GetOdeSystem());

                TS_ASSERT(p_ode_system != NULL);

                // Avoid memory leaks
                delete p_model2;
            }
        }

      ///\todo TestArchiveOwen2011OxygenBasedCellCycleModelForQuiescentCancerCells()
    void TestCellCycleModelOutputParameters()
    {
        std::string output_directory = "TestCellCycleModelOutputParameters";
        OutputFileHandler output_file_handler(output_directory, false);

        // Test with Owen2011OxygenBasedCellCycleModel
        Owen2011OxygenBasedCellCycleModel owen_oxygen_based_cell_cycle_model;
        TS_ASSERT_EQUALS(owen_oxygen_based_cell_cycle_model.GetIdentifier(), "Owen2011OxygenBasedCellCycleModel");

        out_stream owen_oxygen_based_parameter_file = output_file_handler.OpenOutputFile("owen_oxygen_based_results.parameters");
        owen_oxygen_based_cell_cycle_model.OutputCellCycleModelParameters(owen_oxygen_based_parameter_file);
        owen_oxygen_based_parameter_file->close();

        std::string owen_oxygen_based_results_dir = output_file_handler.GetOutputDirectoryFullPath();
        TS_ASSERT_EQUALS(system(("diff " + owen_oxygen_based_results_dir + "owen_oxygen_based_results.parameters projects/VascularTumour/test/data/TestCellCycleModel/owen_oxygen_based_results.parameters").c_str()), 0);

     }
};

#endif /*TESTODEBASEDCELLCYCLEMODELS_HPP_*/
