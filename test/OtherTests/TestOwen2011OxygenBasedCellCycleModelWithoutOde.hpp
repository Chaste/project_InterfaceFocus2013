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

#ifndef TESTOWEN2011OXYGENBASEDCELLCYCLEMODELWITHOUTODE_HPP_
#define TESTOWEN2011OXYGENBASEDCELLCYCLEMODELWITHOUTODE_HPP_


#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <fstream>
#include <boost/shared_ptr.hpp>

#include "Owen2011OxygenBasedCellCycleModelWithoutOde.hpp"

#include "AbstractCellMutationState.hpp"
#include "WildTypeCellMutationState.hpp"
#include "CancerCellMutationState.hpp"
#include "StemCellProliferativeType.hpp"

#include "OutputFileHandler.hpp"
#include "CheckReadyToDivideAndPhaseIsUpdated.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "SmartPointers.hpp"

class TestOwen2011OxygenBasedCellCycleModelWithoutOde : public AbstractCellBasedTestSuite
{
public:

    void TestOwen2011CellCycleModelWithoutOdeForNormalCells() throw(Exception)
{
        SimulationTime* p_simulation_time = SimulationTime::Instance();
        p_simulation_time->SetEndTimeAndNumberOfTimeSteps(75.0, 150);

        // Create cell-cycle models
        Owen2011OxygenBasedCellCycleModelWithoutOde* p_model_1d = new Owen2011OxygenBasedCellCycleModelWithoutOde();
        p_model_1d->SetDimension(1);

        Owen2011OxygenBasedCellCycleModelWithoutOde* p_model_2d = new Owen2011OxygenBasedCellCycleModelWithoutOde();
        p_model_2d->SetDimension(2);

        Owen2011OxygenBasedCellCycleModelWithoutOde* p_model_3d = new Owen2011OxygenBasedCellCycleModelWithoutOde();
        p_model_3d->SetDimension(3);

        // Set up oxygen_concentration
        double oxygen_concentration = 1.0;

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

        TS_ASSERT_DELTA(p_model_2d->GetPhi(), 0.00000000, 1e-5);
        TS_ASSERT_EQUALS(p_model_2d->ReadyToDivide(), false);

        // Divide the cells
        Owen2011OxygenBasedCellCycleModelWithoutOde* p_model_1d_2 = static_cast<Owen2011OxygenBasedCellCycleModelWithoutOde*> (p_model_1d->CreateCellCycleModel());
        CellPtr p_cell_1d_2(new Cell(p_state, p_model_1d_2));
        p_cell_1d_2->GetCellData()->SetItem("oxygen", oxygen_concentration);
        p_cell_1d_2->SetCellProliferativeType(p_stem_type);

        Owen2011OxygenBasedCellCycleModelWithoutOde* p_model_2d_2 = static_cast<Owen2011OxygenBasedCellCycleModelWithoutOde*> (p_model_2d->CreateCellCycleModel());
        CellPtr p_cell_2d_2(new Cell(p_state, p_model_2d_2));
        p_cell_2d_2->GetCellData()->SetItem("oxygen", oxygen_concentration);
        p_cell_2d_2->SetCellProliferativeType(p_stem_type);

        Owen2011OxygenBasedCellCycleModelWithoutOde* p_model_3d_2 = static_cast<Owen2011OxygenBasedCellCycleModelWithoutOde*> (p_model_3d->CreateCellCycleModel());
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
        TS_ASSERT_DELTA(p_model_1d->GetPhi(), 1.00000000, 1e-2);
        TS_ASSERT_EQUALS(p_model_2d->ReadyToDivide(), true);
        TS_ASSERT_EQUALS(p_model_3d->ReadyToDivide(), true);

        TS_ASSERT_THROWS_NOTHING(p_model_2d->ResetForDivision());
        TS_ASSERT_DELTA(p_model_2d->GetPhi(), 0.00000000, 1e-5);

        // For coverage, create a 1D model
        Owen2011OxygenBasedCellCycleModelWithoutOde* p_cell_model3 = new Owen2011OxygenBasedCellCycleModelWithoutOde();
        p_cell_model3->SetDimension(1);
        CellPtr p_cell3(new Cell(p_state, p_cell_model3));
        p_cell3->GetCellData()->SetItem("oxygen", oxygen_concentration);
        p_cell3->SetCellProliferativeType(p_stem_type);
        p_cell3->InitialiseCellCycleModel();

        TS_ASSERT_EQUALS(p_cell_model3->ReadyToDivide(), false);
        p_simulation_time->IncrementTimeOneStep();
        TS_ASSERT_EQUALS(p_cell_model3->ReadyToDivide(), false);
}

    void TestOwen2011CellCycleModelWithoutOdeForCancerCells() throw(Exception)
{
          SimulationTime* p_simulation_time = SimulationTime::Instance();
          p_simulation_time->SetEndTimeAndNumberOfTimeSteps(75.0, 150);

          // Create cell-cycle models
          Owen2011OxygenBasedCellCycleModelWithoutOde* p_model_1d = new Owen2011OxygenBasedCellCycleModelWithoutOde();
          p_model_1d->SetDimension(1);

          Owen2011OxygenBasedCellCycleModelWithoutOde* p_model_2d = new Owen2011OxygenBasedCellCycleModelWithoutOde();
          p_model_2d->SetDimension(2);

          Owen2011OxygenBasedCellCycleModelWithoutOde* p_model_3d = new Owen2011OxygenBasedCellCycleModelWithoutOde();
          p_model_3d->SetDimension(3);

          // Set up oxygen_concentration
          double oxygen_concentration = 10.0;

          // Create cells
          MAKE_PTR(CancerCellMutationState, p_state);
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

          TS_ASSERT_DELTA(p_model_2d->GetPhi(), 0.00000000, 1e-5);
          TS_ASSERT_EQUALS(p_model_2d->ReadyToDivide(), false);

          // Divide the cells
          Owen2011OxygenBasedCellCycleModelWithoutOde* p_model_1d_2 = static_cast<Owen2011OxygenBasedCellCycleModelWithoutOde*> (p_model_1d->CreateCellCycleModel());
          CellPtr p_cell_1d_2(new Cell(p_state, p_model_1d_2));
          p_cell_1d_2->GetCellData()->SetItem("oxygen", oxygen_concentration);
          p_cell_1d_2->SetCellProliferativeType(p_stem_type);

          Owen2011OxygenBasedCellCycleModelWithoutOde* p_model_2d_2 = static_cast<Owen2011OxygenBasedCellCycleModelWithoutOde*> (p_model_2d->CreateCellCycleModel());
          CellPtr p_cell_2d_2(new Cell(p_state, p_model_2d_2));
          p_cell_2d_2->GetCellData()->SetItem("oxygen", oxygen_concentration);
          p_cell_2d_2->SetCellProliferativeType(p_stem_type);

          Owen2011OxygenBasedCellCycleModelWithoutOde* p_model_3d_2 = static_cast<Owen2011OxygenBasedCellCycleModelWithoutOde*> (p_model_3d->CreateCellCycleModel());
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
          TS_ASSERT_DELTA(p_model_1d->GetPhi(), 1.00000000, 1e-1);
          TS_ASSERT_EQUALS(p_model_2d->ReadyToDivide(), true);
          TS_ASSERT_EQUALS(p_model_3d->ReadyToDivide(), true);

          TS_ASSERT_THROWS_NOTHING(p_model_2d->ResetForDivision());
          TS_ASSERT_EQUALS(p_model_2d->ReadyToDivide(), false);
          TS_ASSERT_DELTA(p_model_2d->GetPhi(), 0.00000000, 1e-5);

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
          Owen2011OxygenBasedCellCycleModelWithoutOde* p_cell_model3 = new Owen2011OxygenBasedCellCycleModelWithoutOde();
          p_cell_model3->SetDimension(1);

          CellPtr p_cell3(new Cell(p_state, p_cell_model3));
          p_cell3->GetCellData()->SetItem("oxygen", oxygen_concentration);
          p_cell3->SetCellProliferativeType(p_stem_type);
          p_cell3->InitialiseCellCycleModel();

          TS_ASSERT_EQUALS(p_cell_model3->ReadyToDivide(), false);
          p_simulation_time->IncrementTimeOneStep();
          TS_ASSERT_EQUALS(p_cell_model3->ReadyToDivide(), false);
}

    void TestOwen2011OxygenBasedCellCycleModelWithoutOdeForQuiescentCancerCells() throw(Exception)
{
        // Check that mCurrentQuiescenceOnsetTime and mCurrentQuiescentDuration are updated correctly
    	SimulationTime* p_simulation_time = SimulationTime::Instance();
    	p_simulation_time->SetEndTimeAndNumberOfTimeSteps(3.0, 3);

    	Owen2011OxygenBasedCellCycleModelWithoutOde* p_model1 = new Owen2011OxygenBasedCellCycleModelWithoutOde();
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
        Owen2011OxygenBasedCellCycleModelWithoutOde* p_model = new Owen2011OxygenBasedCellCycleModelWithoutOde();
    	p_model->SetDimension(2);
    	p_model->SetBirthTime(0.0);

        CellPtr p_cell(new Cell(p_state, p_model));
        p_cell->GetCellData()->SetItem("oxygen", hi_oxygen);
        p_cell->SetCellProliferativeType(p_stem_type);
        p_cell->InitialiseCellCycleModel();

    	// Check that the cell cycle phase and ready to divide are updated correctly
        TS_ASSERT_EQUALS(p_model->ReadyToDivide(), false);
        TS_ASSERT_EQUALS(p_model->GetCurrentCellCyclePhase(),G_ONE_PHASE);
        TS_ASSERT_LESS_THAN(p_model->GetPhi(),1.0);

    	for (unsigned i=0; i<53; i++)
    	{
    	  p_simulation_time->IncrementTimeOneStep();
    	  CheckReadyToDivideAndPhaseIsUpdated(p_model, 26.675);
    	  TS_ASSERT_LESS_THAN(p_model->GetPhi(), 1.0);
        }

    	TS_ASSERT_DELTA(p_model->GetAge(), p_simulation_time->GetTime(), 1e-9);

    	// Check that cell division correctly resets the cell cycle phase
        p_simulation_time->IncrementTimeOneStep();
    	TS_ASSERT_EQUALS(p_cell->ReadyToDivide(), true);
        TS_ASSERT_LESS_THAN_EQUALS(1.0, p_model->GetPhi());
        CellPtr p_cell2 = p_cell->Divide();
        Owen2011OxygenBasedCellCycleModelWithoutOde* p_model2 = static_cast <Owen2011OxygenBasedCellCycleModelWithoutOde*>(p_cell2->GetCellCycleModel());

        TS_ASSERT_EQUALS(p_model2->GetCurrentCellCyclePhase(), G_ONE_PHASE);
        TS_ASSERT_EQUALS(p_model2->ReadyToDivide(), false);
        TS_ASSERT_LESS_THAN(p_model2->GetPhi(), 1.0);
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
        Owen2011OxygenBasedCellCycleModelWithoutOde* p_model3 = static_cast <Owen2011OxygenBasedCellCycleModelWithoutOde*>(p_cell3->GetCellCycleModel());

        TS_ASSERT_EQUALS(p_model3->ReadyToDivide(), false);
        TS_ASSERT_LESS_THAN(p_model->GetPhi(), 1.0);
        TS_ASSERT_EQUALS(p_model2->GetCurrentCellCyclePhase(), G_ONE_PHASE);

       	// For coverage, create a 1D model

        Owen2011OxygenBasedCellCycleModelWithoutOde* p_cell_model1d = new Owen2011OxygenBasedCellCycleModelWithoutOde();
    	p_cell_model1d->SetDimension(1);

    	CellPtr p_cell1d(new Cell(p_state, p_cell_model1d));
        p_cell1d->GetCellData()->SetItem("oxygen", hi_oxygen);
        p_cell1d->SetCellProliferativeType(p_stem_type);
        p_cell1d->InitialiseCellCycleModel();

    	TS_ASSERT_EQUALS(p_cell_model1d->ReadyToDivide(), false);

        // For coverage, create a 3D model

        Owen2011OxygenBasedCellCycleModelWithoutOde* p_cell_model3d = new Owen2011OxygenBasedCellCycleModelWithoutOde();
    	p_cell_model3d->SetDimension(3);

    	CellPtr p_cell3d(new Cell(p_state, p_cell_model3d));
        p_cell3d->GetCellData()->SetItem("oxygen", hi_oxygen);
        p_cell3d->SetCellProliferativeType(p_stem_type);
        p_cell3d->InitialiseCellCycleModel();

        TS_ASSERT_EQUALS(p_cell_model3d->ReadyToDivide(), false);
        TS_ASSERT_LESS_THAN(p_model->GetPhi(), 1.0);
    }

      ///\todo TestArchiveOwen2011OxygenBasedCellCycleModeWithoutOde()
//    void TestArchiveSimpleOxygenBasedCellCycleModel() throw (Exception)
//    {
//        OutputFileHandler handler("archive", false);
//        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "SimpleOxygenBasedCellCycleModel.arch";
//
//        std::vector<double> oxygen_concentration;
//        oxygen_concentration.push_back(1.0);
//        CellwiseData<3>::Instance()->SetConstantDataForTesting(oxygen_concentration);
//
//        {
//            // We must set up SimulationTime to avoid memory leaks
//            SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(1.0, 1);
//
//            // As usual, we archive via a pointer to the most abstract class possible
//            AbstractCellCycleModel* const p_model = new SimpleOxygenBasedCellCycleModel;
//            p_model->SetDimension(3);
//            p_model->SetCellProliferativeType(STEM);
//
//            static_cast<SimpleOxygenBasedCellCycleModel*>(p_model)->SetHypoxicConcentration(0.8);
//            static_cast<SimpleOxygenBasedCellCycleModel*>(p_model)->SetQuiescentConcentration(0.7);
//            static_cast<SimpleOxygenBasedCellCycleModel*>(p_model)->SetCriticalHypoxicDuration(2.5);
//            static_cast<SimpleOxygenBasedCellCycleModel*>(p_model)->SetCurrentHypoxiaOnsetTime(3.1);
//
//            std::ofstream ofs(archive_filename.c_str());
//            boost::archive::text_oarchive output_arch(ofs);
//
//            output_arch << p_model;
//
//            delete p_model;
//            SimulationTime::Destroy();
//        }
//
//        {
//            // We must set SimulationTime::mStartTime here to avoid tripping an assertion
//            SimulationTime::Instance()->SetStartTime(0.0);
//
//            AbstractCellCycleModel* p_model2;
//
//            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
//            boost::archive::text_iarchive input_arch(ifs);
//
//            input_arch >> p_model2;
//
//            TS_ASSERT_DELTA(static_cast<SimpleOxygenBasedCellCycleModel*>(p_model2)->GetHypoxicConcentration(), 0.8, 1e-6);
//            TS_ASSERT_DELTA(static_cast<SimpleOxygenBasedCellCycleModel*>(p_model2)->GetQuiescentConcentration(), 0.7, 1e-6);
//            TS_ASSERT_DELTA(static_cast<SimpleOxygenBasedCellCycleModel*>(p_model2)->GetCriticalHypoxicDuration(), 2.5, 1e-6);
//            TS_ASSERT_DELTA(static_cast<SimpleOxygenBasedCellCycleModel*>(p_model2)->GetCurrentHypoxiaOnsetTime(), 3.1, 1e-6);
//
//            // Avoid memory leaks
//            delete p_model2;
//            CellwiseData<3>::Destroy();
//        }
//}
//
//
//    void TestCellCycleModelOutputParameters()
//    {
//        std::string output_directory = "TestCellCycleModelOutputParameters";
//        OutputFileHandler output_file_handler(output_directory, false);
//
//
//        // Test with SimpleOxygenBasedCellCycleModel
//        SimpleOxygenBasedCellCycleModel simple_oxygen_based_cell_cycle_model;
//        TS_ASSERT_EQUALS(simple_oxygen_based_cell_cycle_model.GetIdentifier(), "SimpleOxygenBasedCellCycleModel");
//
//        out_stream simple_oxygen_based_parameter_file = output_file_handler.OpenOutputFile("simple_oxygen_based_results.parameters");
//        simple_oxygen_based_cell_cycle_model.OutputCellCycleModelParameters(simple_oxygen_based_parameter_file);
//        simple_oxygen_based_parameter_file->close();
//
//        std::string simple_oxygen_based_results_dir = output_file_handler.GetOutputDirectoryFullPath();
//        TS_ASSERT_EQUALS(system(("diff " + simple_oxygen_based_results_dir + "simple_oxygen_based_results.parameters cell_based/test/data/TestCellCycleModels/simple_oxygen_based_results.parameters").c_str()), 0);
//
//     }
};

#endif /*TESTOWEN2011OXYGENBASEDCELLCYCLEMODELWITHOUTODE_HPP_*/
