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

#ifndef TESTOWEN2011OXYGENBASEDCELLCYCLEODESYSTEM_HPP_
#define TESTOWEN2011OXYGENBASEDCELLCYCLEODESYSTEM_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <stdio.h>
#include <ctime>
#include <vector>
#include <iostream>

#include "OutputFileHandler.hpp"

#include "Owen2011OxygenBasedCellCycleOdeSystem.hpp"
#include "RungeKutta4IvpOdeSolver.hpp"
#include "RungeKuttaFehlbergIvpOdeSolver.hpp"
#include "BackwardEulerIvpOdeSolver.hpp"

/**
 * This class contains tests for Owen2011OxygenBasedCellCycleOdeSystem,
 * a system of ODEs that are used by the cell-cycle model
 * Owen2011OxygenBasedCellCycleModel to determine when a cell is ready
 * to divide.
 */
class TestOwen2011OxygenBasedCellCycleOdeSystem : public CxxTest::TestSuite
{
public:

    /**
     * Test derivative calculations (correct values calculated using Matlab).
     */
    void TestOwen2011Equations()
    {
        // Set up
        double time = 0.0;
        double oxygen_concentration = 1.0;

        Owen2011OxygenBasedCellCycleOdeSystem normal_system(oxygen_concentration, false);
//      Owen2011OxygenBasedCellCycleOdeSystem cancer_system(oxygen_concentration, true);

        std::vector<double> initial_conditions = normal_system.GetInitialConditions();

        std::vector<double> normal_derivs(initial_conditions.size());
//        std::vector<double> cancer_derivs(initial_conditions.size());
        normal_system.EvaluateYDerivatives(time, initial_conditions, normal_derivs);
//        cancer_system.EvaluateYDerivatives(time, initial_conditions, cancer_derivs);

        /**
         * Test derivatives are correct initially
         * (correct values calculated using Matlab code)
         */

        // Normal cell
        TS_ASSERT_DELTA(normal_derivs[0], 0.0198659050, 1e-5);
        TS_ASSERT_DELTA(normal_derivs[1], 0.0000000000, 1e-5);

//        // Cancer cell
//        TS_ASSERT_DELTA(cancer_derivs[0], 455.630699088, 1e-5);
//        TS_ASSERT_DELTA(cancer_derivs[1], 1.83600000000, 1e-5);
//        TS_ASSERT_DELTA(cancer_derivs[2], 0.42000000000, 1e-5);
//        TS_ASSERT_DELTA(cancer_derivs[3], 1.50000000000, 1e-5);
//        TS_ASSERT_DELTA(cancer_derivs[4], -5.4060000000, 1e-5);
//        TS_ASSERT_DELTA(cancer_derivs[5], 0.00000000000, 1e-5);

        /**
         * Again test derivatives are correct initially, but for
         * different initial conditions (corresponding to a low
         * oxygen concentration). The usual initial condition for
         * z is zero, so we need to change it to see any difference.
         */
        oxygen_concentration = 0.1;

        Owen2011OxygenBasedCellCycleOdeSystem normal_system2(oxygen_concentration, false);
//      Owen2011OxygenBasedCellCycleOdeSystem cancer_system2(oxygen_concentration, true);

        std::vector<double> normal_derivs2(initial_conditions.size());
        normal_system2.SetDefaultInitialCondition(1, 0.1);

//        std::vector<double> cancer_derivs2(initial_conditions.size());
//        cancer_system2.SetDefaultInitialCondition(2, 0.1);

        std::vector<double> initial_conditions2 = normal_system2.GetInitialConditions();

        normal_system2.EvaluateYDerivatives(time, initial_conditions2, normal_derivs2);
//        cancer_system2.EvaluateYDerivatives(time, initial_conditions2, cancer_derivs2);

        // Normal cell
        TS_ASSERT_DELTA(normal_derivs2[0], 0.0187353630, 1e-5);
        TS_ASSERT_DELTA(normal_derivs2[1], 0.0000000000, 1e-5);

//        // Cancer cell
//        TS_ASSERT_DELTA(cancer_derivs2[0], 455.630699088, 1e-5);
//        TS_ASSERT_DELTA(cancer_derivs2[1], 1.82100000000, 1e-5);
//        TS_ASSERT_DELTA(cancer_derivs2[2], 0.36545454545, 1e-5);
//        TS_ASSERT_DELTA(cancer_derivs2[3], 1.50000000000, 1e-5);
//        TS_ASSERT_DELTA(cancer_derivs2[4], -5.4060000000, 1e-5);
//        TS_ASSERT_DELTA(cancer_derivs2[5], 0.00000000000, 1e-5);
    }

    /**
     * Test two ODE solvers with this ODE system (correct values calculated using the Matlab solver ode15s).
     *
     */
    void TestOwen2011Solver() throw(Exception)
    {
        // Set up
        double oxygen_concentration = 1.0;
        Owen2011OxygenBasedCellCycleOdeSystem owen_system(oxygen_concentration, false);

        // Create ODE solvers
        RungeKutta4IvpOdeSolver rk4_solver;
        RungeKuttaFehlbergIvpOdeSolver rkf_solver;
        BackwardEulerIvpOdeSolver back_solver(2);

        // Set up for solver
        OdeSolution solutions;
        std::vector<double> initial_conditions = owen_system.GetInitialConditions();
        double start_time = 0.0;
        double end_time = 0.0;
        double elapsed_time = 0.0;
        double h_value = 1e-4; // maximum tolerance

        // Solve the ODE system using a Runge Kutta fourth order solver
        start_time = std::clock();
        solutions = rk4_solver.Solve(&owen_system, initial_conditions, 0.0, 60.0, h_value, h_value);
        end_time = std::clock();
        elapsed_time = (end_time - start_time)/(CLOCKS_PER_SEC);
        std::cout << "1. Runge-Kutta Elapsed time = " << elapsed_time << "\n";

//        // Reset maximum tolerance for Runge Kutta Fehlber solver
//        h_value = 1e-1;

        // Solve the ODE system using a Runge Kutta Fehlber solver
        initial_conditions = owen_system.GetInitialConditions();
        start_time = std::clock();
        solutions = rkf_solver.Solve(&owen_system, initial_conditions, 0.0, 60.0, h_value, h_value);
        end_time = std::clock();
        elapsed_time = (end_time - start_time)/(CLOCKS_PER_SEC);
        std::cout << "2. Runge-Kutta-Fehlberg Elapsed time = " << elapsed_time << "\n";

        // Test that solutions are accurate for a small time increase
        int end = solutions.rGetSolutions().size() - 1;

        // Test that the solver stops at the right time
        TS_ASSERT_DELTA(solutions.rGetTimes()[end], 50.3375, 1e-2);


        // Test solution - note the high tolerances
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][0], 1.000000000000000, 1e-3);
        TS_ASSERT_DELTA(solutions.rGetSolutions()[end][1], 1.000000000000000, 1e-3);

    }

    void TestArchiving()
    {
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "owen_ode.arch";

        {
            double oxygen_concentration = 0.7;
            bool is_labelled = true;

            Owen2011OxygenBasedCellCycleOdeSystem ode_system(oxygen_concentration, is_labelled);

            TS_ASSERT_DELTA(ode_system.GetOxygenConcentration(), 0.70, 1e-6);
            TS_ASSERT_EQUALS(ode_system.IsLabelled(), true);

            std::vector<double> initial_conditions = ode_system.GetInitialConditions();
            TS_ASSERT_EQUALS(initial_conditions.size(), 2u);
            TS_ASSERT_DELTA(initial_conditions[0], 0.00, 1e-6);
            TS_ASSERT_DELTA(initial_conditions[1], 0.70, 1e-6);

            // Create an output archive
            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Archive ODE system
            AbstractOdeSystem* const p_const_ode_system = &ode_system;
            output_arch << p_const_ode_system;
        }

        {
            AbstractOdeSystem* p_ode_system;

            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            // Restore from the archive
            input_arch >> p_ode_system;

            // Check that archiving worked correctly
            TS_ASSERT_DELTA(static_cast<Owen2011OxygenBasedCellCycleOdeSystem*>(p_ode_system)->GetOxygenConcentration(), 0.70, 1e-6);
            TS_ASSERT_EQUALS(static_cast<Owen2011OxygenBasedCellCycleOdeSystem*>(p_ode_system)->IsLabelled(), true);

            std::vector<double> initial_conditions = p_ode_system->GetInitialConditions();
            TS_ASSERT_EQUALS(initial_conditions.size(), 2u);
            TS_ASSERT_DELTA(initial_conditions[0], 0.00, 1e-6);
            TS_ASSERT_DELTA(initial_conditions[1], 0.70, 1e-6);

            // Tidy up
            delete p_ode_system;
        }
    }
};

#endif /*TESTOWEN2011OXYGENBASEDCELLCYCLEODESYSTEM_HPP_*/
