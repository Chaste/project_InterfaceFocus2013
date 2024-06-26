/*

Copyright (c) 2005-2013, University of Oxford.
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

#include <cmath>
#include <ctime>
#include <iostream>
#include <fstream>
#include <set>

#include "AbstractCellBasedSimulationInterfaceFocus.hpp"
#include "CellBasedEventHandler.hpp"
#include "LogFile.hpp"
#include "Version.hpp"
#include "ExecutableSupport.hpp"
#include "Exception.hpp"
#include <typeinfo>

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
AbstractCellBasedSimulationInterfaceFocus<ELEMENT_DIM,SPACE_DIM>::AbstractCellBasedSimulationInterfaceFocus(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation,
                                              bool deleteCellPopulationInDestructor,
                                              bool initialiseCells)
    : mDt(DOUBLE_UNSET),
      mEndTime(DOUBLE_UNSET),  // hours - this is set later on
      mrCellPopulation(rCellPopulation),
      mDeleteCellPopulationInDestructor(deleteCellPopulationInDestructor),
      mInitialiseCells(initialiseCells),
      mNoBirth(false),
      mUpdateCellPopulation(true),
      mOutputDirectory(""),
      mSimulationOutputDirectory(mOutputDirectory),
      mNumBirths(0),
      mNumDeaths(0),
      mOutputDivisionLocations(false),
      mSamplingTimestepMultiple(1),
      mpCellBasedPdeHandler(NULL)
{
    // Set a random seed of 0 if it wasn't specified earlier
    RandomNumberGenerator::Instance();

    if (mInitialiseCells)
    {
        mrCellPopulation.InitialiseCells();
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
AbstractCellBasedSimulationInterfaceFocus<ELEMENT_DIM,SPACE_DIM>::~AbstractCellBasedSimulationInterfaceFocus()
{
    if (mDeleteCellPopulationInDestructor)
    {
        delete &mrCellPopulation;
        delete mpCellBasedPdeHandler;
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractCellBasedSimulationInterfaceFocus<ELEMENT_DIM,SPACE_DIM>::SetCellBasedPdeHandler(CellBasedPdeHandlerInterfaceFocus<SPACE_DIM>* pCellBasedPdeHandler)
{
    mpCellBasedPdeHandler = pCellBasedPdeHandler;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
CellBasedPdeHandlerInterfaceFocus<SPACE_DIM>* AbstractCellBasedSimulationInterfaceFocus<ELEMENT_DIM,SPACE_DIM>::GetCellBasedPdeHandler()
{
    return mpCellBasedPdeHandler;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned AbstractCellBasedSimulationInterfaceFocus<ELEMENT_DIM,SPACE_DIM>::DoCellBirth()
{
    if (mNoBirth)
    {
        return 0;
    }

    unsigned num_births_this_step = 0;


    // Iterate over all cells, seeing if each one can be divided

    std::vector<CellPtr> vec_CellPopulation;

    for (typename AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>::Iterator cell_iter = mrCellPopulation.Begin();
         cell_iter != mrCellPopulation.End();
         ++cell_iter)
    {
    	vec_CellPopulation.push_back(*cell_iter);
    }
    std::random_shuffle(vec_CellPopulation.begin(),vec_CellPopulation.end());

    for (unsigned i = 0; i < vec_CellPopulation.size(); i++)
    {
    	CellPtr cell_iter = vec_CellPopulation[i];

        // Check if this cell is ready to divide
        if (cell_iter->GetAge() > 0.0)
        {
            if (cell_iter->ReadyToDivide())
            {
                // Check if there is room into which the cell may divide
                if (mrCellPopulation.IsRoomToDivide(cell_iter))
                {
                    // Store age before division
                    double cell_age = cell_iter->GetAge();

                    // Create a new cell
                    CellPtr p_new_cell = cell_iter->Divide();

                    // Call method that determines how cell division occurs and returns a vector
                    c_vector<double, SPACE_DIM> new_location = CalculateCellDivisionVector(cell_iter);

                    // If required, output this location to file
                    if (mOutputDivisionLocations)
                    {
                        *mpDivisionLocationFile << SimulationTime::Instance()->GetTime() << "\t";
                        for (unsigned i=0; i<SPACE_DIM; i++)
                        {
                            *mpDivisionLocationFile << new_location[i] << "\t";
                        }
                        *mpDivisionLocationFile << "\t" << cell_age << "\n";
                    }

                    // Add new cell to the cell population
                    mrCellPopulation.AddCell(p_new_cell, new_location, cell_iter);

                    // Update counter
                    num_births_this_step++;
                }
            }
        }
    }
    return num_births_this_step;
}
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned AbstractCellBasedSimulationInterfaceFocus<ELEMENT_DIM,SPACE_DIM>::DoCellRemoval()
{
    unsigned num_deaths_this_step = 0;

    /*
     * This labels cells as dead or apoptosing. It does not actually remove the cells,
     * mrCellPopulation.RemoveDeadCells() needs to be called for this.
     */
    for (typename std::vector<boost::shared_ptr<AbstractCellKiller<SPACE_DIM> > >::iterator killer_iter = mCellKillers.begin();
         killer_iter != mCellKillers.end();
         ++killer_iter)
    {
        (*killer_iter)->CheckAndLabelCellsForApoptosisOrDeath();
    }

    num_deaths_this_step += mrCellPopulation.RemoveDeadCells();

    return num_deaths_this_step;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractCellBasedSimulationInterfaceFocus<ELEMENT_DIM,SPACE_DIM>::SetDt(double dt)
{
    assert(dt > 0);
    mDt = dt;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double AbstractCellBasedSimulationInterfaceFocus<ELEMENT_DIM,SPACE_DIM>::GetDt()
{
    return mDt;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned AbstractCellBasedSimulationInterfaceFocus<ELEMENT_DIM,SPACE_DIM>::GetNumBirths()
{
    return mNumBirths;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned AbstractCellBasedSimulationInterfaceFocus<ELEMENT_DIM,SPACE_DIM>::GetNumDeaths()
{
    return mNumDeaths;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractCellBasedSimulationInterfaceFocus<ELEMENT_DIM,SPACE_DIM>::SetEndTime(double endTime)
{
    assert(endTime > 0);
    mEndTime = endTime;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractCellBasedSimulationInterfaceFocus<ELEMENT_DIM,SPACE_DIM>::SetOutputDirectory(std::string outputDirectory)
{
    mOutputDirectory = outputDirectory;
    mSimulationOutputDirectory = mOutputDirectory;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::string AbstractCellBasedSimulationInterfaceFocus<ELEMENT_DIM,SPACE_DIM>::GetOutputDirectory()
{
    return mOutputDirectory;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractCellBasedSimulationInterfaceFocus<ELEMENT_DIM,SPACE_DIM>::SetSamplingTimestepMultiple(unsigned samplingTimestepMultiple)
{
    assert(samplingTimestepMultiple > 0);
    mSamplingTimestepMultiple = samplingTimestepMultiple;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& AbstractCellBasedSimulationInterfaceFocus<ELEMENT_DIM,SPACE_DIM>::rGetCellPopulation()
{
    return mrCellPopulation;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
const AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& AbstractCellBasedSimulationInterfaceFocus<ELEMENT_DIM,SPACE_DIM>::rGetCellPopulation() const
{
    return mrCellPopulation;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractCellBasedSimulationInterfaceFocus<ELEMENT_DIM,SPACE_DIM>::SetUpdateCellPopulationRule(bool updateCellPopulation)
{
    mUpdateCellPopulation = updateCellPopulation;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool AbstractCellBasedSimulationInterfaceFocus<ELEMENT_DIM,SPACE_DIM>::GetUpdateCellPopulationRule()
{
    return mUpdateCellPopulation;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractCellBasedSimulationInterfaceFocus<ELEMENT_DIM,SPACE_DIM>::SetNoBirth(bool noBirth)
{
    mNoBirth = noBirth;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractCellBasedSimulationInterfaceFocus<ELEMENT_DIM,SPACE_DIM>::AddCellKiller(boost::shared_ptr<AbstractCellKiller<SPACE_DIM> > pCellKiller)
{
    mCellKillers.push_back(pCellKiller);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractCellBasedSimulationInterfaceFocus<ELEMENT_DIM,SPACE_DIM>::RemoveAllCellKillers()
{
    mCellKillers.clear();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::vector<double> AbstractCellBasedSimulationInterfaceFocus<ELEMENT_DIM,SPACE_DIM>::GetNodeLocation(const unsigned& rNodeIndex)
{
    std::vector<double> location;
    for (unsigned i=0; i<SPACE_DIM; i++)
    {
        location.push_back(mrCellPopulation.GetNode(rNodeIndex)->rGetLocation()[i]);
    }
    return location;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractCellBasedSimulationInterfaceFocus<ELEMENT_DIM,SPACE_DIM>::Solve()
{
    CellBasedEventHandler::BeginEvent(CellBasedEventHandler::EVERYTHING);
    CellBasedEventHandler::BeginEvent(CellBasedEventHandler::SETUP);

    // Set up the simulation time
    SimulationTime* p_simulation_time = SimulationTime::Instance();
    double current_time = p_simulation_time->GetTime();

    assert(mDt != DOUBLE_UNSET);  //Subclass constructors take care of this

    if (mEndTime == DOUBLE_UNSET)
    {
        EXCEPTION("SetEndTime has not yet been called.");
    }

    /*Note that mDt is used here for "ideal time step".  If this step doesn't divide the time remaining then
     * a *different* time step will be taken by the time-stepper.  The real time-step (used in the SimulationTime
     * singleton is currently not available to this class!
     * Should we over-write the value of mDt, or should we change this behaviour? \todo #2159
     */
    unsigned num_time_steps = (unsigned) ((mEndTime-current_time)/mDt+0.5);
    if (current_time > 0) // use the reset function if necessary
    {
        p_simulation_time->ResetEndTimeAndNumberOfTimeSteps(mEndTime, num_time_steps);
    }
    else
    {
        if (p_simulation_time->IsEndTimeAndNumberOfTimeStepsSetUp())
        {
            EXCEPTION("End time and number of timesteps already setup. You should not use SimulationTime::SetEndTimeAndNumberOfTimeSteps in cell-based tests.");
        }
        else
        {
            p_simulation_time->SetEndTimeAndNumberOfTimeSteps(mEndTime, num_time_steps);
        }
    }

    if (mOutputDirectory == "")
    {
        EXCEPTION("OutputDirectory not set");
    }

    double time_now = p_simulation_time->GetTime();
    std::ostringstream time_string;
    time_string << time_now;

    std::string results_directory = mOutputDirectory +"/results_from_time_" + time_string.str();
    mSimulationOutputDirectory = results_directory;

    // Set up simulation

    // Create output files for the visualizer
    OutputFileHandler output_file_handler(results_directory+"/", true);

    mrCellPopulation.CreateOutputFiles(results_directory+"/", false);

    if (mOutputDivisionLocations)
    {
        mpDivisionLocationFile = output_file_handler.OpenOutputFile("divisions.dat");
    }

    if (PetscTools::AmMaster())
    {
        mpVizSetupFile = output_file_handler.OpenOutputFile("results.vizsetup");
    }

    // If any PDEs have been defined, set up results files to store their solution
    if (mpCellBasedPdeHandler != NULL)
    {
        mpCellBasedPdeHandler->OpenResultsFiles(this->mSimulationOutputDirectory);
        if (PetscTools::AmMaster())
        {
            *this->mpVizSetupFile << "PDE \n";
        }

        // If any PDEs have been defined, solve them here before updating cells and store their solution in results files.
        // This also initializes the relevant CellData. NOTE that this works as the PDEs are elliptic
        CellBasedEventHandler::BeginEvent(CellBasedEventHandler::PDE);
        mpCellBasedPdeHandler->SolvePdeAndWriteResultsToFile(this->mSamplingTimestepMultiple);
        CellBasedEventHandler::EndEvent(CellBasedEventHandler::PDE);
    }
    SetupSolve();

    // Age the cells to the correct time. Note that cells are created with
    // negative birth times so that some are initially almost ready to divide.
    LOG(1, "Setting up cells...");
    for (typename AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>::Iterator cell_iter = mrCellPopulation.Begin();
         cell_iter != mrCellPopulation.End();
         ++cell_iter)
    {
        // We don't use the result; this call is just to force the cells to age
        // to the current time running their cell-cycle models to get there
        cell_iter->ReadyToDivide();
    }
    LOG(1, "\tdone\n");

    // Write initial conditions to file for the visualizer
    WriteVisualizerSetupFile();

    if (PetscTools::AmMaster())
    {
        *mpVizSetupFile << std::flush;
    }

    mrCellPopulation.WriteResultsToFiles();

    OutputSimulationSetup();

    CellBasedEventHandler::EndEvent(CellBasedEventHandler::SETUP);

    // Enter main time loop
    while (!( p_simulation_time->IsFinished() || StoppingEventHasOccurred() ) )
    {
        LOG(1, "--TIME = " << p_simulation_time->GetTime() << "\n");

        // This function calls DoCellRemoval(), DoCellBirth() and CellPopulation::Update()
        UpdateCellPopulation();

        // Update cell locations and topology
        UpdateCellLocationsAndTopology();

        // Increment simulation time here, so results files look sensible
        p_simulation_time->IncrementTimeOneStep();

        // If any PDEs have been defined, solve them and store their solution in results files
        if (mpCellBasedPdeHandler != NULL)
        {
            CellBasedEventHandler::BeginEvent(CellBasedEventHandler::PDE);
            mpCellBasedPdeHandler->SolvePdeAndWriteResultsToFile(this->mSamplingTimestepMultiple);
            CellBasedEventHandler::EndEvent(CellBasedEventHandler::PDE);
        }

        // Call UpdateAtEndOfTimeStep(), which may be implemented by child classes
        CellBasedEventHandler::BeginEvent(CellBasedEventHandler::UPDATESIMULATION);
        UpdateAtEndOfTimeStep();
        CellBasedEventHandler::EndEvent(CellBasedEventHandler::UPDATESIMULATION);

        // Output current results to file
        CellBasedEventHandler::BeginEvent(CellBasedEventHandler::OUTPUT);
        if (p_simulation_time->GetTimeStepsElapsed()%mSamplingTimestepMultiple == 0)
        {
            mrCellPopulation.WriteResultsToFiles();
        }
        CellBasedEventHandler::EndEvent(CellBasedEventHandler::OUTPUT);
    }

    LOG(1, "--END TIME = " << p_simulation_time->GetTime() << "\n");

    // Carry out a final update so that cell population is coherent with new cell positions.
    // NB cell birth/death still need to be checked because they may be spatially-dependent.
    UpdateCellPopulation();

    // If any PDEs have been defined, close the results files storing their solution
    if (mpCellBasedPdeHandler != NULL)
    {
        mpCellBasedPdeHandler->CloseResultsFiles();
    }

    CellBasedEventHandler::BeginEvent(CellBasedEventHandler::UPDATESIMULATION);
    UpdateAtEndOfSolve();
    CellBasedEventHandler::EndEvent(CellBasedEventHandler::UPDATESIMULATION);

    CellBasedEventHandler::BeginEvent(CellBasedEventHandler::OUTPUT);
    mrCellPopulation.CloseOutputFiles();

    if (mOutputDivisionLocations)
    {
        mpDivisionLocationFile->close();
    }

    if (PetscTools::AmMaster())
    {
        *mpVizSetupFile << "Complete\n";
        mpVizSetupFile->close();
    }

    CellBasedEventHandler::EndEvent(CellBasedEventHandler::OUTPUT);
    CellBasedEventHandler::EndEvent(CellBasedEventHandler::EVERYTHING);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool AbstractCellBasedSimulationInterfaceFocus<ELEMENT_DIM,SPACE_DIM>::StoppingEventHasOccurred()
{
    return false;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractCellBasedSimulationInterfaceFocus<ELEMENT_DIM,SPACE_DIM>::UpdateCellPopulation()
{
    // Remove dead cells
    CellBasedEventHandler::BeginEvent(CellBasedEventHandler::DEATH);
    unsigned deaths_this_step = DoCellRemoval();
    mNumDeaths += deaths_this_step;
    LOG(1, "\tNum deaths = " << mNumDeaths << "\n");
    CellBasedEventHandler::EndEvent(CellBasedEventHandler::DEATH);

    // Divide cells
    CellBasedEventHandler::BeginEvent(CellBasedEventHandler::BIRTH);
    unsigned births_this_step = DoCellBirth();
    mNumBirths += births_this_step;
    LOG(1, "\tNum births = " << mNumBirths << "\n");
    CellBasedEventHandler::EndEvent(CellBasedEventHandler::BIRTH);

    // This allows NodeBasedCellPopulation::Update() to do the minimum amount of work
    bool births_or_death_occurred = ((births_this_step>0) || (deaths_this_step>0));

    // Update topology of cell population
    CellBasedEventHandler::BeginEvent(CellBasedEventHandler::UPDATECELLPOPULATION);
    if (mUpdateCellPopulation)
    {
        LOG(1, "\tUpdating cell population...");
        mrCellPopulation.Update(births_or_death_occurred);
        LOG(1, "\tdone.\n");
    }
    else if (births_or_death_occurred)
    {
        EXCEPTION("CellPopulation has had births or deaths but mUpdateCellPopulation is set to false, please set it to true.");
    }
    CellBasedEventHandler::EndEvent(CellBasedEventHandler::UPDATECELLPOPULATION);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool AbstractCellBasedSimulationInterfaceFocus<ELEMENT_DIM,SPACE_DIM>::GetOutputDivisionLocations()
{
    return mOutputDivisionLocations;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractCellBasedSimulationInterfaceFocus<ELEMENT_DIM,SPACE_DIM>::SetOutputDivisionLocations(bool outputDivisionLocations)
{
    mOutputDivisionLocations = outputDivisionLocations;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractCellBasedSimulationInterfaceFocus<ELEMENT_DIM,SPACE_DIM>::OutputSimulationSetup()
{
    OutputFileHandler output_file_handler(this->mSimulationOutputDirectory + "/", false);

    // Output machine information
    ExecutableSupport::SetOutputDirectory(output_file_handler.GetOutputDirectoryFullPath());
    ExecutableSupport::WriteMachineInfoFile("system_info");
    // Output Chaste provenance information
    out_stream build_info_file = output_file_handler.OpenOutputFile("build.info");
    std::string build_info;
    ExecutableSupport::GetBuildInfo(build_info);
    *build_info_file << build_info;
    build_info_file->close();
    // Output simulation parameter and setup details
    out_stream parameter_file = output_file_handler.OpenOutputFile("results.parameters");

    // Output simulation details
    std::string simulation_type = "Test";

    *parameter_file << "<Chaste>\n";
    *parameter_file << "\n\t<" << simulation_type << ">\n";
    OutputSimulationParameters(parameter_file);
    *parameter_file << "\t</" << simulation_type << ">\n";
    *parameter_file << "\n";

    // Output cell population details (includes cell-cycle model details)
    mrCellPopulation.OutputCellPopulationInfo(parameter_file);
    // Loop over cell killers
    *parameter_file << "\n\t<CellKillers>\n";
    for (typename std::vector<boost::shared_ptr<AbstractCellKiller<SPACE_DIM> > >::iterator iter = mCellKillers.begin();
         iter != mCellKillers.end();
         ++iter)
    {
        // Output cell killer details
        (*iter)->OutputCellKillerInfo(parameter_file);
    }
    *parameter_file << "\t</CellKillers>\n";
    // This is used to output information about subclasses
    OutputAdditionalSimulationSetup(parameter_file);
    *parameter_file << "\n</Chaste>\n";
    parameter_file->close();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractCellBasedSimulationInterfaceFocus<ELEMENT_DIM,SPACE_DIM>::OutputSimulationParameters(out_stream& rParamsFile)
{

	if (mpCellBasedPdeHandler != NULL)
    {
        //mpCellBasedPdeHandler->OutputParameters(rParamsFile);
    }

    *rParamsFile << "\t\t<Dt>" << mDt << "</Dt>\n";
    *rParamsFile << "\t\t<EndTime>" << mEndTime << "</EndTime>\n";
    *rParamsFile << "\t\t<SamplingTimestepMultiple>" << mSamplingTimestepMultiple << "</SamplingTimestepMultiple>\n";
    *rParamsFile << "\t\t<OutputDivisionLocations>" << mOutputDivisionLocations << "</OutputDivisionLocations>\n";
}

////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
////////////////////////////////////////////////////////////////////////////

template class AbstractCellBasedSimulationInterfaceFocus<1,1>;
template class AbstractCellBasedSimulationInterfaceFocus<1,2>;
template class AbstractCellBasedSimulationInterfaceFocus<2,2>;
template class AbstractCellBasedSimulationInterfaceFocus<1,3>;
template class AbstractCellBasedSimulationInterfaceFocus<2,3>;
template class AbstractCellBasedSimulationInterfaceFocus<3,3>;
