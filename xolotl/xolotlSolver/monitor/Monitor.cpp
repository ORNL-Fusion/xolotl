// Includes
#include "PetscSolver.h"
#include <xolotlPerf/xolotlPerf.h>
#include <VizHandlerRegistryFactory.h>
#include <CvsXDataProvider.h>
#include <LabelProvider.h>
#include <petscts.h>
#include <petscsys.h>
#include <sstream>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <memory>
#include "xolotlCore/io/XFile.h"
#include "xolotlSolver/monitor/Monitor.h"

namespace xolotlSolver {

//! The pointer to the plot that will be used to visualize performance data.
std::shared_ptr<xolotlViz::IPlot> perfPlot;

//! The variable to store the time at the previous time step.
double previousTime = 0.0;
//! The variable to store the threshold on time step defined by the user.
double timeStepThreshold = 0.0;

#undef __FUNCT__
#define __FUNCT__ Actual__FUNCT__("xolotlSolver", "checkTimeStep")
/**
 * This is a method that decides when to extend the network
 */
PetscErrorCode checkTimeStep(TS ts) {
	// Initial declarations
	PetscErrorCode ierr;

	PetscFunctionBeginUser;

	// Get the time step from ts
	PetscReal timestep;
	ierr = TSGetTimeStep(ts, &timestep);
	CHKERRQ(ierr);

	// Stop when the time step is lower than the user defined threshold
	if (timestep < timeStepThreshold) {
		ierr = TSSetConvergedReason(ts, TS_CONVERGED_EVENT);
		CHKERRQ(ierr);
	}

	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ Actual__FUNCT__("xolotlSolver", "monitorTime")
/**
 * This is a monitoring method set the previous time to the time. This is needed here
 * because multiple monitors need the previous time value from the previous timestep.
 * This monitor must be called last when needed.
 */
PetscErrorCode monitorTime(TS, PetscInt, PetscReal time, Vec, void *) {
	PetscFunctionBeginUser;

	// Set the previous time to the current time for the next timestep
	previousTime = time;

	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ Actual__FUNCT__("xolotlSolver", "computeHeliumFluence")
/**
 * This is a monitoring method that will compute the total helium fluence
 */
PetscErrorCode computeFluence(TS, PetscInt, PetscReal time, Vec, void *) {
	PetscFunctionBeginUser;

	// Get the solver handler and the flux handler
	auto& solverHandler = PetscSolver::getSolverHandler();
	auto fluxHandler = solverHandler.getFluxHandler();

	// The length of the time step
	double dt = time - previousTime;

	// Increment the fluence with the value at this current timestep
	fluxHandler->incrementFluence(dt);

	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ Actual__FUNCT__("xolotlSolver", "monitorPerf")
/**
 * This is a monitoring method that will save 1D plots of one performance timer
 */
PetscErrorCode monitorPerf(TS ts, PetscInt timestep, PetscReal time, Vec,
		void *) {
	// To check PETSc errors
	PetscInt ierr;

	PetscFunctionBeginUser;

	// Get the number of processes
	int cwSize;
	int cwRank;
	MPI_Comm_size(PETSC_COMM_WORLD, &cwSize);
	MPI_Comm_rank(PETSC_COMM_WORLD, &cwRank);

	// Print a warning if only one process
	if (cwSize == 1) {
		std::cout
				<< "You are trying to plot things that don't have any sense!! "
				<< "\nRemove -plot_perf or run in parallel." << std::endl;
		PetscFunctionReturn(0);
	}

	// Obtain the current value of the solve timer.
	//
	// Note that the solve timer keeps a cumulative time,
	// not a per-timestep time.   If you need a per-timestep
	// time, you will want to keep a static or global variable
	// with the last known timer value, and subtract it from
	// the current timer value each time this monitor function is called.
	//
	// Note also that we restart the timer immediately after sampling
	// its value.  If you feel it is "unfair" to charge the time
	// required for the rank 0 process to produce the output plot
	// against the solve timer, then you should move the start()
	// call after the code that produces the plot (and probably also
	// put in an MPI_Barrier before starting the timer so that
	// all processes avoid including the time required for rank 0
	// to produce the plot).  We probably don't want to reset the
	// timer here since the main function is using it to get an
	// overall elapsed time measurement of the solve.
	//
	auto solverTimer = xolotlPerf::getHandlerRegistry()->getTimer("solve");
	solverTimer->stop();
	double solverTimerValue = solverTimer->getValue();
	solverTimer->start();

	// Collect all sampled timer values to rank 0.
	double* allTimerValues = (cwRank == 0) ? new double[cwSize] : NULL;
	MPI_Gather(&solverTimerValue,  // send buffer
			1,// number of values to send
			MPI_DOUBLE,// type of items in send buffer
			allTimerValues,// receive buffer (only valid at root)
			1,// number of values to receive from each process
			MPI_DOUBLE,// type of items in receive buffer
			0,// root of MPI collective operation
			PETSC_COMM_WORLD);// communicator defining processes involved in the operation

	if (cwRank == 0) {
		auto allPoints = std::make_shared<std::vector<xolotlViz::Point> >();

		for (int i = 0; i < cwSize; ++i) {
			xolotlViz::Point aPoint;
			aPoint.value = allTimerValues[i];
			aPoint.x = i;
			aPoint.t = time;
			allPoints->push_back(aPoint);
		}

		// Clean up the receive buffer (only valid at root)
		delete[] allTimerValues;

		// Provide the data provider the points.
		perfPlot->getDataProvider()->setPoints(allPoints);
		perfPlot->getDataProvider()->setDataName("SolverTimer");

		// Change the title of the plot
		std::ostringstream title;
		title << "Solver timer (s)";
		perfPlot->plotLabelProvider->titleLabel = title.str();
		// Give the time to the label provider
		std::ostringstream timeLabel;
		timeLabel << "time: " << std::setprecision(4) << time << "s";
		perfPlot->plotLabelProvider->timeLabel = timeLabel.str();
		// Get the current time step
		PetscReal currentTimeStep;
		ierr = TSGetTimeStep(ts, &currentTimeStep);
		CHKERRQ(ierr);
		// Give the timestep to the label provider
		std::ostringstream timeStepLabel;
		timeStepLabel << "dt: " << std::setprecision(4) << currentTimeStep
				<< "s";
		perfPlot->plotLabelProvider->timeStepLabel = timeStepLabel.str();

		// Render and save in file
		std::ostringstream fileName;
		fileName << "timer_TS" << timestep << ".pnm";
		perfPlot->write(fileName.str());
	}

	PetscFunctionReturn(0);
}

void writeNetwork(MPI_Comm _comm, std::string srcFileName,
		std::string targetFileName, experimental::IReactionNetwork& network) {

	int procId;
	MPI_Comm_rank(_comm, &procId);

	// Check if we are supposed to copy the network from
	// another object into our new checkpoint file.
	if (procId == 0) {
		if (not srcFileName.empty()) {

			// Copy the network from the given file.
			// Note that we do this using a single-process
			// communicator because the HDF5 copy operation
			// is not parallelized and gives very poor performance
			// if used with a file opened for parallel access.
			xolotlCore::XFile srcFile(srcFileName,
			MPI_COMM_SELF, xolotlCore::XFile::AccessMode::OpenReadOnly);

			// Check if given file even has a network group.
			auto srcNetGroup =
					srcFile.getGroup<xolotlCore::XFile::NetworkGroup>();
			if (srcNetGroup) {
				// Given file has a network group.  Copy it.
				// First open the checkpoint file using a single-process
				// communicator...
				xolotlCore::XFile checkpointFile(targetFileName,
				MPI_COMM_SELF, xolotlCore::XFile::AccessMode::OpenReadWrite);

				// ...then do the copy.
				srcNetGroup->copyTo(checkpointFile);
			}
		} else {
			// Write from scratch
			xolotlCore::XFile checkpointFile(targetFileName,
			MPI_COMM_SELF, xolotlCore::XFile::AccessMode::OpenReadWrite);
			xolotlCore::XFile::NetworkGroup netGroup(checkpointFile, network);
		}
	}
}

}
/* end namespace xolotlSolver */
