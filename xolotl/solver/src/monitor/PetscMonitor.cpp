#include <xolotl/io/XFile.h>
#include <xolotl/perf/PerfHandlerRegistry.h>
#include <xolotl/solver/monitor/PetscMonitor.h>
#include <xolotl/util/MPIUtils.h>

namespace xolotl
{
namespace solver
{
namespace monitor
{
//! The variable to store the threshold on time step defined by the user.
double timeStepThreshold = 0.0;

/**
 * This is a method that decides when to extend the network
 */
PetscErrorCode
checkTimeStep(TS ts)
{
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

PetscErrorCode
monitorTime(TS ts, PetscInt timestep, PetscReal time, Vec solution, void* ictx)
{
	PetscFunctionBeginUser;
	PetscErrorCode ierr = static_cast<IPetscMonitor*>(ictx)->monitorTime(
		ts, timestep, time, solution);
	CHKERRQ(ierr);
	PetscFunctionReturn(0);
}

PetscErrorCode
computeFluence(
	TS ts, PetscInt timestep, PetscReal time, Vec solution, void* ictx)
{
	PetscFunctionBeginUser;
	PetscErrorCode ierr = static_cast<IPetscMonitor*>(ictx)->computeFluence(
		ts, timestep, time, solution);
	CHKERRQ(ierr);
	PetscFunctionReturn(0);
}

PetscErrorCode
monitorPerf(TS ts, PetscInt timestep, PetscReal time, Vec solution, void* ictx)
{
	PetscFunctionBeginUser;
	PetscErrorCode ierr = static_cast<IPetscMonitor*>(ictx)->monitorPerf(
		ts, timestep, time, solution);
	CHKERRQ(ierr);
	PetscFunctionReturn(0);
}

PetscErrorCode
monitorLargest(
	TS ts, PetscInt timestep, PetscReal time, Vec solution, void* ictx)
{
	PetscFunctionBeginUser;
	PetscErrorCode ierr = static_cast<IPetscMonitor*>(ictx)->monitorLargest(
		ts, timestep, time, solution);
	CHKERRQ(ierr);
	PetscFunctionReturn(0);
}

PetscErrorCode
startStop(TS ts, PetscInt timestep, PetscReal time, Vec solution, void* ictx)
{
	PetscFunctionBeginUser;
	PetscErrorCode ierr = static_cast<IPetscMonitor*>(ictx)->startStop(
		ts, timestep, time, solution);
	CHKERRQ(ierr);
	PetscFunctionReturn(0);
}

PetscErrorCode
computeHeliumRetention(
	TS ts, PetscInt timestep, PetscReal time, Vec solution, void* ictx)
{
	PetscFunctionBeginUser;
	PetscErrorCode ierr =
		static_cast<IPetscMonitor*>(ictx)->computeHeliumRetention(
			ts, timestep, time, solution);
	CHKERRQ(ierr);
	PetscFunctionReturn(0);
}

PetscErrorCode
computeXenonRetention(
	TS ts, PetscInt timestep, PetscReal time, Vec solution, void* ictx)
{
	PetscFunctionBeginUser;
	PetscErrorCode ierr =
		static_cast<IPetscMonitor*>(ictx)->computeXenonRetention(
			ts, timestep, time, solution);
	CHKERRQ(ierr);
	PetscFunctionReturn(0);
}

PetscErrorCode
computeAlloy(TS ts, PetscInt timestep, PetscReal time, Vec solution, void* ictx)
{
	PetscFunctionBeginUser;
	PetscErrorCode ierr = static_cast<IPetscMonitor*>(ictx)->computeAlloy(
		ts, timestep, time, solution);
	CHKERRQ(ierr);
	PetscFunctionReturn(0);
}

PetscErrorCode
monitorScatter(
	TS ts, PetscInt timestep, PetscReal time, Vec solution, void* ictx)
{
	PetscFunctionBeginUser;
	PetscErrorCode ierr = static_cast<IPetscMonitor*>(ictx)->monitorScatter(
		ts, timestep, time, solution);
	CHKERRQ(ierr);
	PetscFunctionReturn(0);
}

PetscErrorCode
eventFunction(
	TS ts, PetscReal time, Vec solution, PetscScalar* fvalue, void* ctx)
{
	PetscFunctionBeginUser;
	PetscErrorCode ierr = static_cast<IPetscMonitor*>(ctx)->eventFunction(
		ts, time, solution, fvalue);
	CHKERRQ(ierr);
	PetscFunctionReturn(0);
}

PetscErrorCode
postEventFunction(TS ts, PetscInt nevents, PetscInt eventList[], PetscReal time,
	Vec solution, PetscBool fwd, void* ctx)
{
	PetscFunctionBeginUser;
	PetscErrorCode ierr = static_cast<IPetscMonitor*>(ctx)->postEventFunction(
		ts, nevents, eventList, time, solution, fwd);
	CHKERRQ(ierr);
	PetscFunctionReturn(0);
}

PetscMonitor::PetscMonitor() = default;

PetscMonitor::PetscMonitor(
	TS ts, const std::shared_ptr<handler::ISolverHandler>& solverHandler) :
	_ts(ts),
	_solverHandler(solverHandler)
{
}

PetscMonitor::~PetscMonitor()
{
}

void
PetscMonitor::writeNetwork(MPI_Comm comm, const std::string& targetFileName,
	const std::string& srcFileName)
{
	int procId;
	MPI_Comm_rank(comm, &procId);

	// Check if we are supposed to copy the network from
	// another object into our new checkpoint file.
	if (procId == 0) {
		if (not srcFileName.empty()) {
			// Copy the network from the given file.
			// Note that we do this using a single-process
			// communicator because the HDF5 copy operation
			// is not parallelized and gives very poor performance
			// if used with a file opened for parallel access.
			io::XFile srcFile(srcFileName, MPI_COMM_SELF,
				io::XFile::AccessMode::OpenReadOnly);

			// Check if given file even has a network group.
			auto srcNetGroup = srcFile.getGroup<io::XFile::NetworkGroup>();
			if (srcNetGroup) {
				// Given file has a network group.  Copy it.
				// First open the checkpoint file using a single-process
				// communicator...
				io::XFile checkpointFile(targetFileName, MPI_COMM_SELF,
					io::XFile::AccessMode::OpenReadWrite);

				// ...then do the copy.
				srcNetGroup->copyTo(checkpointFile);
			}
		}
		else {
			// Write from scratch
			io::XFile checkpointFile(targetFileName, MPI_COMM_SELF,
				io::XFile::AccessMode::OpenReadWrite);
			io::XFile::NetworkGroup netGroup(
				checkpointFile, _solverHandler->getNetwork());
		}
	}
}

PetscErrorCode
PetscMonitor::monitorTime(
	TS ts, PetscInt timestep, PetscReal time, Vec solution)
{
	PetscFunctionBeginUser;

	// Set the previous time to the current time for the next timestep
	_solverHandler->setPreviousTime(time);

	PetscFunctionReturn(0);
}

PetscErrorCode
PetscMonitor::computeFluence(
	TS ts, PetscInt timestep, PetscReal time, Vec solution)
{
	PetscFunctionBeginUser;

	// Get the flux handler
	auto fluxHandler = _solverHandler->getFluxHandler();

	// The length of the time step
	double dt = time - _solverHandler->getPreviousTime();

	// Increment the fluence with the value at this current timestep
	fluxHandler->incrementFluence(dt);

	PetscFunctionReturn(0);
}

PetscErrorCode
PetscMonitor::monitorPerf(
	TS ts, PetscInt timestep, PetscReal time, Vec solution)
{
	// To check PETSc errors
	PetscErrorCode ierr;

	PetscFunctionBeginUser;

	// Get the number of processes
	auto xolotlComm = util::getMPIComm();
	int cwSize;
	int cwRank;
	MPI_Comm_size(xolotlComm, &cwSize);
	MPI_Comm_rank(xolotlComm, &cwRank);

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
	auto solverTimer = perf::PerfHandlerRegistry::get()->getTimer("solve");
	solverTimer->stop();
	double solverTimerValue = solverTimer->getValue();
	solverTimer->start();

	// Collect all sampled timer values to rank 0.
	std::vector<double> allTimerValues((cwRank == 0) ? cwSize : 0);
	MPI_Gather(&solverTimerValue, // send buffer
		1, // number of values to send
		MPI_DOUBLE, // type of items in send buffer
		allTimerValues.data(), // receive buffer (only valid at root)
		1, // number of values to receive from each process
		MPI_DOUBLE, // type of items in receive buffer
		0, // root of MPI collective operation
		xolotlComm); // communicator defining processes involved in the
					 // operation

	if (cwRank == 0) {
		auto allPoints =
			std::make_shared<std::vector<viz::dataprovider::DataPoint>>();

		for (auto i = 0; i < cwSize; ++i) {
			viz::dataprovider::DataPoint aPoint;
			aPoint.value = allTimerValues[i];
			aPoint.x = i;
			aPoint.t = time;
			allPoints->push_back(aPoint);
		}

		// Provide the data provider the points.
		_perfPlot->getDataProvider()->setDataPoints(allPoints);
		_perfPlot->getDataProvider()->setDataName("SolverTimer");

		// Change the title of the plot
		std::ostringstream title;
		title << "Solver timer (s)";
		_perfPlot->plotLabelProvider->titleLabel = title.str();
		// Give the time to the label provider
		std::ostringstream timeLabel;
		timeLabel << "time: " << std::setprecision(4) << time << "s";
		_perfPlot->plotLabelProvider->timeLabel = timeLabel.str();
		// Get the current time step
		PetscReal currentTimeStep;
		ierr = TSGetTimeStep(ts, &currentTimeStep);
		CHKERRQ(ierr);
		// Give the timestep to the label provider
		std::ostringstream timeStepLabel;
		timeStepLabel << "dt: " << std::setprecision(4) << currentTimeStep
					  << "s";
		_perfPlot->plotLabelProvider->timeStepLabel = timeStepLabel.str();

		// Render and save in file
		std::ostringstream fileName;
		fileName << "timer_TS" << timestep << ".pnm";
		_perfPlot->write(fileName.str());
	}

	PetscFunctionReturn(0);
}

PetscErrorCode
PetscMonitor::computeHeliumRetention(
	TS ts, PetscInt timestep, PetscReal time, Vec solution)
{
	PetscFunctionReturn(0);
}

PetscErrorCode
PetscMonitor::computeAlloy(
	TS ts, PetscInt timestep, PetscReal time, Vec solution)
{
	PetscFunctionReturn(0);
}

PetscErrorCode
PetscMonitor::monitorScatter(
	TS ts, PetscInt timestep, PetscReal time, Vec solution)
{
	PetscFunctionReturn(0);
}

PetscErrorCode
PetscMonitor::eventFunction(
	TS ts, PetscReal time, Vec solution, PetscScalar* fvalue)
{
	PetscFunctionReturn(0);
}

PetscErrorCode
PetscMonitor::postEventFunction(TS ts, PetscInt nevents, PetscInt eventList[],
	PetscReal time, Vec solution, PetscBool)
{
	PetscFunctionReturn(0);
}
} // namespace monitor
} // namespace solver
} // namespace xolotl
