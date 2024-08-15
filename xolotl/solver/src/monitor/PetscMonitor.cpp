#include <xolotl/io/XFile.h>
#include <xolotl/perf/ScopedTimer.h>
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
	PetscFunctionBeginUser;

	// Get the time step from ts
	PetscReal timestep;
	PetscCall(TSGetTimeStep(ts, &timestep));

	// Stop when the time step is lower than the user defined threshold
	if (timestep < timeStepThreshold) {
		PetscCall(TSSetConvergedReason(ts, TS_CONVERGED_EVENT));
	}

	PetscFunctionReturn(0);
}

PetscErrorCode
monitorTime(TS ts, PetscInt timestep, PetscReal time, Vec solution, void* ictx)
{
	PetscFunctionBeginUser;
	PetscCall(static_cast<IPetscMonitor*>(ictx)->monitorTime(
		ts, timestep, time, solution));
	PetscFunctionReturn(0);
}

PetscErrorCode
computeFluence(
	TS ts, PetscInt timestep, PetscReal time, Vec solution, void* ictx)
{
	PetscFunctionBeginUser;
	PetscCall(static_cast<IPetscMonitor*>(ictx)->computeFluence(
		ts, timestep, time, solution));
	PetscFunctionReturn(0);
}

PetscErrorCode
monitorPerf(TS ts, PetscInt timestep, PetscReal time, Vec solution, void* ictx)
{
	PetscFunctionBeginUser;
	PetscCall(static_cast<IPetscMonitor*>(ictx)->monitorPerf(
		ts, timestep, time, solution));
	PetscFunctionReturn(0);
}

PetscErrorCode
monitorLargest(
	TS ts, PetscInt timestep, PetscReal time, Vec solution, void* ictx)
{
	PetscFunctionBeginUser;
	PetscCall(static_cast<IPetscMonitor*>(ictx)->monitorLargest(
		ts, timestep, time, solution));
	PetscFunctionReturn(0);
}

PetscErrorCode
startStop(TS ts, PetscInt timestep, PetscReal time, Vec solution, void* ictx)
{
	PetscFunctionBeginUser;
	PetscCall(static_cast<IPetscMonitor*>(ictx)->startStop(
		ts, timestep, time, solution));
	PetscFunctionReturn(0);
}

PetscErrorCode
computeHeliumRetention(
	TS ts, PetscInt timestep, PetscReal time, Vec solution, void* ictx)
{
	PetscFunctionBeginUser;
	PetscCall(static_cast<IPetscMonitor*>(ictx)->computeHeliumRetention(
		ts, timestep, time, solution));
	PetscFunctionReturn(0);
}

PetscErrorCode
computeXenonRetention(
	TS ts, PetscInt timestep, PetscReal time, Vec solution, void* ictx)
{
	PetscFunctionBeginUser;
	PetscCall(static_cast<IPetscMonitor*>(ictx)->computeXenonRetention(
		ts, timestep, time, solution));
	PetscFunctionReturn(0);
}

/**
 * This is a monitoring method that will compute average density and diameter
 * of defects.
 */
PetscErrorCode
computeAlphaZr(
	TS ts, PetscInt timestep, PetscReal time, Vec solution, void* ictx)
{
	PetscFunctionBeginUser;
	PetscCall(static_cast<IPetscMonitor*>(ictx)->computeAlphaZr(
		ts, timestep, time, solution));
	PetscFunctionReturn(0);
}

PetscErrorCode
computeAlloy(TS ts, PetscInt timestep, PetscReal time, Vec solution, void* ictx)
{
	PetscFunctionBeginUser;
	PetscCall(static_cast<IPetscMonitor*>(ictx)->computeAlloy(
		ts, timestep, time, solution));
	PetscFunctionReturn(0);
}

PetscErrorCode
monitorScatter(
	TS ts, PetscInt timestep, PetscReal time, Vec solution, void* ictx)
{
	PetscFunctionBeginUser;
	PetscCall(static_cast<IPetscMonitor*>(ictx)->monitorScatter(
		ts, timestep, time, solution));
	PetscFunctionReturn(0);
}

PetscErrorCode
eventFunction(
	TS ts, PetscReal time, Vec solution, PetscScalar* fvalue, void* ctx)
{
	PetscFunctionBeginUser;
	PetscCall(static_cast<IPetscMonitor*>(ctx)->eventFunction(
		ts, time, solution, fvalue));
	PetscFunctionReturn(0);
}

PetscErrorCode
postEventFunction(TS ts, PetscInt nevents, PetscInt eventList[], PetscReal time,
	Vec solution, PetscBool fwd, void* ctx)
{
	PetscFunctionBeginUser;
	PetscCall(static_cast<IPetscMonitor*>(ctx)->postEventFunction(
		ts, nevents, eventList, time, solution, fwd));
	PetscFunctionReturn(0);
}

PetscMonitor::PetscMonitor() = default;

PetscMonitor::PetscMonitor(TS ts,
	const std::shared_ptr<handler::ISolverHandler>& solverHandler,
	const std::string& checkpointFileName) :
	_ts(ts),
	_solverHandler(solverHandler),
	_hdf5OutputName(checkpointFileName)
{
	// Create the checkpoint file if necessary
	if (this->checkForCreatingCheckpoint()) {
		io::XFile checkpointFile(_hdf5OutputName, 1, util::getMPIComm());
	}
}

PetscMonitor::~PetscMonitor()
{
}

bool
PetscMonitor::checkForCreatingCheckpoint() const
{
	return (not _hdf5OutputName.empty()) and (not fs::exists(_hdf5OutputName));
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
PetscMonitor::startStop(TS ts, PetscInt timestep, PetscReal time, Vec solution)
{
	PetscFunctionBeginUser;
	perf::ScopedTimer myTimer(_startStopTimer);

	// Compute the dt
	double previousTime = _solverHandler->getPreviousTime();
	double dt = time - previousTime;

	// Don't do anything if it is not on the stride
	if (((PetscInt)((time + dt / 10.0) / _hdf5Stride) <= _hdf5Previous) &&
		(time > 0.0)) {
		PetscFunctionReturn(0);
	}

	// Update the previous time
	if ((PetscInt)((time + dt / 10.0) / _hdf5Stride) > _hdf5Previous) {
		_hdf5Previous++;
	}

	// Gets MPI comm
	auto xolotlComm = util::getMPIComm();

	// Open the existing HDF5 file
	io::XFile checkpointFile(
		_hdf5OutputName, xolotlComm, io::XFile::AccessMode::OpenReadWrite);

	// Get the current time step
	double currentTimeStep;
	PetscCall(TSGetTimeStep(ts, &currentTimeStep));

	// Add a concentration time step group for the current time step.
	auto concGroup = checkpointFile.getGroup<io::XFile::ConcentrationGroup>();
	assert(concGroup);
	auto tsGroup = concGroup->addTimestepGroup(
		_ctrlStep, _loopNumber, timestep, time, previousTime, currentTimeStep);

	// Save the fluence
	auto fluxHandler = _solverHandler->getFluxHandler();
	auto fluence = fluxHandler->getFluence();
	tsGroup->writeFluence(fluence);

	// Get the names of the species in the network
	auto& network = _solverHandler->getNetwork();
	auto numSpecies = network.getSpeciesListSize();
	std::vector<std::string> speciesNames;
	for (auto id = core::network::SpeciesId(numSpecies); id; ++id) {
		speciesNames.push_back(network.getSpeciesName(id));
	}

	this->startStopImpl(ts, timestep, time, solution, checkpointFile,
		tsGroup.get(), speciesNames);

	PetscFunctionReturn(0);
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

	// If it's a restart
	bool hasConcentrations = false;
	if (_solverHandler->checkForRestart()) {
		auto restartFilePath = _solverHandler->getRestartFilePath();
		auto networkFile = std::make_unique<io::XFile>(restartFilePath);
		auto concGroup = networkFile->getGroup<io::XFile::ConcentrationGroup>();
		hasConcentrations = (concGroup and concGroup->hasTimesteps());
	}

	if (timestep == 0 and hasConcentrations)
		PetscFunctionReturn(0);

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
	auto solverTimer = _solverHandler->getPerfHandler()->getTimer("solve");
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
		PetscCall(TSGetTimeStep(ts, &currentTimeStep));
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
PetscMonitor::computeAlphaZr(
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
