// Includes
#include <ctime>
#include <iostream>

#include <xolotl/core/network/INetworkHandler.h>
#include <xolotl/factory/perf/PerfHandlerFactory.h>
#include <xolotl/factory/solver/SolverFactory.h>
#include <xolotl/factory/viz/VizHandlerFactory.h>
#include <xolotl/interface/Interface.h>
#include <xolotl/options/Options.h>
#include <xolotl/perf/PerfHandlerRegistry.h>
#include <xolotl/solver/Solver.h>
#include <xolotl/solver/handler/ISolverHandler.h>
#include <xolotl/util/MPIUtils.h>
#include <xolotl/version.h>
#include <xolotl/viz/VizHandlerRegistry.h>

namespace xolotl
{
namespace interface
{
class Context
{
public:
	Context(int argc, const char* argv[]) :
		_kokkosContext(argc, const_cast<char**>(argv))
	{
		if (!initialized()) {
			util::mpiInit(argc, argv);
			_mpiInitializedHere = true;
		}
	}

	~Context()
	{
		if (_mpiInitializedHere) {
			if (!finalized()) {
				MPI_Finalize();
			}
		}
	}

	static bool
	initialized()
	{
		int flag;
		MPI_Initialized(&flag);
		return flag != 0;
	}

	static bool
	finalized()
	{
		int flag;
		MPI_Finalized(&flag);
		return flag != 0;
	}

private:
	bool _mpiInitializedHere{false};
	Kokkos::ScopeGuard _kokkosContext;
};

void
reportException(const std::exception& e)
{
	std::cerr << e.what() << "\nAborting." << std::endl;
}

XolotlInterface::XolotlInterface() = default;

XolotlInterface::XolotlInterface(int argc, const char* argv[], MPI_Comm mpiComm)
{
	initializeXolotl(argc, argv, mpiComm);
	initializedHere = true;
}

XolotlInterface::~XolotlInterface()
{
	if (initializedHere) {
		finalizeXolotl();
	}
}

void
XolotlInterface::printSomething()
{
	std::cout << "I'm in Xolotl !!!" << std::endl;
}

void
XolotlInterface::initializeXolotl(int argc, const char* argv[], MPI_Comm comm)
try {
	context = std::make_unique<Context>(argc, argv);

	// Initialize the MPI communicator to use
	util::setMPIComm(comm);
	auto xolotlComm = util::getMPIComm();

	// Get the MPI rank
	int rank;
	MPI_Comm_rank(xolotlComm, &rank);

	if (rank == 0) {
		// Print the start message
		std::cout << "Starting Xolotl (" << getExactVersionString() << ")\n";
		// TODO! Print copyright message
		// Print date and time
		std::time_t currentTime = std::time(NULL);
		std::cout << std::asctime(std::localtime(&currentTime)) << std::flush;
	}

	options::Options opts;
	opts.readParams(argc, argv);
	if (!opts.shouldRun()) {
		throw std::runtime_error("Unable to read the options.");
	}

	// Set up our performance data infrastructure.
	perf::PerfHandlerRegistry::set(
		factory::perf::PerfHandlerFactory::get().generate(opts));

	// Initialize the visualization
	viz::VizHandlerRegistry::set(
		factory::viz::VizHandlerFactory::get().generate(opts));

	// Setup the solver
	solver = factory::solver::SolverFactory::get().generate(opts);
	assert(solver);
	// Initialize the solver
	solver->initialize();
}
catch (const std::exception& e) {
	reportException(e);
	throw;
}

void
XolotlInterface::setTimes(double finalTime, double dt)
try {
	// Set the time in the solver
	solver->setTimes(finalTime, dt);
}
catch (const std::exception& e) {
	reportException(e);
	throw;
}

void
XolotlInterface::solveXolotl()
try {
	// Launch the PetscSolver
	solver->solve();
}
catch (const std::exception& e) {
	reportException(e);
	throw;
}

std::vector<std::vector<std::vector<std::array<double, 4>>>>
XolotlInterface::getLocalNE()
try {
	// Get the solver handler and return the rate vector
	return solver::Solver::getSolverHandler().getLocalNE();
}
catch (const std::exception& e) {
	reportException(e);
	throw;
}

void
XolotlInterface::setLocalNE(
	const std::vector<std::vector<std::vector<std::array<double, 4>>>>&
		rateVector)
try {
	// Get the solver handler
	auto& solverHandler = solver::Solver::getSolverHandler();
	// Set the rate vector
	solverHandler.setLocalNE(rateVector);
}
catch (const std::exception& e) {
	reportException(e);
	throw;
}

void
XolotlInterface::getLocalCoordinates(IdType& xs, IdType& xm, IdType& Mx,
	IdType& ys, IdType& ym, IdType& My, IdType& zs, IdType& zm, IdType& Mz)
try {
	// Get the solver handler
	auto& solverHandler = solver::Solver::getSolverHandler();
	// Get the local coordinates
	solverHandler.getLocalCoordinates(xs, xm, Mx, ys, ym, My, zs, zm, Mz);
}
catch (const std::exception& e) {
	reportException(e);
	throw;
}

void
XolotlInterface::setGBLocation(IdType i, IdType j, IdType k)
try {
	// Get the solver handler
	auto& solverHandler = solver::Solver::getSolverHandler();
	// Set the coordinate of the GB
	solverHandler.setGBLocation(i, j, k);
}
catch (const std::exception& e) {
	reportException(e);
	throw;
}

void
XolotlInterface::resetGBVector()
try {
	// Get the solver handler
	auto& solverHandler = solver::Solver::getSolverHandler();
	// Reset the location
	solverHandler.resetGBVector();
}
catch (const std::exception& e) {
	reportException(e);
	throw;
}

std::vector<std::vector<std::vector<std::vector<std::pair<IdType, double>>>>>
XolotlInterface::getConcVector()
try {
	// Get the vector
	return solver->getConcVector();
}
catch (const std::exception& e) {
	reportException(e);
	throw;
}

void
XolotlInterface::setConcVector(std::vector<
	std::vector<std::vector<std::vector<std::pair<IdType, double>>>>>
		concVector)
try {
	// Set the vector
	solver->setConcVector(concVector);
}
catch (const std::exception& e) {
	reportException(e);
	throw;
}

double
XolotlInterface::getPreviousTime()
try {
	// Get the solver handler
	auto& solverHandler = solver::Solver::getSolverHandler();
	return solverHandler.getPreviousTime();
}
catch (const std::exception& e) {
	reportException(e);
	throw;
}

void
XolotlInterface::setPreviousTime(double time)
try {
	// Get the solver handler
	auto& solverHandler = solver::Solver::getSolverHandler();
	solverHandler.setPreviousTime(time, true); // Update the fluence from here
}
catch (const std::exception& e) {
	reportException(e);
	throw;
}

double
XolotlInterface::getCurrentDt()
try {
	return solver->getCurrentDt();
}
catch (const std::exception& e) {
	reportException(e);
	throw;
}

void
XolotlInterface::setCurrentTimes(double time, double dt)
try {
	solver->setCurrentTimes(time, dt);
}
catch (const std::exception& e) {
	reportException(e);
	throw;
}

double
XolotlInterface::getNXeGB()
try {
	// Get the solver handler
	auto& solverHandler = solver::Solver::getSolverHandler();
	return solverHandler.getNXeGB();
}
catch (const std::exception& e) {
	reportException(e);
	throw;
}

void
XolotlInterface::setNXeGB(double nXe)
try {
	// Get the solver handler
	auto& solverHandler = solver::Solver::getSolverHandler();
	solverHandler.setNXeGB(nXe);
}
catch (const std::exception& e) {
	reportException(e);
	throw;
}

TS&
XolotlInterface::getTS()
{
	return solver->getTS();
}

std::vector<double>
XolotlInterface::getGridInfo(double& hy, double& hz)
try {
	// Get the solver handler
	auto& solverHandler = solver::Solver::getSolverHandler();

	// Get the step size
	hy = solverHandler.getStepSizeY();
	hz = solverHandler.getStepSizeZ();

	return solverHandler.getXGrid();
}
catch (const std::exception& e) {
	reportException(e);
	throw;
}

bool
XolotlInterface::getConvergenceStatus()
try {
	return solver->getConvergenceStatus();
}
catch (const std::exception& e) {
	reportException(e);
	throw;
}

void
XolotlInterface::finalizeXolotl()
try {
	// Call solver finalize
	solver->finalize();

	// Report statistics about the performance data collected during
	// the run we just completed.
	auto handlerRegistry = perf::PerfHandlerRegistry::get();
	perf::PerfObjStatsMap<perf::ITimer::ValType> timerStats;
	perf::PerfObjStatsMap<perf::IEventCounter::ValType> counterStats;
	perf::PerfObjStatsMap<perf::IHardwareCounter::CounterType> hwCtrStats;
	handlerRegistry->collectStatistics(timerStats, counterStats, hwCtrStats);

	auto xolotlComm = util::getMPIComm();

	// Get the MPI rank
	int rank;
	MPI_Comm_rank(xolotlComm, &rank);
	if (rank == 0) {
		handlerRegistry->reportStatistics(
			std::cout, timerStats, counterStats, hwCtrStats);
	}

	solver.reset();
}
catch (const std::exception& e) {
	reportException(e);
	throw;
}

} /* namespace interface */
} /* namespace xolotl */
