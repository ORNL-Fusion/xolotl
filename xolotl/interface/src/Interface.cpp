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
#include <xolotl/viz/VizHandlerRegistry.h>

namespace xolotl
{
namespace interface
{
class Context
{
public:
	Context(int argc, char* argv[]) : _kokkosContext(argc, argv)
	{
		if (!initialized()) {
			MPI_Init(&argc, &argv);
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

XolotlInterface::XolotlInterface() = default;

XolotlInterface::XolotlInterface(int argc, char* argv[], MPI_Comm mpiComm)
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
XolotlInterface::initializeXolotl(int argc, char* argv[], MPI_Comm comm)
{
	context = std::make_unique<Context>(argc, argv);

	// Initialize the MPI communicator to use
	util::setMPIComm(comm);
	auto xolotlComm = util::getMPIComm();

	// Get the MPI rank
	int rank;
	MPI_Comm_rank(xolotlComm, &rank);

	if (rank == 0) {
		// Print the start message
		std::cout << "Starting Xolotl Plasma-Surface Interactions Simulator"
				  << std::endl;
		// TODO! Print copyright message
		// Print date and time
		std::time_t currentTime = std::time(NULL);
		std::cout << std::asctime(std::localtime(&currentTime));
	}

	try {
		options::Options opts;
		opts.readParams(argc, argv);
		if (!opts.shouldRun()) {
			std::cerr << "Unable to read the options.  Aborting" << std::endl;
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
		std::cerr << e.what() << std::endl;
		std::cerr << "Aborting." << std::endl;
	}
	catch (const std::string& error) {
		std::cerr << error << std::endl;
		std::cerr << "Aborting." << std::endl;
	}
	catch (...) {
		std::cerr << "Unrecognized exception seen." << std::endl;
		std::cerr << "Aborting." << std::endl;
	}
}

void
XolotlInterface::setTimes(double finalTime, double dt)
{
	try {
		// Set the time in the solver
		solver->setTimes(finalTime, dt);
	}
	catch (const std::exception& e) {
		std::cerr << e.what() << std::endl;
		std::cerr << "Aborting." << std::endl;
	}
	catch (const std::string& error) {
		std::cerr << error << std::endl;
		std::cerr << "Aborting." << std::endl;
	}
	catch (...) {
		std::cerr << "Unrecognized exception seen." << std::endl;
		std::cerr << "Aborting." << std::endl;
	}
}

void
XolotlInterface::solveXolotl()
{
	try {
		// Launch the PetscSolver
		solver->solve();
	}
	catch (const std::exception& e) {
		std::cerr << e.what() << std::endl;
		std::cerr << "Aborting." << std::endl;
	}
	catch (const std::string& error) {
		std::cerr << error << std::endl;
		std::cerr << "Aborting." << std::endl;
	}
	catch (...) {
		std::cerr << "Unrecognized exception seen." << std::endl;
		std::cerr << "Aborting." << std::endl;
	}
}

std::vector<std::vector<std::vector<std::array<double, 4>>>>
XolotlInterface::getLocalNE()
{
	std::vector<std::vector<std::vector<std::array<double, 4>>>> toReturn;
	try {
		// Get the solver handler
		auto& solverHandler = solver::Solver::getSolverHandler();
		// Get the rate vector
		toReturn = solverHandler.getLocalNE();
	}
	catch (const std::exception& e) {
		std::cerr << e.what() << std::endl;
		std::cerr << "Aborting." << std::endl;
	}
	catch (const std::string& error) {
		std::cerr << error << std::endl;
		std::cerr << "Aborting." << std::endl;
	}
	catch (...) {
		std::cerr << "Unrecognized exception seen." << std::endl;
		std::cerr << "Aborting." << std::endl;
	}

	return toReturn;
}

void
XolotlInterface::setLocalNE(
	const std::vector<std::vector<std::vector<std::array<double, 4>>>>&
		rateVector)
{
	try {
		// Get the solver handler
		auto& solverHandler = solver::Solver::getSolverHandler();
		// Set the rate vector
		solverHandler.setLocalNE(rateVector);
	}
	catch (const std::exception& e) {
		std::cerr << e.what() << std::endl;
		std::cerr << "Aborting." << std::endl;
	}
	catch (const std::string& error) {
		std::cerr << error << std::endl;
		std::cerr << "Aborting." << std::endl;
	}
	catch (...) {
		std::cerr << "Unrecognized exception seen." << std::endl;
		std::cerr << "Aborting." << std::endl;
	}
}

void
XolotlInterface::getLocalCoordinates(int& xs, int& xm, int& Mx, int& ys,
	int& ym, int& My, int& zs, int& zm, int& Mz)
{
	try {
		// Get the solver handler
		auto& solverHandler = solver::Solver::getSolverHandler();
		// Get the local coordinates
		solverHandler.getLocalCoordinates(xs, xm, Mx, ys, ym, My, zs, zm, Mz);
	}
	catch (const std::exception& e) {
		std::cerr << e.what() << std::endl;
		std::cerr << "Aborting." << std::endl;
	}
	catch (const std::string& error) {
		std::cerr << error << std::endl;
		std::cerr << "Aborting." << std::endl;
	}
	catch (...) {
		std::cerr << "Unrecognized exception seen." << std::endl;
		std::cerr << "Aborting." << std::endl;
	}
}

void
XolotlInterface::setGBLocation(int i, int j, int k)
{
	try {
		// Get the solver handler
		auto& solverHandler = solver::Solver::getSolverHandler();
		// Set the coordinate of the GB
		solverHandler.setGBLocation(i, j, k);
	}
	catch (const std::exception& e) {
		std::cerr << e.what() << std::endl;
		std::cerr << "Aborting." << std::endl;
	}
	catch (const std::string& error) {
		std::cerr << error << std::endl;
		std::cerr << "Aborting." << std::endl;
	}
	catch (...) {
		std::cerr << "Unrecognized exception seen." << std::endl;
		std::cerr << "Aborting." << std::endl;
	}
}

void
XolotlInterface::resetGBVector()
{
	try {
		// Get the solver handler
		auto& solverHandler = solver::Solver::getSolverHandler();
		// Reset the location
		solverHandler.resetGBVector();
	}
	catch (const std::exception& e) {
		std::cerr << e.what() << std::endl;
		std::cerr << "Aborting." << std::endl;
	}
	catch (const std::string& error) {
		std::cerr << error << std::endl;
		std::cerr << "Aborting." << std::endl;
	}
	catch (...) {
		std::cerr << "Unrecognized exception seen." << std::endl;
		std::cerr << "Aborting." << std::endl;
	}
}

std::vector<std::vector<std::vector<std::vector<std::pair<int, double>>>>>
XolotlInterface::getConcVector()
{
	std::vector<std::vector<std::vector<std::vector<std::pair<int, double>>>>>
		toReturn;
	try {
		// Get the vector
		toReturn = solver->getConcVector();
	}
	catch (const std::exception& e) {
		std::cerr << e.what() << std::endl;
		std::cerr << "Aborting." << std::endl;
	}
	catch (const std::string& error) {
		std::cerr << error << std::endl;
		std::cerr << "Aborting." << std::endl;
	}
	catch (...) {
		std::cerr << "Unrecognized exception seen." << std::endl;
		std::cerr << "Aborting." << std::endl;
	}

	return toReturn;
}

void
XolotlInterface::setConcVector(
	std::vector<std::vector<std::vector<std::vector<std::pair<int, double>>>>>
		concVector)
{
	try {
		// Set the vector
		solver->setConcVector(concVector);
	}
	catch (const std::exception& e) {
		std::cerr << e.what() << std::endl;
		std::cerr << "Aborting." << std::endl;
	}
	catch (const std::string& error) {
		std::cerr << error << std::endl;
		std::cerr << "Aborting." << std::endl;
	}
	catch (...) {
		std::cerr << "Unrecognized exception seen." << std::endl;
		std::cerr << "Aborting." << std::endl;
	}
}

double
XolotlInterface::getPreviousTime()
{
	double toReturn = 0.0;
	try {
		// Get the solver handler
		auto& solverHandler = solver::Solver::getSolverHandler();
		toReturn = solverHandler.getPreviousTime();
	}
	catch (const std::exception& e) {
		std::cerr << e.what() << std::endl;
		std::cerr << "Aborting." << std::endl;
	}
	catch (const std::string& error) {
		std::cerr << error << std::endl;
		std::cerr << "Aborting." << std::endl;
	}
	catch (...) {
		std::cerr << "Unrecognized exception seen." << std::endl;
		std::cerr << "Aborting." << std::endl;
	}

	return toReturn;
}

void
XolotlInterface::setPreviousTime(double time)
{
	try {
		// Get the solver handler
		auto& solverHandler = solver::Solver::getSolverHandler();
		solverHandler.setPreviousTime(
			time, true); // Update the fluence from here
	}
	catch (const std::exception& e) {
		std::cerr << e.what() << std::endl;
		std::cerr << "Aborting." << std::endl;
	}
	catch (const std::string& error) {
		std::cerr << error << std::endl;
		std::cerr << "Aborting." << std::endl;
	}
	catch (...) {
		std::cerr << "Unrecognized exception seen." << std::endl;
		std::cerr << "Aborting." << std::endl;
	}
}

double
XolotlInterface::getCurrentDt()
{
	double toReturn = 0.0;
	try {
		toReturn = solver->getCurrentDt();
	}
	catch (const std::exception& e) {
		std::cerr << e.what() << std::endl;
		std::cerr << "Aborting." << std::endl;
	}
	catch (const std::string& error) {
		std::cerr << error << std::endl;
		std::cerr << "Aborting." << std::endl;
	}
	catch (...) {
		std::cerr << "Unrecognized exception seen." << std::endl;
		std::cerr << "Aborting." << std::endl;
	}

	return toReturn;
}

void
XolotlInterface::setCurrentTimes(double time, double dt)
{
	try {
		solver->setCurrentTimes(time, dt);
	}
	catch (const std::exception& e) {
		std::cerr << e.what() << std::endl;
		std::cerr << "Aborting." << std::endl;
	}
	catch (const std::string& error) {
		std::cerr << error << std::endl;
		std::cerr << "Aborting." << std::endl;
	}
	catch (...) {
		std::cerr << "Unrecognized exception seen." << std::endl;
		std::cerr << "Aborting." << std::endl;
	}
}

double
XolotlInterface::getNXeGB()
{
	double toReturn = 0.0;
	try {
		// Get the solver handler
		auto& solverHandler = solver::Solver::getSolverHandler();
		toReturn = solverHandler.getNXeGB();
	}
	catch (const std::exception& e) {
		std::cerr << e.what() << std::endl;
		std::cerr << "Aborting." << std::endl;
	}
	catch (const std::string& error) {
		std::cerr << error << std::endl;
		std::cerr << "Aborting." << std::endl;
	}
	catch (...) {
		std::cerr << "Unrecognized exception seen." << std::endl;
		std::cerr << "Aborting." << std::endl;
	}

	return toReturn;
}

void
XolotlInterface::setNXeGB(double nXe)
{
	try {
		// Get the solver handler
		auto& solverHandler = solver::Solver::getSolverHandler();
		solverHandler.setNXeGB(nXe);
	}
	catch (const std::exception& e) {
		std::cerr << e.what() << std::endl;
		std::cerr << "Aborting." << std::endl;
	}
	catch (const std::string& error) {
		std::cerr << error << std::endl;
		std::cerr << "Aborting." << std::endl;
	}
	catch (...) {
		std::cerr << "Unrecognized exception seen." << std::endl;
		std::cerr << "Aborting." << std::endl;
	}
}

TS&
XolotlInterface::getTS()
{
	return solver->getTS();
}

std::vector<double>
XolotlInterface::getGridInfo(double& hy, double& hz)
{
	// The vector to return
	std::vector<double> toReturn;
	try {
		// Get the solver handler
		auto& solverHandler = solver::Solver::getSolverHandler();
		// Get the grid
		toReturn = solverHandler.getXGrid();
		// Get the step size
		hy = solverHandler.getStepSizeY();
		hz = solverHandler.getStepSizeZ();
	}
	catch (const std::exception& e) {
		std::cerr << e.what() << std::endl;
		std::cerr << "Aborting." << std::endl;
	}
	catch (const std::string& error) {
		std::cerr << error << std::endl;
		std::cerr << "Aborting." << std::endl;
	}
	catch (...) {
		std::cerr << "Unrecognized exception seen." << std::endl;
		std::cerr << "Aborting." << std::endl;
	}

	return toReturn;
}

bool
XolotlInterface::getConvergenceStatus()
{
	bool toReturn = true;
	try {
		toReturn = solver->getConvergenceStatus();
	}
	catch (const std::exception& e) {
		std::cerr << e.what() << std::endl;
		std::cerr << "Aborting." << std::endl;
	}
	catch (const std::string& error) {
		std::cerr << error << std::endl;
		std::cerr << "Aborting." << std::endl;
	}
	catch (...) {
		std::cerr << "Unrecognized exception seen." << std::endl;
		std::cerr << "Aborting." << std::endl;
	}

	return toReturn;
}

void
XolotlInterface::finalizeXolotl()
{
	try {
		// Call solver finalize
		solver->finalize();

		// Report statistics about the performance data collected during
		// the run we just completed.
		auto handlerRegistry = perf::PerfHandlerRegistry::get();
		perf::PerfObjStatsMap<perf::ITimer::ValType> timerStats;
		perf::PerfObjStatsMap<perf::IEventCounter::ValType> counterStats;
		perf::PerfObjStatsMap<perf::IHardwareCounter::CounterType> hwCtrStats;
		handlerRegistry->collectStatistics(
			timerStats, counterStats, hwCtrStats);

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
		std::cerr << e.what() << std::endl;
		std::cerr << "Aborting." << std::endl;
	}
	catch (const std::string& error) {
		std::cerr << error << std::endl;
		std::cerr << "Aborting." << std::endl;
	}
	catch (...) {
		std::cerr << "Unrecognized exception seen." << std::endl;
		std::cerr << "Aborting." << std::endl;
	}
}

} /* namespace interface */
} /* namespace xolotl */
