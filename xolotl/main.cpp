/**
 * Main.c, currently only able to load clusters
 */
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cassert>
#include <Reactant.h>
#include <PetscSolver.h>
#include <mpi.h>
#include <MPIUtils.h>
#include <Options.h>
#include <xolotlPerf.h>
#include <IMaterialFactory.h>
#include <TemperatureHandlerFactory.h>
#include <VizHandlerRegistryFactory.h>
#include <INetworkLoader.h>
#include <IReactionNetwork.h>
#include <SolverHandlerFactory.h>
#include <ISolverHandler.h>
#include <IReactionHandlerFactory.h>
#include <ctime>

using namespace std;
using std::shared_ptr;
namespace xperf = xolotlPerf;

//! This operation prints the start message
void printStartMessage() {
	std::cout << "Starting Xolotl Plasma-Surface Interactions Simulator"
			<< std::endl;
	// TODO! Print copyright message
	// Print date and time
	std::time_t currentTime = std::time(NULL);
	std::cout << std::asctime(std::localtime(&currentTime)); // << std::endl;
}

std::shared_ptr<xolotlFactory::IMaterialFactory> initMaterial(
		const Options &options) {
	// Create the material factory
	auto materialFactory =
			xolotlFactory::IMaterialFactory::createMaterialFactory(options);

	// Initialize it with the options
	materialFactory->initializeMaterial(options);

	return materialFactory;
}

bool initTemp(const Options &options) {

	bool tempInitOK = xolotlFactory::initializeTempHandler(options);
	if (!tempInitOK) {
		std::cerr << "Unable to initialize requested temperature.  Aborting"
				<< std::endl;
		return EXIT_FAILURE;
	} else
		return tempInitOK;
}

bool initViz(bool opts) {

	bool vizInitOK = xolotlFactory::initializeVizHandler(opts);
	if (!vizInitOK) {
		std::cerr
				<< "Unable to initialize requested visualization infrastructure. "
				<< "Aborting" << std::endl;
		return EXIT_FAILURE;
	} else
		return vizInitOK;
}

std::unique_ptr<xolotlSolver::PetscSolver> setUpSolver(
		std::shared_ptr<xolotlPerf::IHandlerRegistry> handlerRegistry,
		std::shared_ptr<xolotlFactory::IMaterialFactory> material,
		std::shared_ptr<xolotlCore::ITemperatureHandler> tempHandler,
		xolotlSolver::ISolverHandler &solvHandler, const Options &options) {
	// Initialize the solver handler
	solvHandler.initializeHandlers(material, tempHandler, options);

	// Setup the solver
	auto solverInitTimer = handlerRegistry->getTimer("initSolver");
	solverInitTimer->start();
	// Once we have widespread C++14 support, use std::make_unique.
	std::unique_ptr<xolotlSolver::PetscSolver> solver(
			new xolotlSolver::PetscSolver(solvHandler, handlerRegistry));
	solver->setCommandLineOptions(options.getPetscArg());
	solver->initialize();
	solverInitTimer->stop();

	return solver;
}

void launchPetscSolver(xolotlSolver::PetscSolver &solver,
		std::shared_ptr<xolotlPerf::IHandlerRegistry> handlerRegistry) {

	xperf::IHardwareCounter::SpecType hwctrSpec;
	hwctrSpec.push_back(xperf::IHardwareCounter::FPOps);
	hwctrSpec.push_back(xperf::IHardwareCounter::Cycles);
	hwctrSpec.push_back(xperf::IHardwareCounter::L3CacheMisses);

	// Launch the PetscSolver
	auto solverTimer = handlerRegistry->getTimer("solve");
	auto solverHwctr = handlerRegistry->getHardwareCounter("solve", hwctrSpec);
	solverTimer->start();
	solverHwctr->start();
	solver.solve();
	solverHwctr->stop();
	solverTimer->stop();
}

//! Run the Xolotl simulation.
int runXolotl(const Options &opts) {

	// Set up our performance data infrastructure.
	xperf::initialize(opts.getPerfHandlerType());

	// Get the MPI rank
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if (rank == 0) {
		// Print the start message
		printStartMessage();
	}

	// Set up the material infrastructure that is used to calculate flux
	auto material = initMaterial(opts);
	// Set up the temperature infrastructure
	bool tempInitOK = initTemp(opts);
	if (!tempInitOK) {
		throw std::runtime_error("Unable to initialize temperature.");
	}
	// Set up the visualization infrastructure.
	bool vizInitOK = initViz(opts.useVizStandardHandlers());
	if (!vizInitOK) {
		throw std::runtime_error(
				"Unable to initialize visualization infrastructure.");
	}

	// Access the temperature handler registry to get the temperature
	auto tempHandler = xolotlFactory::getTemperatureHandler();

	// Access our performance handler registry to obtain a Timer
	// measuring the runtime of the entire program.
	auto handlerRegistry = xolotlPerf::getHandlerRegistry();
	auto totalTimer = handlerRegistry->getTimer("total");
	totalTimer->start();

	// Create the network handler factory
	auto networkFactory =
			xolotlFactory::IReactionHandlerFactory::createNetworkFactory(
					opts.getMaterial());

	// Build a reaction network
	auto networkLoadTimer = handlerRegistry->getTimer("loadNetwork");
	networkLoadTimer->start();
	networkFactory->initializeReactionNetwork(opts, handlerRegistry);
	networkLoadTimer->stop();
	if (rank == 0) {
		std::time_t currentTime = std::time(NULL);
		std::cout << std::asctime(std::localtime(&currentTime));
	}
	auto &network = networkFactory->getNetworkHandler();

	// Initialize and get the solver handler
	bool dimOK = xolotlFactory::initializeDimension(opts, network);
	if (!dimOK) {
		throw std::runtime_error("Unable to initialize dimension from inputs.");
	}
	auto &solvHandler = xolotlFactory::getSolverHandler();

	// Setup the solver
	auto solver = setUpSolver(handlerRegistry, material, tempHandler,
			solvHandler, opts);

	// Launch the PetscSolver
	launchPetscSolver(*solver, handlerRegistry);

	// Finalize our use of the solver.
	auto solverFinalizeTimer = handlerRegistry->getTimer("solverFinalize");
	solverFinalizeTimer->start();
	solver->finalize();
	solverFinalizeTimer->stop();

	totalTimer->stop();

	// Report statistics about the performance data collected during
	// the run we just completed.
	xperf::PerfObjStatsMap<xperf::ITimer::ValType> timerStats;
	xperf::PerfObjStatsMap<xperf::IEventCounter::ValType> counterStats;
	xperf::PerfObjStatsMap<xperf::IHardwareCounter::CounterType> hwCtrStats;
	handlerRegistry->collectStatistics(timerStats, counterStats, hwCtrStats);
	if (rank == 0) {
		handlerRegistry->reportStatistics(std::cout, timerStats, counterStats,
				hwCtrStats);
	}

	return 0;
}

//! Main program
int main(int argc, char **argv) {

	// Local Declarations
	int ret = EXIT_SUCCESS;

	// Initialize MPI. We do this instead of leaving it to some
	// other package (e.g., PETSc), because we want to avoid problems
	// with overlapping Timer scopes.
	// We do this before our own parsing of the command line,
	// because it may change the command line.
	MPI_Init(&argc, &argv);

	try {
		// Check the command line arguments.
		Options opts;
		opts.readParams(argc, argv);
		if (opts.shouldRun()) {
			// Skip the name of the parameter file that was just used.
			// The arguments should be empty now.
			// TODO is this needed?
			argc -= 2;
			argv += 2;

			// Run the simulation.
			ret = runXolotl(opts);
		} else {
			ret = opts.getExitCode();
		}
	} catch (const std::exception &e) {
		std::cerr << e.what() << std::endl;
		std::cerr << "Aborting." << std::endl;
		ret = EXIT_FAILURE;
	} catch (const std::string &error) {
		std::cerr << error << std::endl;
		std::cerr << "Aborting." << std::endl;
		ret = EXIT_FAILURE;
	} catch (...) {
		std::cerr << "Unrecognized exception seen." << std::endl;
		std::cerr << "Aborting." << std::endl;
		ret = EXIT_FAILURE;
	}

	// Clean up.
	MPI_Finalize();

	return ret;
}
