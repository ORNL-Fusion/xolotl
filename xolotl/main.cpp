/**
 * Main.c, currently only able to load clusters
 */
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cassert>
#include <Reactant.h>
#include <PSIClusterNetworkLoader.h>
#include <PetscSolver.h>
#include <mpi.h>
#include <MPIUtils.h>
#include <XolotlOptions.h>
#include "xolotlPerf.h"
#include <MaterialHandlerFactory.h>
#include <TemperatureHandlerFactory.h>
#include <VizHandlerRegistryFactory.h>
#include <HDF5NetworkLoader.h>
#include <IVizHandlerRegistry.h>

using namespace std;
using std::shared_ptr;
namespace xperf = xolotlPerf;


//! This operation prints the start message
void printStartMessage() {
	cout << "Starting Xolotl Plasma-Surface Interactions Simulator" << endl;
	// TODO! Print copyright message
	// TODO! Print date and time
}

bool initMaterial() {

	bool materialInitOK = xolotlSolver::initializeMaterial();
	if (!materialInitOK) {
		std::cerr << "Unable to initialize requested material.  Aborting"
				<< std::endl;
		return EXIT_FAILURE;
	} else
		return materialInitOK;
}

bool initTemp(bool opts, bool opts1, XolotlOptions &options) {

	bool tempInitOK = xolotlSolver::initializeTempHandler(opts, opts1, options);
	if (!tempInitOK) {
		std::cerr << "Unable to initialize requested temperature.  Aborting"
				<< std::endl;
		return EXIT_FAILURE;
	} else
		return tempInitOK;
}


bool initViz(bool opts) {

	bool vizInitOK = xolotlViz::initialize(opts);
	if (!vizInitOK) {
		std::cerr
				<< "Unable to initialize requested visualization infrastructure. "
				<< "Aborting" << std::endl;
		return EXIT_FAILURE;
	} else
		return vizInitOK;
}

std::shared_ptr<xolotlSolver::PetscSolver>
setUpSolver( std::shared_ptr<xolotlPerf::IHandlerRegistry> handlerRegistry, 
            int argc, char **argv) {
	// Setup the solver
	//auto solverInitTimer = handlerRegistry->getTimer("initSolver");
	//solverInitTimer->start();
	std::shared_ptr<xolotlSolver::PetscSolver> solver = 
        std::make_shared<xolotlSolver::PetscSolver>(handlerRegistry);
	solver->setCommandLineOptions(argc, argv);
	solver->initialize();
	//solverInitTimer->stop();

	return solver;
}

void launchPetscSolver(std::shared_ptr<xolotlSolver::PetscSolver> solver,
		std::shared_ptr<xolotlPerf::IHandlerRegistry> handlerRegistry,
		std::shared_ptr<xolotlSolver::IFluxHandler> materialHandler,
		std::shared_ptr<xolotlSolver::ITemperatureHandler> tempHandler) {

    xperf::IHardwareCounter::SpecType hwctrSpec;
    hwctrSpec.push_back( xperf::IHardwareCounter::FPOps );
    hwctrSpec.push_back( xperf::IHardwareCounter::Cycles );
    hwctrSpec.push_back( xperf::IHardwareCounter::L3CacheMisses );

	// Launch the PetscSolver
	auto solverTimer = handlerRegistry->getTimer("solve");
    auto solverHwctr = handlerRegistry->getHardwareCounter( "solve", hwctrSpec );
	solverTimer->start();
    solverHwctr->start();
	solver->solve(materialHandler, tempHandler);
    solverHwctr->stop();
	solverTimer->stop();
}

std::shared_ptr<PSIClusterNetworkLoader> setUpNetworkLoader(int rank,
		MPI_Comm comm, std::string networkFilename,
		std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) {

	// Create a HDF5NetworkLoader
	std::shared_ptr<HDF5NetworkLoader> networkLoader;
	networkLoader = std::make_shared<HDF5NetworkLoader>(registry);
	// Give the networkFilename to the network loader
	networkLoader->setFilename(networkFilename);

	return networkLoader;
}


//! Main program
int main(int argc, char **argv) {

	// Local Declarations
	int rank;

	// Check the command line arguments.
	// Skip the executable name before parsing.
	argc -= 1; // one for the executable name
	argv += 1; // one for the executable name

	XolotlOptions xopts;
	xopts.readParams(argc, argv);
	if (!xopts.shouldRun()) {
		return xopts.getExitCode();
	}

	// Skip the name of the parameter file that was just used.
	// The arguments should be empty now.
	argc -= 1;
	argv += 1;

	// Extract the argument for the file name
	std::string networkFilename = xopts.getNetworkFilename();
	assert(!networkFilename.empty());

	try {
		auto materialInitOK = initMaterial();
		auto tempInitOK = initTemp(xopts.useConstTemperatureHandlers(),
				xopts.useTemperatureProfileHandlers(), xopts);

		// Set up our performance data infrastructure.
        xperf::initialize(xopts.getPerfRegistryType());

		// Set up the visualization infrastructure.
		auto vizInitOK = initViz(xopts.useVizStandardHandlers());

		// Initialize MPI. We do this instead of leaving it to some
		// other package (e.g., PETSc), because we want to avoid problems
		// with overlapping Timer scopes.
		MPI_Init(&argc, &argv);

		auto materialHandler = xolotlSolver::getMaterialHandler();
		auto tempHandler = xolotlSolver::getTemperatureHandler(xopts);

		// Access our handler registry to obtain a Timer
		// measuring the runtime of the entire program.
		// NOTE: these long template types could be replaced with 'auto'
		auto handlerRegistry = xolotlPerf::getHandlerRegistry();
		auto totalTimer = handlerRegistry->getTimer("total");
		totalTimer->start();

		// Setup the solver
		auto solver = setUpSolver(handlerRegistry,
				xopts.getPetscArgc(), xopts.getPetscArgv());

		// Load the network
		auto networkLoadTimer = handlerRegistry->getTimer("loadNetwork");
		networkLoadTimer->start();

		// Get the MPI rank
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);

		// Set up the network loader
		auto networkLoader = setUpNetworkLoader(rank, MPI_COMM_WORLD,
				networkFilename, handlerRegistry);

		// Give the network loader to PETSc as input
		solver->setNetworkLoader(networkLoader);
		networkLoadTimer->stop();

		// Launch the PetscSolver
		launchPetscSolver(solver, handlerRegistry, materialHandler,
				tempHandler);

		// Finalize our use of the solver.
		solver->finalize();

		totalTimer->stop();

        // Report the performance data about the run we just completed.
        handlerRegistry->reportStatistics(std::cout);

    } catch (std::exception& e ) {

        std::cerr << e.what() << std::endl;
        std::cerr << "Aborting." << std::endl;
        return EXIT_FAILURE;

	} catch (std::string & error) {
		std::cout << error << std::endl;
		std::cout << "Aborting." << std::endl;
		return EXIT_FAILURE;
	}

	// finalize our use of MPI
	MPI_Finalize();

	return EXIT_SUCCESS;
}
