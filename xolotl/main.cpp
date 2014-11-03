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
#include <Options.h>
#include <xolotlPerf.h>
#include <IMaterialFactory.h>
#include <TemperatureHandlerFactory.h>
#include <VizHandlerRegistryFactory.h>
#include <HDF5NetworkLoader.h>
#include <ctime>

using namespace std;
using std::shared_ptr;
namespace xperf = xolotlPerf;


//! This operation prints the start message
void printStartMessage() {
	std::cout << "Starting Xolotl Plasma-Surface Interactions Simulator" << std::endl;
	// TODO! Print copyright message
	// Print date and time
	std::time_t currentTime = std::time(NULL);
	std::cout << std::asctime(std::localtime(&currentTime)); // << std::endl;
}

std::shared_ptr<xolotlFactory::IMaterialFactory> initMaterial(Options &options) {
	// Create the material factory
	auto materialFactory = xolotlFactory::IMaterialFactory::createMaterialFactory(options.getMaterial());

	// Initialize it with the options
	materialFactory->initializeMaterial(options);

	return materialFactory;
}

bool initTemp(Options &options) {

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

std::shared_ptr<xolotlSolver::PetscSolver> setUpSolver(
		std::shared_ptr<xolotlPerf::IHandlerRegistry> handlerRegistry, 
        const char* argv0,
        int argc, char **argv) {
	// Setup the solver
	auto solverInitTimer = handlerRegistry->getTimer("initSolver");
	solverInitTimer->start();
	std::shared_ptr<xolotlSolver::PetscSolver> solver = std::make_shared<
			xolotlSolver::PetscSolver>(handlerRegistry);

    // PETSc assumes that argv[0] in the arguments it is given is the
    // program name.  But our parsing of the PETSc arguments from
    // the input parameter file gives us only the PETSc arguments without
    // the program name as argv[0].  So - we adjust the arguments array
    // so that it has the right argv[0].
    int petscArgc = argc + 1;
    char** petscArgv = new char*[petscArgc+1];
    petscArgv[0] = new char[strlen(argv0)+1];
    strcpy( petscArgv[0], argv0 );
    for( int idx = 0; idx < petscArgc; ++idx )
    {
        petscArgv[idx+1] = argv[idx];
    }
    petscArgv[petscArgc] = NULL;
	solver->setCommandLineOptions(petscArgc, petscArgv);

	solver->initialize();
	solverInitTimer->stop();

	return solver;
}

void launchPetscSolver(std::shared_ptr<xolotlSolver::PetscSolver> solver,
		std::shared_ptr<xolotlPerf::IHandlerRegistry> handlerRegistry,
		std::shared_ptr<xolotlFactory::IMaterialFactory> material,
		std::shared_ptr<xolotlCore::ITemperatureHandler> tempHandler,
		double stepSize) {

    xperf::IHardwareCounter::SpecType hwctrSpec;
    hwctrSpec.push_back( xperf::IHardwareCounter::FPOps );
    hwctrSpec.push_back( xperf::IHardwareCounter::Cycles );
    hwctrSpec.push_back( xperf::IHardwareCounter::L3CacheMisses );

	// Launch the PetscSolver
	auto solverTimer = handlerRegistry->getTimer("solve");
    auto solverHwctr = handlerRegistry->getHardwareCounter( "solve", hwctrSpec );
	solverTimer->start();
    solverHwctr->start();
	solver->solve(material, tempHandler, stepSize);
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
	// Skip the executable name before parsing (but save it, 
    // because PETSc expects it in its argument list).
    char* progName = argv[0];
	argc -= 1; // one for the executable name
	argv += 1; // one for the executable name
	Options opts;
	opts.readParams(argc, argv);
	if (!opts.shouldRun()) {
		return opts.getExitCode();
	}

	// Skip the name of the parameter file that was just used.
	// The arguments should be empty now.
	argc -= 1;
	argv += 1;

	// Extract the argument for the file name
	std::string networkFilename = opts.getNetworkFilename();
	assert(!networkFilename.empty());

	try {
		// Set up our performance data infrastructure.
        xperf::initialize(opts.getPerfHandlerType());

		// Initialize MPI. We do this instead of leaving it to some
		// other package (e.g., PETSc), because we want to avoid problems
		// with overlapping Timer scopes.
		MPI_Init(&argc, &argv);

		// Get the MPI rank
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);

		if (rank == 0) {
			// Print the start message
			printStartMessage();
		}

		// Set up the material infrastructure that is used to calculate flux
		auto material = initMaterial(opts);
		// Set up the temperature infrastructure
		auto tempInitOK = initTemp(opts);
		// Set up the visualization infrastructure.
		auto vizInitOK = initViz(opts.useVizStandardHandlers());

		// Access the temperature handler registry to get the temperature
		auto tempHandler = xolotlFactory::getTemperatureHandler();

		// Access our performance handler registry to obtain a Timer
		// measuring the runtime of the entire program.
		// NOTE: these long template types could be replaced with 'auto'
		auto handlerRegistry = xolotlPerf::getHandlerRegistry();
		auto totalTimer = handlerRegistry->getTimer("total");
		totalTimer->start();

		// Setup the solver
		auto solver = setUpSolver(handlerRegistry, progName,
                                    opts.getPetscArgc(), opts.getPetscArgv());

		// Load the network
		auto networkLoadTimer = handlerRegistry->getTimer("loadNetwork");
		networkLoadTimer->start();

		// Set up the network loader
		auto networkLoader = setUpNetworkLoader(rank, MPI_COMM_WORLD,
				networkFilename, handlerRegistry);

		// Give the network loader to PETSc as input
		solver->setNetworkLoader(networkLoader);
		networkLoadTimer->stop();

		// Launch the PetscSolver
		launchPetscSolver(solver, handlerRegistry, material,
				tempHandler, opts.getStepSize());

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
        handlerRegistry->collectStatistics( timerStats, counterStats, hwCtrStats );
        if( rank == 0 )
        {
            handlerRegistry->reportStatistics( std::cout, 
                                                timerStats, 
                                                counterStats, 
                                                hwCtrStats );
        }

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
