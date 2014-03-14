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
#include "xolotlCore/io/MPIUtils.h"
#include "xolotlCore/commandline/XolotlOptions.h"
#include "xolotlPerf/HandlerRegistryFactory.h"
#include "xolotlPerf/HardwareQuantities.h"


using namespace std;
using std::shared_ptr;

//! This operation prints the start message
void printStartMessage() {
	cout << "Starting Xolotl Plasma-Surface Interactions Simulator" << endl;
	// TODO! Print copyright message
	// TODO! Print date and time
}


//! Main program
int main(int argc, char **argv) {

	// Local Declarations
	shared_ptr<std::istream> networkStream;
	std::shared_ptr<PSIClusterNetworkLoader> networkLoader;
	int rank;

	// Check the command line arguments.
    // Skip the executable name before parsing.
    argc -= 1;  // one for the executable name
    argv += 1;  // one for the executable name
    XolotlOptions xopts;
    int nOptsUsed = xopts.parseCommandLine( argc, argv );
    if( !xopts.shouldRun() )
    {
        return xopts.getExitCode();
    }
    argc -= nOptsUsed;
    argv += nOptsUsed;

	// Extract the argument for the file name
    std::string networkFilename = xopts.getNetworkFilename();
    assert( !networkFilename.empty() );

	try {
        // Set up our performance data infrastructure.
        // Indicate we want to monitor some important hardware counters.
        std::vector<xolotlPerf::HardwareQuantities> hwq;
        hwq.push_back( xolotlPerf::FP_OPS );
        hwq.push_back( xolotlPerf::L3_CACHE_MISS );
        bool perfInitOK = xolotlPerf::initialize( xopts.useStandardHandlers(), hwq );
        if( !perfInitOK )
        {
            std::cerr << "Unable to initialize requested performance data infrastructure.  Aborting" << std::endl;
            return EXIT_FAILURE;
        }

        // Initialize MPI.  We do this instead of leaving it to some 
        // other package (e.g., PETSc), because we want to avoid problems 
        // with overlapping Timer scopes.
        MPI_Init( &argc, &argv );

        // Access our handler registry to obtain a Timer 
        // measuring the runtime of the entire program.
        // NOTE: these long template types could be replaced with 'auto'
        auto handlerRegistry = xolotlPerf::getHandlerRegistry();
        auto totalTimer = handlerRegistry->getTimer( "total" );
        totalTimer->start();

		// Setup the solver
        auto solverInitTimer = handlerRegistry->getTimer( "initSolver" );
        solverInitTimer->start();
		xolotlSolver::PetscSolver solver;
		solver.setCommandLineOptions(argc, argv);
		solver.initialize();
        solverInitTimer->stop();

        auto networkLoadTimer =  handlerRegistry->getTimer( "loadNetwork" );
        networkLoadTimer->start();

		// Get the MPI rank
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		// Setup the master
		if (rank == 0) {
			// Say hello
			printStartMessage();
			// Set the input stream on the master
			networkStream = make_shared<std::ifstream>(networkFilename);
		}

		// Broadcast the stream to all worker tasks
		networkLoader = std::make_shared<
					PSIClusterNetworkLoader>();
		networkStream = xolotlCore::MPIUtils::broadcastStream(networkStream, 0,
				MPI_COMM_WORLD );

		// Create a network loader and set the stream on every MPI task
		networkLoader->setInputstream(networkStream);
		// Give the network loader to PETSc as input
		solver.setNetworkLoader(networkLoader);
        networkLoadTimer->stop();

		// Launch the PetscSolver
        auto solverTimer = handlerRegistry->getTimer( "solve" );
        solverTimer->start();
		solver.solve();
        solverTimer->stop();

        // Finalize our use of the solver.
		solver.finalize();
        totalTimer->stop();

        // Report the performance data about the run we just completed
        // TODO Currently, this call writes EventCounter data to the
        // given stream, but Timer and any hardware counter data is
        // written by the underlying timing library to files, one per process.
        if( rank == 0 ) {
            handlerRegistry->dump( std::cout );
        }

	} catch (std::string & error) {
		std::cout << error << std::endl;
		std::cout << "Aborting." << std::endl;
		return EXIT_FAILURE;
	}

    // finalize our use of MPI
    MPI_Finalize();

	return EXIT_SUCCESS;
}
