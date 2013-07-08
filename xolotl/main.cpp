/**
 * Main.c, currently only able to load clusters
 */
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <Reactant.h>
#include <PSIClusterNetworkLoader.h>
#include <PetscSolver.h>
#include <mpi.h>
#include "xolotlCore/io/MPIUtils.h"

using namespace std;
using std::shared_ptr;


//! This operation prints the start message
void printStartMessage() {
	cout << "Starting Xolotl Plasma-Surface Interactions Simulator" << endl;
	// TODO! Print copyright message
	// TODO! Print date and time
}

//! This operation prints proper usage instructions for Xolotl
void printUsage() {
	cout << "Usage:" << endl;
	cout << "\txolotl <network file name>" << endl;
}

//! Main program
int main(int argc, char **argv) {

	printStartMessage();
	
	// Check the arguments
	if (argc < 2) {
		cout << "Insufficient input provided! Aborting!" << std::endl;
		printUsage();
		return EXIT_FAILURE;
	}
	
	// Extract the argument values
	const char *networkFilename = argv[1];
	
	try {
		// Setup the solver
		xolotlSolver::PetscSolver solver;
		solver.setCommandLineOptions(argc, argv);
		solver.initialize();
		
		// Load the input file from the master task
		shared_ptr<std::istream> networkStream;
		
		int rank;
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		
		if (rank == 0) {
			networkStream.reset(new std::ifstream(networkFilename));
		}
		
		// Broadcast the stream to all worker tasks
		networkStream = xolotlCore::MPIUtils::broadcastStream(
			networkStream, 0, MPI_COMM_WORLD);
		
		// Create a network loader and set the istream on every MPI task
		std::shared_ptr<PSIClusterNetworkLoader> networkLoader(
			new PSIClusterNetworkLoader());
		networkLoader->setInputstream(networkStream);
		
		solver.setNetworkLoader(networkLoader);
		solver.setupSolver();
		
		// Launch the PetscSolver
		/* solver.solve(); */
		solver.finalize();
	}
	catch (std::string error) {
		std::cout << error << std::endl;
		std::cout << "Aborting." << std::endl;
		return EXIT_FAILURE;
	}
	
	return EXIT_SUCCESS;
}
